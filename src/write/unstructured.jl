# UnstructuredGrid writer.

mutable struct UnstructuredState <: DatasetKind
    total_points::Int
    total_cells::Int
    total_conn::Int
    total_parts::Int
    # polyhedron bookkeeping
    has_polyhedra::Bool
    total_faces::Int
    total_face_conn::Int
    total_p2f::Int
    part_cells::Vector{Int}  # cells per partition, for lazy polyhedron backfill
    # current step
    step_parts::Int
    step_points::Int
    step_cells::Int
    snap::NTuple{7, Int}  # (parts, points, cells, conn, faces, face_conn, p2f) at step begin
    # geometry of the most recently written step (reused when a step adds none)
    geom::Union{Nothing, NamedTuple}
    # cumulative expected data rows over completed steps
    cum_points::Int
    cum_cells::Int
end

new_unstructured_state() = UnstructuredState(
    0, 0, 0, 0,
    false, 0, 0, 0, Int[],
    0, 0, 0, (0, 0, 0, 0, 0, 0, 0), nothing, 0, 0
)

const UnstructuredCells = AbstractVector{<:Union{MeshCell{VTKCellType}, VTKPolyhedron}}

function init_unstructured(
        dest, points, cells; temporal = false,
        compress = false, chunk_size = 0
    )
    file, root = open_dest(dest)
    write_ascii_attribute(root, "Type", "UnstructuredGrid")
    write_version_attribute(root, (2, 0))
    vtk = make_vtkfile(
        file, root, new_unstructured_state();
        temporal, compress, chunk_size, version = (2, 0)
    )
    if points !== nothing
        add_partition(vtk, points, cells)
    end
    return vtk
end

totals_snapshot(k::UnstructuredState) =
    (
    k.total_parts, k.total_points, k.total_cells, k.total_conn,
    k.total_faces, k.total_face_conn, k.total_p2f,
)

"""
    add_partition(vtk, points, cells; pointdata = (), celldata = ())        # UnstructuredGrid
    add_partition(vtk, points, cellvecs...; pointdata = (), celldata = ())  # PolyData

Append one partition of geometry (and optionally its data) to an
unstructured grid or poly data file.

Such a file is usually created without geometry, via
`vtkhdf_grid(VTKUnstructuredGrid(), filename)` /
`vtkhdf_grid(VTKPolyData(), filename)`. For PolyData, pass one homogeneous
cell vector per category, as in the constructor. For temporal files this
must be called inside [`write_timestep`](@ref) (all partitions of a step
must be supplied when its geometry changes); a file created without
geometry can also add its initial partitions before the first step.

`pointdata`/`celldata` are iterables of `name => data` pairs whose lengths are
validated against this partition. Alternatively, data for all partitions of a
step can be written afterwards with `vtk[name] = data` using the concatenated
arrays.
"""
function add_partition(
        vtk::VTKHDFFile{UnstructuredState}, points, cells::UnstructuredCells;
        pointdata = (), celldata = ()
    )
    check_partition_allowed(vtk)
    kind = vtk.kind
    npoints, pts = prepare_points(points)
    topo = PartitionTopology(cells, npoints)
    # preflight everything before the first HDF5 mutation
    preflight_partition_data(vtk, pointdata, npoints, VTKPointData(), "point")
    preflight_partition_data(vtk, celldata, n_cells(topo), VTKCellData(), "cell")
    root = vtk.root
    check_points_eltype(root, eltype(pts))
    mutating(vtk) do
        append_rows(appendable(vtk, root, "Points", eltype(pts), (3,)), pts)
        append_rows(appendable(vtk, root, "Connectivity", Int64, ()), topo.connectivity)
        append_rows(appendable(vtk, root, "Offsets", Int64, ()), topo.offsets)
        append_rows(appendable(vtk, root, "Types", UInt8, ()), topo.types)
        append_rows(appendable(vtk, root, "NumberOfPoints", Int64, ()), Int64(npoints))
        append_rows(appendable(vtk, root, "NumberOfCells", Int64, ()), Int64(n_cells(topo)))
        append_rows(appendable(vtk, root, "NumberOfConnectivityIds", Int64, ()), Int64(length(topo.connectivity)))
        if topo.n_polyhedra > 0 && !kind.has_polyhedra
            enable_polyhedra!(vtk)
        end
        if kind.has_polyhedra
            append_polyhedron_partition!(vtk, topo)
        end
        kind.total_points += npoints
        kind.total_cells += n_cells(topo)
        kind.total_conn += length(topo.connectivity)
        kind.total_parts += 1
        push!(kind.part_cells, n_cells(topo))
        if vtk.in_step
            kind.step_parts += 1
            kind.step_points += npoints
            kind.step_cells += n_cells(topo)
        else
            # geometry added at construction time; usable (and reused) by steps
            kind.geom = geometry_ref(kind, (0, 0, 0, 0, 0, 0, 0))
        end
        write_partition_data(vtk, pointdata, VTKPointData(), npoints)
        write_partition_data(vtk, celldata, VTKCellData(), n_cells(topo))
    end
    return vtk
end

function preflight_partition_data(vtk, pairs, expected::Int, loc, what::String)
    for (name, data) in pairs
        check_name(String(name))
        n = tuple_count(data)
        n == expected || error("partition $what data $name has $n tuples, expected $expected")
        check_tuple_data(vtk, location_group(loc), String(name), data)
    end
    return nothing
end

function check_points_eltype(root, ::Type{T}) where {T}
    if haskey(root, "Points")
        E = eltype(get_dataset(root, "Points"))
        E == T || error("points change element type ($E -> $T); convert the coordinates explicitly")
    end
    return nothing
end

function check_partition_allowed(vtk::VTKHDFFile)
    vtk.isopen || error("file is closed")
    vtk.failed && error("a previous write to this file failed; the file is incomplete")
    if vtk.temporal && !vtk.in_step && vtk.nsteps > 0
        error("temporal file: partitions must be added inside write_timestep")
    end
    return nothing
end

function write_partition_data(vtk, pairs, loc, expected::Int)
    for (name, data) in pairs
        set_data!(vtk, data, String(name), loc; attribute = nothing)
    end
    return nothing
end

# An unstructured file closed without any partition still needs the full
# (empty) layout: write one zero-sized partition.
function finalize_kind!(vtk, kind::UnstructuredState)
    if kind.total_parts == 0
        add_partition(vtk, zeros(Float64, 3, 0), Union{MeshCell{VTKCellType}, VTKPolyhedron}[])
    end
    return nothing
end

# Geometry offsets describing the partitions appended after `snap` (a totals
# snapshot taken before they were written).
function geometry_ref(kind::UnstructuredState, snap::NTuple{7, Int})
    return (
        part_offset = snap[1], nparts = kind.total_parts - snap[1],
        point_offset = snap[2], npoints = kind.total_points - snap[2],
        cell_offset = snap[3], ncells = kind.total_cells - snap[3],
        conn_offset = snap[4],
        face_offset = snap[5], face_conn_offset = snap[6], p2f_offset = snap[7],
    )
end

function enable_polyhedra!(vtk::VTKHDFFile{UnstructuredState})
    kind = vtk.kind
    kind.has_polyhedra = true
    bump_version!(vtk, (2, 5))
    root = vtk.root
    fc = create_appendable(vtk, root, "FaceConnectivity", Int64, ())
    fo = create_appendable(vtk, root, "FaceOffsets", Int64, ())
    p2f = create_appendable(vtk, root, "PolyhedronToFaces", Int64, ())
    po = create_appendable(vtk, root, "PolyhedronOffsets", Int64, ())
    nf = create_appendable(vtk, root, "NumberOfFaces", Int64, ())
    np2f = create_appendable(vtk, root, "NumberOfPolyhedronToFaceIds", Int64, ())
    nfc = create_appendable(vtk, root, "NumberOfFaceConnectivityIds", Int64, ())
    # backfill empty polyhedron information for the partitions written so far
    for ncells in kind.part_cells
        append_rows(fo, Int64[0])
        append_rows(po, zeros(Int64, ncells + 1))
        append_rows(nf, Int64(0))
        append_rows(np2f, Int64(0))
        append_rows(nfc, Int64(0))
    end
    # and zero step offsets for already-completed time steps
    if vtk.temporal && vtk.nsteps > 0
        sg = steps_group(vtk)
        zs = zeros(Int64, vtk.nsteps)
        append_rows(appendable(vtk, sg, "FaceConnectivityOffsets", Int64, ()), zs)
        append_rows(appendable(vtk, sg, "FaceOffsetsOffsets", Int64, ()), zs)
        append_rows(appendable(vtk, sg, "PolyhedronToFaceIdOffsets", Int64, ()), zs)
    end
    return nothing
end

function append_polyhedron_partition!(vtk::VTKHDFFile{UnstructuredState}, topo::PartitionTopology)
    kind = vtk.kind
    root = vtk.root
    append_rows(get_dataset(root, "FaceConnectivity"), topo.face_connectivity)
    append_rows(get_dataset(root, "FaceOffsets"), topo.face_offsets)
    append_rows(get_dataset(root, "PolyhedronToFaces"), topo.polyhedron_to_faces)
    append_rows(get_dataset(root, "PolyhedronOffsets"), topo.polyhedron_offsets)
    append_rows(get_dataset(root, "NumberOfFaces"), Int64(n_faces(topo)))
    append_rows(get_dataset(root, "NumberOfPolyhedronToFaceIds"), Int64(length(topo.polyhedron_to_faces)))
    append_rows(get_dataset(root, "NumberOfFaceConnectivityIds"), Int64(length(topo.face_connectivity)))
    kind.total_faces += n_faces(topo)
    kind.total_face_conn += length(topo.face_connectivity)
    kind.total_p2f += length(topo.polyhedron_to_faces)
    return nothing
end

# ---- data ----

step_points(kind::UnstructuredState) =
    kind.step_parts > 0 ? kind.step_points :
    kind.geom === nothing ? 0 : kind.geom.npoints
step_cells(kind::UnstructuredState) =
    kind.step_parts > 0 ? kind.step_cells :
    kind.geom === nothing ? 0 : kind.geom.ncells

function resolve_location(vtk, kind::UnstructuredState, data)
    n = tuple_count(data)
    p = vtk.in_step ? step_points(kind) : kind.total_points
    c = vtk.in_step ? step_cells(kind) : kind.total_cells
    p == 0 && c == 0 && error("no geometry has been written yet; add points/cells first or pass an explicit data location")
    if n == p && n == c
        error("data length $n matches both the number of points and of cells; pass VTKPointData() or VTKCellData() explicitly")
    end
    n == p && return VTKPointData()
    n == c && return VTKCellData()
    error("data length $n matches neither the number of points ($p) nor of cells ($c); for field data pass VTKFieldData()")
end

write_array!(vtk::VTKHDFFile, kind::UnstructuredState, loc::AbstractFieldData, name::AbstractString, data) =
    append_tuple_data!(vtk, location_group(loc), name, data)

expected_totals(vtk, kind::UnstructuredState) =
    Dict("PointData" => kind.total_points, "CellData" => kind.total_cells)

# ---- temporal hooks ----

uses_data_offsets(::UnstructuredState) = true

function kind_begin_step!(vtk, kind::UnstructuredState)
    kind.snap = totals_snapshot(kind)
    kind.step_parts = 0
    kind.step_points = 0
    kind.step_cells = 0
    return nothing
end

function supply_step_geometry!(vtk, kind::UnstructuredState; points = nothing, cells = nothing)
    (points === nothing || cells === nothing) &&
        throw(ArgumentError("both points and cells must be given to update unstructured geometry"))
    add_partition(vtk, points, cells)
    return nothing
end

function finish_step_geometry!(vtk, kind::UnstructuredState)
    if kind.step_parts > 0
        kind.geom = geometry_ref(kind, kind.snap)
    end
    geom = kind.geom
    geom === nothing &&
        error("no geometry for this time step: pass points/cells to write_timestep or to the file constructor")
    kind.cum_points += geom.npoints
    kind.cum_cells += geom.ncells
    return geom
end

step_expected_totals(vtk, kind::UnstructuredState, geom) =
    Dict("PointData" => kind.cum_points, "CellData" => kind.cum_cells)

function append_step_offsets!(vtk, kind::UnstructuredState, sg, geom)
    append_rows(appendable(vtk, sg, "PartOffsets", Int64, ()), Int64(geom.part_offset))
    append_rows(appendable(vtk, sg, "NumberOfParts", Int64, ()), Int64(geom.nparts))
    append_rows(appendable(vtk, sg, "PointOffsets", Int64, ()), Int64(geom.point_offset))
    append_rows(appendable(vtk, sg, "CellOffsets", Int64, (1,)), reshape(Int64[geom.cell_offset], 1, 1))
    append_rows(appendable(vtk, sg, "ConnectivityIdOffsets", Int64, (1,)), reshape(Int64[geom.conn_offset], 1, 1))
    if kind.has_polyhedra
        append_rows(appendable(vtk, sg, "FaceConnectivityOffsets", Int64, ()), Int64(geom.face_conn_offset))
        append_rows(appendable(vtk, sg, "FaceOffsetsOffsets", Int64, ()), Int64(geom.face_offset))
        append_rows(appendable(vtk, sg, "PolyhedronToFaceIdOffsets", Int64, ()), Int64(geom.p2f_offset))
    end
    return nothing
end
