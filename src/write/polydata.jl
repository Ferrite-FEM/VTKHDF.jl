# PolyData writer. Cells are grouped in four fixed categories written in the
# spec-mandated order Vertices, Lines, Polygons, Strips; cell data must be
# supplied in that concatenated order.

const POLY_CATEGORIES = ("Vertices", "Lines", "Polygons", "Strips")

category_index(::Type{PolyData.Verts}) = 1
category_index(::Type{PolyData.Lines}) = 2
category_index(::Type{PolyData.Polys}) = 3
category_index(::Type{PolyData.Strips}) = 4

const PolyDataCell = MeshCell{<:PolyData.CellType}

mutable struct PolyDataState <: DatasetKind
    total_points::Int
    total_cells::NTuple{4, Int}
    total_conn::NTuple{4, Int}
    total_parts::Int
    step_parts::Int
    step_points::Int
    step_cells::Int
    snap::NTuple{10, Int}  # (parts, points, cells[4], conn[4]) at step begin
    geom::Union{Nothing, NamedTuple}
    cum_points::Int
    cum_cells::Int
end

new_polydata_state() = PolyDataState(
    0, (0, 0, 0, 0), (0, 0, 0, 0), 0,
    0, 0, 0, ntuple(_ -> 0, 10), nothing, 0, 0
)

function init_polydata(
        dest, points, cellvecs; temporal = false,
        compress = false, chunk_size = 0
    )
    file, root = open_dest(dest)
    write_ascii_attribute(root, "Type", "PolyData")
    write_version_attribute(root, (2, 0))
    vtk = make_vtkfile(
        file, root, new_polydata_state();
        temporal, compress, chunk_size, version = (2, 0)
    )
    if points !== nothing
        add_partition(vtk, points, cellvecs...)
    end
    return vtk
end

totals_snapshot(k::PolyDataState) =
    (k.total_parts, k.total_points, k.total_cells..., k.total_conn...)

# Split user-provided cell vectors into the four categories; each vector must
# be homogeneous in category.
function poly_categorize(cellvecs)
    cats = ntuple(_ -> Vector{PolyDataCell}(), 4)
    for v in cellvecs
        isempty(v) && continue
        C = typeof(first(v).ctype)
        i = category_index(C)
        for c in v
            c.ctype isa C || throw(
                ArgumentError(
                    "mixed PolyData cell categories in one vector; pass separate vectors per category " *
                        "(the on-disk cell order is Vertices, Lines, Polygons, Strips)"
                )
            )
            push!(cats[i], c)
        end
    end
    return cats
end

function add_partition(
        vtk::VTKHDFFile{PolyDataState}, points,
        cellvecs::Vararg{AbstractVector{<:PolyDataCell}};
        pointdata = (), celldata = ()
    )
    check_partition_allowed(vtk)
    kind = vtk.kind
    npoints, pts = prepare_points(points)
    cats = poly_categorize(cellvecs)
    topos = map(cat -> PartitionTopology(cat, npoints), cats)
    ncells_part = sum(n_cells, topos)
    preflight_partition_data(vtk, pointdata, npoints, VTKPointData(), "point")
    preflight_partition_data(vtk, celldata, ncells_part, VTKCellData(), "cell")
    root = vtk.root
    check_points_eltype(root, eltype(pts))
    mutating(vtk) do
        append_rows(appendable(vtk, root, "Points", eltype(pts), (3,)), pts)
        append_rows(appendable(vtk, root, "NumberOfPoints", Int64, ()), Int64(npoints))
        newcells = Int[]
        newconn = Int[]
        for (i, catname) in enumerate(POLY_CATEGORIES)
            grp = get_or_create_group(root, catname)
            topo = topos[i]
            append_rows(appendable(vtk, grp, "Connectivity", Int64, ()), topo.connectivity)
            append_rows(appendable(vtk, grp, "Offsets", Int64, ()), topo.offsets)
            append_rows(appendable(vtk, grp, "NumberOfCells", Int64, ()), Int64(n_cells(topo)))
            append_rows(appendable(vtk, grp, "NumberOfConnectivityIds", Int64, ()), Int64(length(topo.connectivity)))
            push!(newcells, n_cells(topo))
            push!(newconn, length(topo.connectivity))
        end
        kind.total_points += npoints
        kind.total_cells = ntuple(i -> kind.total_cells[i] + newcells[i], 4)
        kind.total_conn = ntuple(i -> kind.total_conn[i] + newconn[i], 4)
        kind.total_parts += 1
        if vtk.in_step
            kind.step_parts += 1
            kind.step_points += npoints
            kind.step_cells += ncells_part
        else
            kind.geom = geometry_ref(kind, ntuple(_ -> 0, 10))
        end
        write_partition_data(vtk, pointdata, VTKPointData(), npoints)
        write_partition_data(vtk, celldata, VTKCellData(), ncells_part)
    end
    return vtk
end

function geometry_ref(kind::PolyDataState, snap::NTuple{10, Int})
    cells_off = ntuple(i -> snap[2 + i], 4)
    conn_off = ntuple(i -> snap[6 + i], 4)
    return (
        part_offset = snap[1], nparts = kind.total_parts - snap[1],
        point_offset = snap[2], npoints = kind.total_points - snap[2],
        cell_offsets = cells_off,
        ncells = sum(kind.total_cells) - sum(cells_off),
        conn_offsets = conn_off,
    )
end

# A polydata file closed without any partition still needs the full (empty)
# layout: write one zero-sized partition.
function finalize_kind!(vtk, kind::PolyDataState)
    if kind.total_parts == 0
        add_partition(vtk, zeros(Float64, 3, 0))
    end
    return nothing
end

# ---- data ----

step_points(kind::PolyDataState) =
    kind.step_parts > 0 ? kind.step_points :
    kind.geom === nothing ? 0 : kind.geom.npoints
step_cells(kind::PolyDataState) =
    kind.step_parts > 0 ? kind.step_cells :
    kind.geom === nothing ? 0 : kind.geom.ncells

function resolve_location(vtk, kind::PolyDataState, data)
    n = tuple_count(data)
    p = vtk.in_step ? step_points(kind) : kind.total_points
    c = vtk.in_step ? step_cells(kind) : sum(kind.total_cells)
    p == 0 && c == 0 && error("no geometry has been written yet; add points/cells first or pass an explicit data location")
    if n == p && n == c
        error("data length $n matches both the number of points and of cells; pass VTKPointData() or VTKCellData() explicitly")
    end
    n == p && return VTKPointData()
    n == c && return VTKCellData()
    error("data length $n matches neither the number of points ($p) nor of cells ($c); for field data pass VTKFieldData()")
end

write_array!(vtk::VTKHDFFile, kind::PolyDataState, loc::AbstractFieldData, name::AbstractString, data) =
    append_tuple_data!(vtk, location_group(loc), name, data)

expected_totals(vtk, kind::PolyDataState) =
    Dict("PointData" => kind.total_points, "CellData" => sum(kind.total_cells))

# ---- temporal hooks ----

uses_data_offsets(::PolyDataState) = true

function kind_begin_step!(vtk, kind::PolyDataState)
    kind.snap = totals_snapshot(kind)
    kind.step_parts = 0
    kind.step_points = 0
    kind.step_cells = 0
    return nothing
end

function supply_step_geometry!(vtk, kind::PolyDataState; points = nothing, cells = nothing)
    (points === nothing || cells === nothing) &&
        throw(ArgumentError("both points and cells must be given to update polydata geometry"))
    cellvecs = cells isa AbstractVector{<:PolyDataCell} ? (cells,) : Tuple(cells)
    add_partition(vtk, points, cellvecs...)
    return nothing
end

function finish_step_geometry!(vtk, kind::PolyDataState)
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

step_expected_totals(vtk, kind::PolyDataState, geom) =
    Dict("PointData" => kind.cum_points, "CellData" => kind.cum_cells)

function append_step_offsets!(vtk, kind::PolyDataState, sg, geom)
    append_rows(appendable(vtk, sg, "PartOffsets", Int64, ()), Int64(geom.part_offset))
    append_rows(appendable(vtk, sg, "NumberOfParts", Int64, ()), Int64(geom.nparts))
    append_rows(appendable(vtk, sg, "PointOffsets", Int64, ()), Int64(geom.point_offset))
    append_rows(
        appendable(vtk, sg, "CellOffsets", Int64, (4,)),
        reshape(Int64[geom.cell_offsets...], 4, 1)
    )
    append_rows(
        appendable(vtk, sg, "ConnectivityIdOffsets", Int64, (4,)),
        reshape(Int64[geom.conn_offsets...], 4, 1)
    )
    return nothing
end
