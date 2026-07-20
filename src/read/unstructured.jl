# UnstructuredGrid reader.
#
# On disk, partitioned datasets are concatenations of per-partition blocks
# (delimited by the NumberOf* count datasets), and connectivity uses
# partition-local 0-based point ids; polyhedron faces likewise use
# partition-local face ids. Reading rebases everything to 1-based indices
# into the concatenated points of the reader (or of the time step).

struct ReadUnstructured <: ReaderKind
    points_prefix::Vector{Int}
    cells_prefix::Vector{Int}
    conn_prefix::Vector{Int}
    has_polyhedra::Bool
    faces_prefix::Vector{Int}
    faceconn_prefix::Vector{Int}
    p2f_prefix::Vector{Int}
end

reader_type_string(::ReadUnstructured) = "UnstructuredGrid"

const POLYHEDRON_DATASETS = (
    "FaceConnectivity", "FaceOffsets", "PolyhedronToFaces", "PolyhedronOffsets",
    "NumberOfFaces", "NumberOfPolyhedronToFaceIds", "NumberOfFaceConnectivityIds",
)

function read_unstructured_kind(root::HDF5.Group)
    pp = count_prefix(root, "NumberOfPoints")
    cp = count_prefix(root, "NumberOfCells")
    np = count_prefix(root, "NumberOfConnectivityIds")
    nparts = npartitions_prefix(pp)
    (npartitions_prefix(cp) == nparts && npartitions_prefix(np) == nparts) || error(
        "inconsistent partition counts: NumberOfPoints/NumberOfCells/NumberOfConnectivityIds " *
            "have $(nparts)/$(npartitions_prefix(cp))/$(npartitions_prefix(np)) entries"
    )
    check_rows(root, "Points", pp[end])
    check_rows(root, "Types", cp[end])
    check_rows(root, "Connectivity", np[end])
    check_rows(root, "Offsets", cp[end] + nparts)
    has_polyhedra = any(name -> haskey(root, name), POLYHEDRON_DATASETS)
    if has_polyhedra
        for name in POLYHEDRON_DATASETS
            haskey(root, name) || error("incomplete polyhedron information: missing dataset $name")
        end
        fp = count_prefix(root, "NumberOfFaces")
        fcp = count_prefix(root, "NumberOfFaceConnectivityIds")
        p2fp = count_prefix(root, "NumberOfPolyhedronToFaceIds")
        (npartitions_prefix(fp) == nparts && npartitions_prefix(fcp) == nparts && npartitions_prefix(p2fp) == nparts) ||
            error("inconsistent polyhedron partition counts")
        check_rows(root, "FaceOffsets", fp[end] + nparts)
        check_rows(root, "FaceConnectivity", fcp[end])
        check_rows(root, "PolyhedronToFaces", p2fp[end])
        check_rows(root, "PolyhedronOffsets", cp[end] + nparts)
    else
        fp = fcp = p2fp = [0]
    end
    return ReadUnstructured(pp, cp, np, has_polyhedra, fp, fcp, p2fp)
end

function check_rows(root::HDF5.Group, name::String, want::Int)
    ds = require_dataset(root, name)
    n = n_rows(ds)
    n == want || error("dataset $name has $n rows, expected $want from the partition counts")
    return nothing
end

# ---- partition/step geometry bookkeeping ----

# The geometry of a read (all partitions for static files, one step's
# partitions for temporal ones): the global partition range plus the offsets
# and totals needed for slicing and data reading.
struct PartitionSlice
    parts::UnitRange{Int}
    point_offset::Int   # 0-based row offset into Points
    npoints::Int
    cell_offset::Int    # 0-based cell offset (cell-data rows, Types rows)
    ncells::Int
end

function static_slice(k::ReadUnstructured)
    nparts = npartitions_prefix(k.points_prefix)
    return PartitionSlice(1:nparts, 0, k.points_prefix[end], 0, k.cells_prefix[end])
end

function step_parts(r::VTKHDFReader, i::Int)
    k = r.kind
    si = steps_info(r)
    nparts_total = npartitions_prefix(k.points_prefix)
    nparts = steps_entry(si, "NumberOfParts", i, -1)
    if nparts < 0
        # spec: optional when the number of parts per step is constant
        (si.nsteps > 0 && nparts_total % si.nsteps == 0) || error(
            "Steps/NumberOfParts is missing and the $nparts_total partitions do not divide " *
                "evenly into $(si.nsteps) steps"
        )
        nparts = nparts_total ÷ si.nsteps
    end
    off = steps_entry(si, "PartOffsets", i, -1)
    if off < 0
        haskey(si.datasets, "NumberOfParts") && error(
            "Steps/PartOffsets is missing; it is required when Steps/NumberOfParts is present"
        )
        off = (i - 1) * nparts
    end
    parts = (off + 1):(off + nparts)
    (first(parts) >= 1 && last(parts) <= nparts_total) ||
        error("step $i references partitions $parts, but the file has $nparts_total")
    return parts
end

# The reader derives all dataset positions from the NumberOf* partition
# structure; explicit Steps offsets, when present, must agree with it (they
# are redundant in a self-consistent file). A mismatch means the offset
# tables and the partition layout tell different stories — error instead of
# silently reading the wrong slabs.
function check_step_offset(what::String, i::Int, stored, derived::Int)
    stored === nothing && return derived
    Int(stored) == derived || error(
        "Steps/$what of step $i is $stored, inconsistent with the partition layout " *
            "(expected $derived from the NumberOf* datasets)"
    )
    return derived
end

function step_slice(r::VTKHDFReader{ReadUnstructured}, i::Int)
    check_open(r)
    k = r.kind
    si = steps_info(r)
    parts = step_parts(r, i)
    p1 = first(parts)
    npts = k.points_prefix[last(parts) + 1] - k.points_prefix[p1]
    ncls = k.cells_prefix[last(parts) + 1] - k.cells_prefix[p1]
    point_offset = check_step_offset(
        "PointOffsets", i,
        haskey(si.datasets, "PointOffsets") ? steps_entry(si, "PointOffsets", i, 0) : nothing,
        k.points_prefix[p1]
    )
    cell_offsets = steps_column(si, "CellOffsets", i, 1)
    cell_offset = check_step_offset(
        "CellOffsets", i, cell_offsets === nothing ? nothing : cell_offsets[1],
        k.cells_prefix[p1]
    )
    conn_offsets = steps_column(si, "ConnectivityIdOffsets", i, 1)
    check_step_offset(
        "ConnectivityIdOffsets", i, conn_offsets === nothing ? nothing : conn_offsets[1],
        k.conn_prefix[p1]
    )
    if k.has_polyhedra
        for (key, prefix) in (
                ("FaceOffsetsOffsets", k.faces_prefix),
                ("FaceConnectivityOffsets", k.faceconn_prefix),
                ("PolyhedronToFaceIdOffsets", k.p2f_prefix),
            )
            haskey(si.datasets, key) || continue
            check_step_offset(key, i, steps_entry(si, key, i, 0), prefix[p1])
        end
    end
    return PartitionSlice(parts, point_offset, npts, cell_offset, ncls)
end

# ---- geometry accessors ----

function read_points(r::VTKHDFReader{ReadUnstructured})
    check_static_geometry(r, "read_points")
    return read_points_slice(r, static_slice(r.kind))
end
function read_points(s::VTKHDFTimeStep{<:VTKHDFReader{ReadUnstructured}})
    slice = step_slice(s.reader, s.index)
    return read_points_slice(s.reader, slice)
end

read_points_slice(r::VTKHDFReader, slice::PartitionSlice) = read_tuple_rows(
    require_dataset(r.root, "Points"),
    (slice.point_offset + 1):(slice.point_offset + slice.npoints)
)

function read_cells(r::VTKHDFReader{ReadUnstructured})
    check_static_geometry(r, "read_cells")
    return build_cells(r, static_slice(r.kind).parts)
end
function read_cells(s::VTKHDFTimeStep{<:VTKHDFReader{ReadUnstructured}})
    return build_cells(s.reader, step_slice(s.reader, s.index).parts)
end

function build_cells(r::VTKHDFReader{ReadUnstructured}, parts::UnitRange{Int})
    k = r.kind
    root = r.root
    conn_ds = require_dataset(root, "Connectivity")
    off_ds = require_dataset(root, "Offsets")
    types_ds = require_dataset(root, "Types")
    cells = Union{MeshCell{VTKCellType}, VTKPolyhedron}[]
    sizehint!(cells, k.cells_prefix[last(parts) + 1] - k.cells_prefix[first(parts)])
    point_base = k.points_prefix[first(parts)]
    for p in parts
        c0, c1 = k.cells_prefix[p], k.cells_prefix[p + 1]
        ncells_p = c1 - c0
        ncells_p == 0 && continue
        offs = off_ds[(c0 + p):(c1 + p)]                       # ncells_p + 1 entries
        types = types_ds[(c0 + 1):c1]
        conn = conn_ds[(k.conn_prefix[p] + 1):k.conn_prefix[p + 1]]
        rebase = k.points_prefix[p] - point_base + 1           # 0-based local -> 1-based slice ids
        poly = k.has_polyhedra ? partition_polyhedron_data(r, p, c0, c1) : nothing
        for i in 1:ncells_p
            ids = Int64[id + rebase for id in view(conn, (offs[i] + 1):offs[i + 1])]
            t = Int(types[i])
            if t == Int(VTKCellTypes.VTK_POLYHEDRON.vtk_id)
                poly === nothing && error(
                    "cell $(c0 + i) is a polyhedron but the file has no polyhedron datasets"
                )
                faceids = view(poly.p2f, (poly.po[i] + 1):poly.po[i + 1])  # 0-based partition-local
                faces = [
                    Int64[id + rebase for id in view(poly.fconn, (poly.fo[fid + 1] + 1):poly.fo[fid + 2])]
                        for fid in faceids
                ]
                push!(cells, VTKPolyhedron(ids, faces...))
            else
                ct = try
                    VTKCellTypes.VTKCellType(t)
                catch
                    error("unknown VTK cell type id $t (cell $(c0 + i) of the file)")
                end
                push!(cells, MeshCell(ct, ids))
            end
        end
    end
    return cells
end

# The polyhedron datasets of partition p, sliced (fo has nfaces_p + 1
# entries, po has ncells_p + 1; both hold partition-local 0-based values).
function partition_polyhedron_data(r::VTKHDFReader{ReadUnstructured}, p::Int, c0::Int, c1::Int)
    k = r.kind
    root = r.root
    f0, f1 = k.faces_prefix[p], k.faces_prefix[p + 1]
    return (
        fo = get_dataset(root, "FaceOffsets")[(f0 + p):(f1 + p)],
        po = get_dataset(root, "PolyhedronOffsets")[(c0 + p):(c1 + p)],
        p2f = get_dataset(root, "PolyhedronToFaces")[(k.p2f_prefix[p] + 1):k.p2f_prefix[p + 1]],
        fconn = get_dataset(root, "FaceConnectivity")[(k.faceconn_prefix[p] + 1):k.faceconn_prefix[p + 1]],
    )
end

# ---- counts / partition structure ----

function npoints(r::VTKHDFReader{ReadUnstructured})
    check_static_geometry(r, "npoints")
    return r.kind.points_prefix[end]
end
function ncells(r::VTKHDFReader{ReadUnstructured})
    check_static_geometry(r, "ncells")
    return r.kind.cells_prefix[end]
end
function npartitions(r::VTKHDFReader{ReadUnstructured})
    check_static_geometry(r, "npartitions")
    return npartitions_prefix(r.kind.points_prefix)
end
npoints(s::VTKHDFTimeStep{<:VTKHDFReader{ReadUnstructured}}) =
    step_slice(s.reader, s.index).npoints
ncells(s::VTKHDFTimeStep{<:VTKHDFReader{ReadUnstructured}}) =
    step_slice(s.reader, s.index).ncells
npartitions(s::VTKHDFTimeStep{<:VTKHDFReader{ReadUnstructured}}) =
    length(step_slice(s.reader, s.index).parts)

function partition_ranges(r::VTKHDFReader{ReadUnstructured})
    check_static_geometry(r, "partition_ranges")
    return partition_ranges_slice(r.kind, static_slice(r.kind).parts)
end
partition_ranges(s::VTKHDFTimeStep{<:VTKHDFReader{ReadUnstructured}}) =
    partition_ranges_slice(s.reader.kind, step_slice(s.reader, s.index).parts)

function partition_ranges_slice(k::ReadUnstructured, parts::UnitRange{Int})
    pbase = k.points_prefix[first(parts)]
    cbase = k.cells_prefix[first(parts)]
    return (
        points = [(k.points_prefix[p] - pbase + 1):(k.points_prefix[p + 1] - pbase) for p in parts],
        cells = [(k.cells_prefix[p] - cbase + 1):(k.cells_prefix[p + 1] - cbase) for p in parts],
    )
end

# ---- per-step data arrays ----

function read_step_array(
        r::VTKHDFReader{ReadUnstructured}, kind::ReadUnstructured,
        loc::Union{VTKPointData, VTKCellData}, name::String, ds::HDF5.Dataset, i::Int
    )
    slice = step_slice(r, i)
    if loc isa VTKPointData
        off = steps_entry(steps_info(r), "PointDataOffsets/" * name, i, slice.point_offset)
        n = slice.npoints
    else
        off = steps_entry(steps_info(r), "CellDataOffsets/" * name, i, slice.cell_offset)
        n = slice.ncells
    end
    return read_tuple_rows(ds, (off + 1):(off + n))
end
