# PolyData reader.
#
# Cells live in the four category groups (Vertices/Lines/Polygons/Strips),
# each a concatenation of per-partition blocks like an unstructured grid.
# Cell data, in contrast, is partition-major: each partition contributes all
# of its cells in category order. partition_ranges exposes that mapping.

const POLY_CATEGORY_FIELDS = (:vertices, :lines, :polygons, :strips)
const POLY_CATEGORY_TYPES = (PolyData.Verts(), PolyData.Lines(), PolyData.Polys(), PolyData.Strips())

struct ReadPolyData <: ReaderKind
    points_prefix::Vector{Int}
    cells_prefix::NTuple{4, Vector{Int}}
    conn_prefix::NTuple{4, Vector{Int}}
end

reader_type_string(::ReadPolyData) = "PolyData"

function read_polydata_kind(root::HDF5.Group)
    pp = count_prefix(root, "NumberOfPoints")
    nparts = npartitions_prefix(pp)
    check_rows(root, "Points", pp[end])
    zeros_prefix = zeros(Int, nparts + 1)
    cps = Vector{Int}[]
    nps = Vector{Int}[]
    for cat in POLY_CATEGORIES
        if haskey(root, cat)
            grp = root[cat]::HDF5.Group
            cp = count_prefix(grp, "NumberOfCells")
            np = count_prefix(grp, "NumberOfConnectivityIds")
            (npartitions_prefix(cp) == nparts && npartitions_prefix(np) == nparts) ||
                error("category $cat has inconsistent partition counts")
            check_rows(grp, "Connectivity", np[end])
            check_rows(grp, "Offsets", cp[end] + nparts)
            push!(cps, cp)
            push!(nps, np)
        else
            # a category with no cells may be omitted entirely
            push!(cps, zeros_prefix)
            push!(nps, zeros_prefix)
        end
    end
    return ReadPolyData(pp, Tuple(cps), Tuple(nps))
end

poly_total_cells(k::ReadPolyData, p::Int) = sum(cp -> cp[p], k.cells_prefix)

function step_slice(r::VTKHDFReader{ReadPolyData}, i::Int)
    check_open(r)
    k = r.kind
    si = steps_info(r)
    parts = step_parts(r, i)
    p1 = first(parts)
    npts = k.points_prefix[last(parts) + 1] - k.points_prefix[p1]
    ncls = poly_total_cells(k, last(parts) + 1) - poly_total_cells(k, p1)
    point_offset = check_step_offset(
        "PointOffsets", i,
        haskey(si.datasets, "PointOffsets") ? steps_entry(si, "PointOffsets", i, 0) : nothing,
        k.points_prefix[p1]
    )
    # explicit per-category topology offsets must agree with the partition
    # layout (see check_step_offset); the scalar cell-data offset is their sum
    for (key, prefixes) in (("CellOffsets", k.cells_prefix), ("ConnectivityIdOffsets", k.conn_prefix))
        offsets = steps_column(si, key, i, 4)
        offsets === nothing && continue
        for ci in 1:4
            check_step_offset("$key[$ci]", i, offsets[ci], prefixes[ci][p1])
        end
    end
    cell_offset = poly_total_cells(k, p1)
    return PartitionSlice(parts, point_offset, npts, cell_offset, ncls)
end

static_slice(k::ReadPolyData) = PartitionSlice(
    1:npartitions_prefix(k.points_prefix), 0, k.points_prefix[end], 0,
    poly_total_cells(k, length(k.points_prefix))
)

# ---- geometry accessors ----

function read_points(r::VTKHDFReader{ReadPolyData})
    check_static_geometry(r, "read_points")
    return read_points_slice(r, static_slice(r.kind))
end
function read_points(s::VTKHDFTimeStep{<:VTKHDFReader{ReadPolyData}})
    return read_points_slice(s.reader, step_slice(s.reader, s.index))
end

function read_cells(r::VTKHDFReader{ReadPolyData})
    check_static_geometry(r, "read_cells")
    return build_poly_cells(r, static_slice(r.kind).parts)
end
function read_cells(s::VTKHDFTimeStep{<:VTKHDFReader{ReadPolyData}})
    return build_poly_cells(s.reader, step_slice(s.reader, s.index).parts)
end

function build_poly_cells(r::VTKHDFReader{ReadPolyData}, parts::UnitRange{Int})
    k = r.kind
    point_base = k.points_prefix[first(parts)]
    cats = map(1:4) do ci
        ctype = POLY_CATEGORY_TYPES[ci]
        cells = MeshCell{typeof(ctype), Vector{Int64}}[]
        cp = k.cells_prefix[ci]
        np = k.conn_prefix[ci]
        cp[last(parts) + 1] - cp[first(parts)] == 0 && return cells
        grp = r.root[POLY_CATEGORIES[ci]]::HDF5.Group
        conn_ds = get_dataset(grp, "Connectivity")
        off_ds = get_dataset(grp, "Offsets")
        for p in parts
            c0, c1 = cp[p], cp[p + 1]
            c1 - c0 == 0 && continue
            offs = off_ds[(c0 + p):(c1 + p)]                   # ncells_p + 1 entries
            conn = conn_ds[(np[p] + 1):np[p + 1]]
            rebase = k.points_prefix[p] - point_base + 1
            for i in 1:(c1 - c0)
                ids = Int64[id + rebase for id in view(conn, (offs[i] + 1):offs[i + 1])]
                push!(cells, MeshCell(ctype, ids))
            end
        end
        return cells
    end
    return NamedTuple{POLY_CATEGORY_FIELDS}(Tuple(cats))
end

# ---- counts / partition structure ----

function npoints(r::VTKHDFReader{ReadPolyData})
    check_static_geometry(r, "npoints")
    return r.kind.points_prefix[end]
end
function ncells(r::VTKHDFReader{ReadPolyData})
    check_static_geometry(r, "ncells")
    return poly_total_cells(r.kind, length(r.kind.points_prefix))
end
function npartitions(r::VTKHDFReader{ReadPolyData})
    check_static_geometry(r, "npartitions")
    return npartitions_prefix(r.kind.points_prefix)
end
npoints(s::VTKHDFTimeStep{<:VTKHDFReader{ReadPolyData}}) =
    step_slice(s.reader, s.index).npoints
ncells(s::VTKHDFTimeStep{<:VTKHDFReader{ReadPolyData}}) =
    step_slice(s.reader, s.index).ncells
npartitions(s::VTKHDFTimeStep{<:VTKHDFReader{ReadPolyData}}) =
    length(step_slice(s.reader, s.index).parts)

function partition_ranges(r::VTKHDFReader{ReadPolyData})
    check_static_geometry(r, "partition_ranges")
    return partition_ranges_slice(r.kind, static_slice(r.kind).parts)
end
partition_ranges(s::VTKHDFTimeStep{<:VTKHDFReader{ReadPolyData}}) =
    partition_ranges_slice(s.reader.kind, step_slice(s.reader, s.index).parts)

function partition_ranges_slice(k::ReadPolyData, parts::UnitRange{Int})
    pbase = k.points_prefix[first(parts)]
    points = [(k.points_prefix[p] - pbase + 1):(k.points_prefix[p + 1] - pbase) for p in parts]
    cells = UnitRange{Int}[]
    by_category = NamedTuple{POLY_CATEGORY_FIELDS, NTuple{4, UnitRange{Int}}}[]
    pos = 0  # partition-major position into the cell-data arrays
    for p in parts
        start = pos
        catranges = UnitRange{Int}[]
        for ci in 1:4
            n = k.cells_prefix[ci][p + 1] - k.cells_prefix[ci][p]
            push!(catranges, (pos + 1):(pos + n))
            pos += n
        end
        push!(cells, (start + 1):pos)
        push!(by_category, NamedTuple{POLY_CATEGORY_FIELDS}(Tuple(catranges)))
    end
    return (points = points, cells = cells, cells_by_category = by_category)
end

# ---- per-step data arrays ----

function read_step_array(
        r::VTKHDFReader{ReadPolyData}, kind::ReadPolyData,
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
