# OverlappingAMR reader (static; temporal AMR is rejected at open).

struct ReadAMR <: ReaderKind
    nlevels::Int
    origin::NTuple{3, Float64}
    grid_description::String
end

reader_type_string(::ReadAMR) = "OverlappingAMR"
supports_temporal(::ReadAMR) = false

function read_amr_kind(root::HDF5.Group)
    origin = read_attr_tuple(root, "Origin", Float64, Val(3))
    origin === nothing && error("$(HDF5.name(root)) has no Origin attribute")
    desc = haskey(HDF5.attrs(root), "GridDescription") ?
        rstrip(String(HDF5.read_attribute(root, "GridDescription")), '\0') : "XYZ"
    nlevels = 0
    while haskey(root, "Level$nlevels")
        nlevels += 1
    end
    return ReadAMR(nlevels, origin, desc)
end

grid_info(r::VTKHDFReader{ReadAMR}) =
    (origin = r.kind.origin, grid_description = r.kind.grid_description)

"""
    nlevels(r) -> Int

Number of refinement levels of an OverlappingAMR file.
"""
nlevels(r::VTKHDFReader{ReadAMR}) = (check_open(r); r.kind.nlevels)

Base.getindex(r::VTKHDFReader{ReadAMR}, name::AbstractString) =
    error("OverlappingAMR data is stored per level; use amr_level(r, l)[$(repr(String(name)))]")
Base.getindex(r::VTKHDFReader{ReadAMR}, name::AbstractString, loc::AbstractFieldData) =
    error("OverlappingAMR data is stored per level; use amr_level(r, l)[$(repr(String(name))), loc]")

"""
    amr_level(r, l) -> level

The refinement level `l` of an OverlappingAMR file, with `l` the on-disk
level number in `0:nlevels(r)-1` (`Level0` is the coarsest). The level
supports [`level_info`](@ref), `keys(level, loc)`, data reading by indexing
(`level["ρ"]`, `level["ρ", VTKCellData()]` — arrays are concatenated over
the level's boxes) and [`partition_ranges`](@ref) with one point/cell range
per box.
"""
function amr_level(r::VTKHDFReader{ReadAMR}, l::Integer)
    check_open(r)
    0 <= l < r.kind.nlevels ||
        error("level $l does not exist; the file has levels 0:$(r.kind.nlevels - 1)")
    group = r.root["Level$l"]::HDF5.Group
    spacing = read_attr_tuple(group, "Spacing", Float64, Val(3))
    spacing === nothing && error("Level$l has no Spacing attribute")
    boxes = NTuple{6, Int}[]
    if haskey(group, "AMRBox")
        B = read(get_dataset(group, "AMRBox"))
        size(B, 1) == 6 || error("Level$l/AMRBox has $(size(B, 1)) columns, expected 6")
        boxes = [NTuple{6, Int}(Int.(B[:, j])) for j in axes(B, 2)]
    end
    return VTKHDFAMRLevel(r, Int(l), group, spacing, boxes)
end

struct VTKHDFAMRLevel
    reader::VTKHDFReader{ReadAMR}
    level::Int
    group::HDF5.Group
    spacing::NTuple{3, Float64}
    boxes::Vector{NTuple{6, Int}}
end

function Base.show(io::IO, lvl::VTKHDFAMRLevel)
    return print(io, "VTKHDFAMRLevel ", lvl.level, " (", length(lvl.boxes), " boxes)")
end

"""
    level_info(level) -> (; spacing, boxes)

Spacing and AMR boxes (inclusive cell extents
`(imin, imax, jmin, jmax, kmin, kmax)`) of an [`amr_level`](@ref).
"""
level_info(lvl::VTKHDFAMRLevel) = (spacing = lvl.spacing, boxes = copy(lvl.boxes))

box_ncells(b::NTuple{6, Int}) = prod(b[2i] - b[2i - 1] + 1 for i in 1:3)
box_npoints(b::NTuple{6, Int}) = prod(b[2i] - b[2i - 1] + 2 for i in 1:3)

npoints(lvl::VTKHDFAMRLevel) = sum(box_npoints, lvl.boxes; init = 0)
ncells(lvl::VTKHDFAMRLevel) = sum(box_ncells, lvl.boxes; init = 0)

function partition_ranges(lvl::VTKHDFAMRLevel)
    points = UnitRange{Int}[]
    cells = UnitRange{Int}[]
    ppos = cpos = 0
    for b in lvl.boxes
        push!(points, (ppos + 1):(ppos + box_npoints(b)))
        push!(cells, (cpos + 1):(cpos + box_ncells(b)))
        ppos = last(points[end])
        cpos = last(cells[end])
    end
    return (points = points, cells = cells)
end

function Base.keys(lvl::VTKHDFAMRLevel, loc::AbstractFieldData)
    check_open(lvl.reader)
    grp = read_group_name(loc)
    haskey(lvl.group, grp) || return String[]
    return sort!(collect(String, keys(lvl.group[grp]::HDF5.Group)))
end

Base.haskey(lvl::VTKHDFAMRLevel, name::AbstractString, loc::AbstractFieldData) =
    (check_open(lvl.reader); haskey(lvl.group, read_group_name(loc) * "/" * (String(name)::String)))

function Base.getindex(lvl::VTKHDFAMRLevel, name::AbstractString)
    check_open(lvl.reader)
    hits = AbstractFieldData[]
    for loc in (VTKPointData(), VTKCellData(), VTKFieldData())
        haskey(lvl, name, loc) && push!(hits, loc)
    end
    isempty(hits) &&
        error("no data array named $(repr(String(name))) in Level$(lvl.level)")
    length(hits) > 1 && error(
        "data array name $(repr(String(name))) is ambiguous in Level$(lvl.level); " *
            "pass the location explicitly"
    )
    return lvl[name, hits[1]]
end

function Base.getindex(lvl::VTKHDFAMRLevel, name::AbstractString, loc::AbstractFieldData)
    check_open(lvl.reader)
    haskey(lvl, name, loc) ||
        error("no data array named $(repr(String(name))) in Level$(lvl.level)/$(read_group_name(loc))")
    ds = get_dataset(lvl.group[read_group_name(loc)]::HDF5.Group, String(name))
    return read_static_array(lvl.reader.kind, loc, ds)
end
