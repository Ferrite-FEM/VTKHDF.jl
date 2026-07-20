# OverlappingAMR writer (static).

mutable struct AMRLevel
    vtk::Any            # ::VTKHDFFile{AMRState} (untyped to avoid circular definition)
    group::HDF5.Group
    total_points::Int
    total_cells::Int
    data_rows::Dict{String, Int}
end

mutable struct AMRState <: DatasetKind
    levels::Vector{AMRLevel}
end

"""
    vtkhdf_amr(filename; origin, grid_description = :XYZ, compress = false)

Create a VTKHDF `OverlappingAMR` file. Levels are added with
[`add_level`](@ref) and boxes with [`add_box`](@ref):

```julia
vtkhdf_amr("amr"; origin = (0, 0, 0)) do amr
    lvl = add_level(amr; spacing = (1.0, 1.0, 1.0))
    add_box(lvl, (0, 4, 0, 4, 0, 4); celldata = ("ρ" => ρ,))
end
```

Temporal AMR is not supported.
"""
function vtkhdf_amr(filename::AbstractString; kwargs...)
    return init_amr(filename; kwargs...)
end

function init_amr(
        dest; origin, grid_description::Symbol = :XYZ,
        compress = false, chunk_size = 0, temporal = false
    )
    temporal && throw(ArgumentError("temporal OverlappingAMR writing is not supported"))
    grid_description === :XYZ ||
        throw(ArgumentError("only grid_description = :XYZ is supported"))
    length(origin) == 3 || throw(ArgumentError("origin must have 3 entries"))
    file, root = open_dest(dest)
    write_ascii_attribute(root, "Type", "OverlappingAMR")
    write_version_attribute(root, (2, 3))
    HDF5.attrs(root)["Origin"] = Float64[origin...]
    write_ascii_attribute(root, "GridDescription", String(grid_description))
    return make_vtkfile(
        file, root, AMRState(AMRLevel[]);
        temporal = false, compress, chunk_size, version = (2, 3)
    )
end

"""
    add_level(amr; spacing) -> level

Append a refinement level (`Level0`, `Level1`, ...) with the given cell
`spacing` to an OverlappingAMR file.
"""
function add_level(vtk::VTKHDFFile{AMRState}; spacing)
    vtk.isopen || error("file is closed")
    length(spacing) == 3 || throw(ArgumentError("spacing must have 3 entries"))
    name = "Level$(length(vtk.kind.levels))"
    grp = HDF5.create_group(vtk.root, name)
    HDF5.attrs(grp)["Spacing"] = Float64[spacing...]
    lvl = AMRLevel(vtk, grp, 0, 0, Dict{String, Int}())
    push!(vtk.kind.levels, lvl)
    return lvl
end

"""
    add_box(level, (imin, imax, jmin, jmax, kmin, kmax); pointdata = (), celldata = (), fielddata = ())

Append one AMR box (inclusive cell extents) to a level. Data arrays are
serialized per box, in the order boxes are added; sizes are validated against
the box extents (`imax - imin + 1` cells and `imax - imin + 2` points per
direction). Data can alternatively be written for a whole level at once with
`level[name, VTKCellData()] = data` after adding its boxes.
"""
function add_box(
        lvl::AMRLevel, extents::NTuple{6, Integer};
        pointdata = (), celldata = (), fielddata = ()
    )
    vtk = lvl.vtk::VTKHDFFile{AMRState}
    vtk.isopen || error("file is closed")
    vtk.failed && error("a previous write to this file failed; the file is incomplete")
    (imin, imax, jmin, jmax, kmin, kmax) = Int64.(extents)
    (imax >= imin && jmax >= jmin && kmax >= kmin) ||
        throw(ArgumentError("invalid AMR box extents $extents"))
    ncells = Base.checked_mul(imax - imin + 1, jmax - jmin + 1, kmax - kmin + 1)
    npoints = Base.checked_mul(imax - imin + 2, jmax - jmin + 2, kmax - kmin + 2)
    # preflight all data before mutating the file
    for (name, data) in pointdata
        check_name(String(name))
        n = tuple_count(data)
        n == npoints || error("box point data $name has $n tuples, expected $npoints")
        prepare_tuples(data)
    end
    for (name, data) in celldata
        check_name(String(name))
        n = tuple_count(data)
        n == ncells || error("box cell data $name has $n tuples, expected $ncells")
        prepare_tuples(data)
    end
    for (name, data) in fielddata
        check_name(String(name))
        prepare_tuples(data)
    end
    mutating(vtk) do
        append_rows(
            appendable(vtk, lvl.group, "AMRBox", Int64, (6,)),
            reshape(Int64[imin, imax, jmin, jmax, kmin, kmax], 6, 1)
        )
        lvl.total_points += npoints
        lvl.total_cells += ncells
        for (name, data) in pointdata
            level_append_data!(lvl, "PointData", String(name), data)
        end
        for (name, data) in celldata
            level_append_data!(lvl, "CellData", String(name), data)
        end
        for (name, data) in fielddata
            level_append_data!(lvl, "FieldData", String(name), data)
        end
    end
    return lvl
end

function level_append_data!(lvl::AMRLevel, groupname::String, name::AbstractString, data)
    vtk = lvl.vtk::VTKHDFFile{AMRState}
    vtk.isopen || error("file is closed")
    vtk.failed && error("a previous write to this file failed; the file is incomplete")
    check_name(name)
    ncomp, n, arr = prepare_tuples(data)
    grp = get_or_create_group(lvl.group, groupname)
    ds = haskey(grp, name) ? get_dataset(grp, name) :
        create_appendable(vtk, grp, name, eltype(arr), ncomp == 1 ? () : (ncomp,))
    append_rows(ds, arr)
    key = groupname * "/" * name
    lvl.data_rows[key] = get(lvl.data_rows, key, 0) + n
    return nothing
end

function Base.setindex!(lvl::AMRLevel, data, name::AbstractString, loc::AbstractFieldData)
    level_append_data!(lvl, location_group(loc), name, data)
    return data
end

function Base.setindex!(lvl::AMRLevel, data, name::AbstractString)
    n = tuple_count(data)
    p, c = lvl.total_points, lvl.total_cells
    if n == p && n == c
        error("data length $n matches both points and cells of the level; pass VTKPointData() or VTKCellData()")
    end
    loc = n == p ? VTKPointData() : n == c ? VTKCellData() :
        error("data length $n matches neither level points ($p) nor cells ($c)")
    return setindex!(lvl, data, name, loc)
end

function finalize_kind!(vtk, kind::AMRState)
    for (i, lvl) in enumerate(kind.levels)
        for (key, rows) in lvl.data_rows
            groupname = first(split(key, '/'))
            want = groupname == "PointData" ? lvl.total_points :
                groupname == "CellData" ? lvl.total_cells : rows
            rows == want || error("Level$(i - 1) array $key has $rows tuples, expected $want")
        end
    end
    return nothing
end
