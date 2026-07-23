# Core file object shared by all dataset kinds.

abstract type DatasetKind end

mutable struct VTKHDFFile{K <: DatasetKind}
    file::Union{HDF5.File, Nothing}  # nothing for composite blocks (collection owns the handle)
    root::HDF5.Group                # /VTKHDF or /VTKHDF/<BlockName>
    kind::K
    temporal::Bool
    compress::Int
    chunk_size::Int
    version::NTuple{2, Int}
    isopen::Bool
    in_step::Bool
    nsteps::Int
    data_rows::Dict{String, Int}        # "PointData/u" => total rows written
    step_start_rows::Dict{String, Int}  # snapshot of data_rows at step begin
    field_ncomp::Dict{String, Int}      # ncomp of FieldData arrays (for FieldDataSizes)
    schema::Union{Nothing, Set{String}}  # frozen array set after first step
    step_values::Vector{Float64}       # time values written so far
    failed::Bool                       # a write failed; the file is incomplete
end

function make_vtkfile(
        file, root, kind::DatasetKind; temporal::Bool,
        compress::Union{Bool, Integer}, chunk_size::Integer, version::NTuple{2, Int}
    )
    level = compress isa Bool ? (compress ? 6 : 0) : Int(compress)
    0 <= level <= 9 || throw(ArgumentError("compress must be false/true or a gzip level 0-9"))
    return VTKHDFFile(
        file, root, kind, temporal, level, Int(chunk_size), version,
        true, false, 0,
        Dict{String, Int}(), Dict{String, Int}(), Dict{String, Int}(), nothing, Float64[], false
    )
end

function bump_version!(vtk::VTKHDFFile, v::NTuple{2, Int})
    if v > vtk.version
        vtk.version = v
        write_version_attribute(vtk.root, v)
    end
    return nothing
end

function add_extension(filename::AbstractString)
    _, ext = splitext(filename)
    return isempty(ext) ? filename * ".vtkhdf" : String(filename)
end

function open_vtkhdf(filename::AbstractString; track_order::Bool = false)
    file = h5open(add_extension(filename), "w")
    root = HDF5.create_group(file, "VTKHDF"; track_order)
    return file, root
end

function Base.close(vtk::VTKHDFFile)
    vtk.isopen || return nothing
    vtk.in_step && error("close called inside write_timestep")
    try
        if vtk.failed
            @warn "closing a VTKHDF file after a failed write; the file is incomplete"
        elseif vtk.temporal
            if vtk.nsteps == 0
                @warn "closing temporal VTKHDF file without any time steps"
                # materialize an empty but complete Steps layout
                sg = steps_group(vtk)
                appendable(vtk, sg, "Values", Float64, ())
            end
            finalize_kind!(vtk, vtk.kind)
        else
            validate_static_data(vtk)
            finalize_kind!(vtk, vtk.kind)
        end
        write_version_attribute(vtk.root, vtk.version)
    finally
        # release the HDF5 handles even when validation/finalization throws
        vtk.isopen = false
        close(vtk.root)
        file = vtk.file
        file === nothing || close(file)
    end
    return nothing
end

function Base.flush(vtk::VTKHDFFile)
    vtk.isopen || return nothing
    file = vtk.file
    file === nothing || flush(file)  # nothing for composite blocks (collection owns the handle)
    return nothing
end

finalize_kind!(vtk, kind::DatasetKind) = nothing

function Base.show(io::IO, vtk::VTKHDFFile)
    status = vtk.isopen ? "open" : "closed"
    return print(
        io, "VTKHDFFile{", nameof(typeof(vtk.kind)), "} (", status,
        vtk.temporal ? ", temporal, $(vtk.nsteps) steps)" : ")"
    )
end

# ---- data writing ----

location_group(::VTKPointData) = "PointData"
location_group(::VTKCellData) = "CellData"
location_group(::VTKFieldData) = "FieldData"

"""
    vtk["name"] = data
    vtk["name", loc] = data
    vtk["name", loc, attribute = :Scalars] = data

Write a data array. Without a location the location is inferred from the
array size; `loc` (`VTKPointData()`, `VTKCellData()`, `VTKFieldData()`,
`VTKRowData()`) gives it explicitly. `attribute` (`:Scalars`, `:Vectors`,
...) marks the array as the location's active attribute. See the manual's
[Data arrays](@ref data-arrays) for the accepted array shapes.
"""
function Base.setindex!(vtk::VTKHDFFile, data, name::AbstractString; attribute::Union{Nothing, Symbol} = nothing)
    return set_data!(vtk, data, name, nothing; attribute)
end
function Base.setindex!(
        vtk::VTKHDFFile, data, name::AbstractString, loc::AbstractFieldData;
        attribute::Union{Nothing, Symbol} = nothing
    )
    return set_data!(vtk, data, name, loc; attribute)
end

function set_data!(vtk::VTKHDFFile, data, name::AbstractString, loc; attribute)
    vtk.isopen || error("file is closed")
    vtk.failed && error("a previous write to this file failed; the file is incomplete")
    check_name(name)
    if vtk.temporal && !vtk.in_step
        error("temporal file: data must be written inside write_timestep")
    end
    locr = loc === nothing ? resolve_location(vtk, vtk.kind, data) : loc
    if attribute !== nothing  # validated before anything is written
        locr isa VTKFieldData && throw(ArgumentError("attribute marking not supported for field data"))
        check_attribute_kind(attribute)
    end
    if locr isa VTKFieldData && (data isa AbstractString || data isa AbstractVector{<:AbstractString})
        write_string_field!(vtk, name, data)
    else
        write_array!(vtk, vtk.kind, locr, name, data)
    end
    attribute === nothing ||
        mark_attribute(vtk, vtk.root[location_group(locr)], name, attribute)
    return data
end

# String field data (allowed for static files only; VTK stores these as
# variable-length ASCII).
function write_string_field!(vtk::VTKHDFFile, name::AbstractString, data)
    vtk.temporal && throw(ArgumentError("string field data is not supported in temporal files"))
    strs = data isa AbstractString ? [String(data)] : String.(data)
    grp = field_group(vtk)
    haskey(grp, name) && throw(ArgumentError("field array $name already written"))
    write_string_dataset(grp, name, strs)
    return nothing
end

field_group(vtk::VTKHDFFile) = get_or_create_group(vtk.root, "FieldData")

# Generic append-style writing of tuple data (used by unstructured-like kinds
# for Point/Cell/RowData, and by all kinds for numeric FieldData).
function append_tuple_data!(vtk::VTKHDFFile, groupname::String, name::AbstractString, data)
    ncomp, n, arr = prepare_tuples(data)
    grp = get_or_create_group(vtk.root, groupname)
    key = groupname * "/" * name
    rowdims = ncomp == 1 ? () : (ncomp,)
    if haskey(grp, name)
        ds = get_dataset(grp, name)
        if vtk.schema !== nothing && !(key in vtk.schema)
            error("array $key was not part of the first time step; the array schema is fixed by the first step")
        end
        eltype(ds) == eltype(arr) ||
            error("array $key changes element type ($(eltype(ds)) -> $(eltype(arr)))")
    else
        if vtk.schema !== nothing
            error("array $key was not part of the first time step; the array schema is fixed by the first step")
        end
        ds = create_appendable(vtk, grp, name, eltype(arr), rowdims)
    end
    if groupname == "FieldData"
        known = get(vtk.field_ncomp, key, ncomp)
        known == ncomp || error("field array $key changes component count ($known -> $ncomp)")
        vtk.field_ncomp[key] = ncomp
    end
    append_rows(ds, arr)
    vtk.data_rows[key] = get(vtk.data_rows, key, 0) + n
    return nothing
end

# Validation half of append_tuple_data! — everything that can throw, without
# touching the file. Used to preflight multi-array operations so a failure
# cannot leave the file partially mutated.
function check_tuple_data(vtk::VTKHDFFile, groupname::String, name::AbstractString, data)
    ncomp, _, arr = prepare_tuples(data)
    key = groupname * "/" * name
    if vtk.schema !== nothing && !(key in vtk.schema)
        error("array $key was not part of the first time step; the array schema is fixed by the first step")
    end
    if haskey(vtk.root, groupname)
        grp = vtk.root[groupname]::HDF5.Group
        if haskey(grp, name)
            ds = get_dataset(grp, name)
            eltype(ds) == eltype(arr) ||
                error("array $key changes element type ($(eltype(ds)) -> $(eltype(arr)))")
        end
    end
    if groupname == "FieldData"
        known = get(vtk.field_ncomp, key, ncomp)
        known == ncomp || error("field array $key changes component count ($known -> $ncomp)")
    end
    return nothing
end

# Run `f`, which appends to the file; a throw midway leaves the on-disk state
# out of sync with the bookkeeping, so mark the file failed (all further
# writes are then rejected).
function mutating(f::Function, vtk::VTKHDFFile)
    try
        return f()
    catch
        vtk.failed = true
        rethrow()
    end
end

# Number of tuples in `data` without materializing the conversion.
tuple_count(data::AbstractVector) = length(data)
tuple_count(data::AbstractMatrix) = size(data, 2)
tuple_count(data) = -1

# Fallbacks; kinds override what they support.
resolve_location(vtk, kind::DatasetKind, data) =
    throw(ArgumentError("cannot determine data location automatically; pass VTKPointData()/VTKCellData()/VTKFieldData() explicitly"))

write_array!(vtk::VTKHDFFile, kind::DatasetKind, loc::AbstractFieldData, name::AbstractString, data) =
    throw(ArgumentError("$(nameof(typeof(loc))) is not supported for $(nameof(typeof(kind)))"))

validate_static_data(vtk::VTKHDFFile) = validate_data_totals(vtk, expected_totals(vtk, vtk.kind))

# `expected` maps group name => expected number of rows (checked for every
# array in that group); groups not listed are unconstrained.
function validate_data_totals(vtk::VTKHDFFile, expected::AbstractDict{String, Int})
    for (key, rows) in vtk.data_rows
        groupname = first(split(key, '/'))
        haskey(expected, groupname) || continue
        want = expected[groupname]
        rows == want || error("array $key has $rows tuples, expected $want")
    end
    return nothing
end

expected_totals(vtk, kind::DatasetKind) = Dict{String, Int}()
