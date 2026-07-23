# Reader core: file opening, type detection, generic data-array reading and
# the temporal Steps cache. Kind-specific geometry reading lives in the other
# files of this directory.

abstract type ReaderKind end

# Temporal metadata, parsed eagerly at open (all Steps datasets are small).
# `datasets` maps the path relative to Steps ("PartOffsets",
# "PointDataOffsets/u", ...) to the dataset read as an Int64 array.
struct StepsInfo
    nsteps::Int
    values::Vector{Float64}
    datasets::Dict{String, Array{Int64}}
end

mutable struct VTKHDFReader{K <: ReaderKind}
    file::Union{HDF5.File, Nothing}  # nothing for composite blocks
    root::HDF5.Group
    kind::K
    version::NTuple{2, Int}
    steps::Union{Nothing, StepsInfo}
    parent::Any                      # ::Union{Nothing, VTKHDFCollectionReader}
    isopen::Bool
end

"""
    vtkhdf_open(filename) -> reader
    vtkhdf_open(f::Function, filename)

Open a VTKHDF file for reading. Returns a `VTKHDFReader` for simple datasets
or a `VTKHDFCollectionReader` for composite files
(`PartitionedDataSetCollection`/`MultiBlockDataSet`); the reader must be
`close`d (the do-block form does so automatically).

Data arrays are read by indexing (the location is found automatically, or
passed explicitly): `r["u"]`, `r["u", VTKCellData()]`. Geometry is read with
[`read_points`](@ref), [`read_cells`](@ref) and [`read_coordinates`](@ref).
Temporal files are read one step at a time via [`read_timestep`](@ref).

See the manual for the full reading API (much of which — `nsteps`,
`time_values`, `grid_info`, ... — is public but not exported).
"""
function vtkhdf_open(filename::AbstractString)
    if !ispath(filename)
        # mirror the writer: an extension-less name means <name>.vtkhdf
        with_ext = add_extension(filename)
        with_ext != filename && ispath(with_ext) || throw(ArgumentError("file $filename does not exist"))
        filename = with_ext
    end
    HDF5.ishdf5(String(filename)) ||
        throw(ArgumentError("$filename is not an HDF5 file"))
    file = h5open(filename, "r")
    try
        haskey(file, "VTKHDF") ||
            error("$filename is not a VTKHDF file: no /VTKHDF group")
        root = file["VTKHDF"]::HDF5.Group
        return open_reader(file, root, nothing)
    catch
        close(file)
        rethrow()
    end
end

function vtkhdf_open(f::Function, filename::AbstractString)
    r = vtkhdf_open(filename)
    try
        return f(r)
    finally
        close(r)
    end
end

# Build the appropriate reader for a /VTKHDF-style group (the file root, or a
# composite block with `parent` set). Composite types are handled by the
# caller (vtkhdf_open) via open_collection.
function open_reader(file, root::HDF5.Group, parent)
    type = read_type_attribute(root)
    version = read_version_attribute(root)
    if type in ("PartitionedDataSetCollection", "MultiBlockDataSet")
        parent === nothing ||
            error("nested composite datasets are not supported (block $(HDF5.name(root)))")
        return open_collection(file, root, type)
    end
    kind = make_reader_kind(type, root)
    steps = read_steps_info(root, kind)
    steps === nothing || validate_reader(kind, root, steps)
    return VTKHDFReader(file, root, kind, version, steps, parent, true)
end

# Kind-specific consistency checks run at open for temporal files.
validate_reader(kind::ReaderKind, root, steps) = nothing

function make_reader_kind(type::String, root::HDF5.Group)
    return type == "UnstructuredGrid" ? read_unstructured_kind(root) :
        type == "PolyData" ? read_polydata_kind(root) :
        type == "ImageData" ? read_image_kind(root) :
        type == "RectilinearGrid" ? read_rectilinear_kind(root) :
        type == "StructuredGrid" ? read_structured_kind(root) :
        type == "Table" ? read_table_kind(root) :
        type == "OverlappingAMR" ? read_amr_kind(root) :
        type == "HyperTreeGrid" ? read_htg_kind(root) :
        error("unsupported VTKHDF dataset type $(repr(type))")
end

function read_type_attribute(root)
    haskey(HDF5.attrs(root), "Type") ||
        error("$(HDF5.name(root)) has no Type attribute; not a valid VTKHDF group")
    type = HDF5.read_attribute(root, "Type")
    type isa AbstractString ||
        error("$(HDF5.name(root)) has a malformed Type attribute")
    return String(rstrip(type, '\0'))
end

function read_version_attribute(root)
    haskey(HDF5.attrs(root), "Version") ||
        error("$(HDF5.name(root)) has no Version attribute; not a valid VTKHDF group")
    v = HDF5.read_attribute(root, "Version")
    (v isa AbstractVector{<:Integer} && length(v) >= 2) ||
        error("$(HDF5.name(root)) has a malformed Version attribute")
    version = (Int(v[1]), Int(v[2]))
    1 <= version[1] <= 2 || error(
        "unsupported VTKHDF version $(version[1]).$(version[2]) in $(HDF5.name(root)); " *
            "this reader supports major versions 1 and 2"
    )
    return version
end

# ---- lifecycle ----

function check_open(r::VTKHDFReader)
    if r.parent !== nothing
        r.parent.isopen || error("the collection this block belongs to is closed")
    end
    r.isopen || error("reader is closed")
    return nothing
end

function Base.close(r::VTKHDFReader)
    r.isopen || return nothing
    r.isopen = false
    if r.parent === nothing
        close(r.root)
        file = r.file
        file === nothing || close(file)
    end
    return nothing
end

Base.isopen(r::VTKHDFReader) = r.isopen && (r.parent === nothing || r.parent.isopen)

function Base.show(io::IO, r::VTKHDFReader)
    status = isopen(r) ? "open" : "closed"
    print(io, "VTKHDFReader{", nameof(typeof(r.kind)), "} (", status)
    r.steps === nothing || print(io, ", temporal, ", r.steps.nsteps, " steps")
    return print(io, ")")
end

"""
    dataset_type(r) -> String

The VTKHDF dataset type of an open reader (the `Type` attribute), e.g.
`"UnstructuredGrid"`. For composite files, `"PartitionedDataSetCollection"`
or `"MultiBlockDataSet"`; for an empty (type-less) composite block,
`nothing`.
"""
dataset_type(r::VTKHDFReader) = (check_open(r); reader_type_string(r.kind))

"""
    file_version(r) -> (major, minor)

The VTKHDF specification version of the file (the `Version` attribute).
"""
file_version(r::VTKHDFReader) = (check_open(r); r.version)

# ---- temporal metadata ----

"""
    is_temporal(r) -> Bool

Whether the file contains time-dependent data (a `Steps` group). Temporal
data is read step by step with [`read_timestep`](@ref).
"""
is_temporal(r::VTKHDFReader) = r.steps !== nothing

# The StepsInfo of a temporal reader (throws for static files); used
# internally so field access is well-typed.
function steps_info(r::VTKHDFReader)
    si = r.steps
    return si === nothing ? error("not a temporal file") : si
end

"""
    nsteps(r) -> Int

Number of time steps in a temporal file ([`is_temporal`](@ref)).
"""
nsteps(r::VTKHDFReader) = (check_open(r); steps_info(r).nsteps)

"""
    time_values(r) -> Vector{Float64}

The time value of every step of a temporal file.
"""
time_values(r::VTKHDFReader) = (check_open(r); copy(steps_info(r).values))

function read_steps_info(root::HDF5.Group, kind::ReaderKind)
    haskey(root, "Steps") || return nothing
    supports_temporal(kind) || error(
        "temporal $(reader_type_string(kind)) reading is not supported"
    )
    sg = root["Steps"]::HDF5.Group
    haskey(HDF5.attrs(sg), "NSteps") || error("Steps group has no NSteps attribute")
    n = Int(HDF5.read_attribute(sg, "NSteps"))
    n >= 0 || error("Steps has negative NSteps ($n)")
    haskey(sg, "Values") || error("Steps group has no Values dataset")
    values = Float64.(read(get_dataset(sg, "Values")))
    length(values) == n ||
        error("Steps/Values has $(length(values)) entries, expected NSteps = $n")
    datasets = Dict{String, Array{Int64}}()
    for name in keys(sg)
        name == "Values" && continue
        obj = sg[name]
        if obj isa HDF5.Group
            for sub in keys(obj)
                store_steps_dataset!(datasets, name * "/" * sub, obj[sub], n)
            end
        else
            store_steps_dataset!(datasets, name, obj, n)
        end
    end
    return StepsInfo(n, values, datasets)
end

function store_steps_dataset!(datasets, key, ds, nsteps)
    ds isa HDF5.Dataset || error("Steps/$key is not a dataset")
    A = Int64.(read(ds))
    size(A, ndims(A)) == nsteps ||
        error("Steps/$key has $(size(A, ndims(A))) entries, expected NSteps = $nsteps")
    datasets[key] = A
    return nothing
end

# Scalar per-step entry with a default; `key` is e.g. "PointOffsets".
function steps_entry(si::StepsInfo, key::String, i::Int, default::Int)
    A = get(si.datasets, key, nothing)
    A === nothing && return default
    A isa Vector || error("Steps/$key has unexpected rank $(ndims(A)); expected one entry per step")
    return Int(A[i])
end

# Per-topology column of a (NTopologies, NSteps) offsets dataset (a 1-D
# dataset is accepted as the single-topology case).
function steps_column(si::StepsInfo, key::String, i::Int, ntopo::Int)
    A = get(si.datasets, key, nothing)
    A === nothing && return nothing
    if A isa Vector
        ntopo == 1 || error("Steps/$key is 1-D but $ntopo topology columns were expected")
        return Int[A[i]]
    end
    size(A, 1) == ntopo ||
        error("Steps/$key has $(size(A, 1)) topology columns, expected $ntopo")
    return Int.(A[:, i])
end

steps_dataset(si::StepsInfo, key::String) = get(si.datasets, key, nothing)

function check_step_index(r::VTKHDFReader, i::Integer)
    check_open(r)
    r.steps === nothing &&
        error("not a temporal file; read data and geometry directly from the reader")
    si = steps_info(r)
    1 <= i <= si.nsteps || throw(BoundsError(si.values, i))
    return si
end

struct VTKHDFTimeStep{R <: VTKHDFReader}
    reader::R
    index::Int
end

"""
    read_timestep(r, i) -> step

A view of time step `i` (1-based) of a temporal file. The step supports the
same data and geometry access as a static reader — `step["u"]`,
`read_points(step)`, `read_cells(step)`, ... — plus [`time_value`](@ref).
It holds no resources of its own (the parent reader must stay open).
"""
function read_timestep(r::VTKHDFReader, i::Integer)
    check_step_index(r, i)
    return VTKHDFTimeStep(r, Int(i))
end

"""
    time_value(step) -> Float64

The time value of a step obtained from [`read_timestep`](@ref).
"""
time_value(s::VTKHDFTimeStep) = (check_open(s.reader); steps_info(s.reader).values[s.index])

function Base.show(io::IO, s::VTKHDFTimeStep)
    si = steps_info(s.reader)
    return print(
        io, "VTKHDFTimeStep ", s.index, "/", si.nsteps, " (t = ", si.values[s.index], ")"
    )
end

# ---- data arrays ----

# Data groups a kind can hold, in auto-location search order.
data_locations(::ReaderKind) = (VTKPointData(), VTKCellData(), VTKFieldData())

# location_group with a concrete return type for abstractly-typed locations.
read_group_name(loc::AbstractFieldData) = location_group(loc)::String

"""
    keys(r, loc) -> Vector{String}

Names of the data arrays stored at `loc` (`VTKPointData()`, `VTKCellData()`,
`VTKFieldData()`, `VTKRowData()`); works on readers and time steps.
"""
function Base.keys(r::VTKHDFReader, loc::AbstractFieldData)
    check_open(r)
    grp = read_group_name(loc)
    haskey(r.root, grp) || return String[]
    return sort!(collect(String, keys(r.root[grp]::HDF5.Group)))
end

Base.haskey(r::VTKHDFReader, name::AbstractString, loc::AbstractFieldData) =
    (check_open(r); haskey(r.root, read_group_name(loc) * "/" * (String(name)::String)))

function resolve_read_location(r::VTKHDFReader, name::AbstractString)
    hits = AbstractFieldData[]
    for loc in data_locations(r.kind)
        haskey(r, name, loc) && push!(hits, loc)
    end
    isempty(hits) && error(
        "no data array named $(repr(String(name))) in " *
            join((location_group(l) for l in data_locations(r.kind)), ", ")
    )
    length(hits) > 1 && error(
        "data array name $(repr(String(name))) is ambiguous (present in " *
            join((location_group(l) for l in hits), " and ") *
            "); pass the location explicitly"
    )
    return hits[1]
end

"""
    r["name"]
    r["name", loc]

Read a data array by name from a reader or a time step. Without a location
the dataset's locations are searched (an ambiguous name throws); `loc`
reads from that location only. Values come back as plain arrays in the
writing shape convention (see [Data arrays](@ref data-arrays)).
"""
function Base.getindex(r::VTKHDFReader, name::AbstractString)
    check_open(r)
    return r[name, resolve_read_location(r, name)]
end

function Base.getindex(r::VTKHDFReader, name::AbstractString, loc::AbstractFieldData)
    check_open(r)
    is_temporal(r) &&
        error("temporal file: read data through read_timestep, e.g. read_timestep(r, i)[$(repr(String(name)))]")
    haskey(r, name, loc) ||
        error("no data array named $(repr(String(name))) in $(read_group_name(loc))")
    ds = get_dataset(r.root[read_group_name(loc)]::HDF5.Group, String(name))
    return read_static_array(r.kind, loc, ds)
end

function Base.getindex(s::VTKHDFTimeStep, name::AbstractString)
    check_open(s.reader)
    return s[name, resolve_read_location(s.reader, name)]
end

function Base.getindex(s::VTKHDFTimeStep, name::AbstractString, loc::AbstractFieldData)
    r = s.reader
    check_step_index(r, s.index)
    haskey(r, name, loc) ||
        error("no data array named $(repr(String(name))) in $(read_group_name(loc))")
    ds = get_dataset(r.root[read_group_name(loc)]::HDF5.Group, String(name))
    loc isa VTKFieldData && return read_step_field_array(r, String(name), ds, s.index)
    return read_step_array(r, r.kind, loc, String(name), ds, s.index)
end

Base.keys(s::VTKHDFTimeStep, loc::AbstractFieldData) = keys(s.reader, loc)
Base.haskey(s::VTKHDFTimeStep, name::AbstractString, loc::AbstractFieldData) =
    haskey(s.reader, name, loc)

# Static read of a whole dataset in the tuple layout used by
# unstructured-like kinds: 1-D => Vector, 2-D (ncomp, N) => Matrix; strings
# (FieldData) => Vector{String}.
function read_static_array(kind::ReaderKind, loc::AbstractFieldData, ds::HDF5.Dataset)
    is_string_dataset(ds) && return read_string_array(ds)
    nd = ndims(ds)
    nd == 1 || nd == 2 ||
        error("unexpected rank $nd of data array $(HDF5.name(ds))")
    return read(ds)
end

function is_string_dataset(ds::HDF5.Dataset)
    dt = HDF5.datatype(ds)
    try
        return HDF5.API.h5t_get_class(dt.id) == HDF5.API.H5T_STRING
    finally
        close(dt)
    end
end

read_string_array(ds::HDF5.Dataset) = String.(read(ds))

# Tuple-slab read of rows `range` (1-based) of an appendable dataset.
function read_tuple_rows(ds::HDF5.Dataset, range::UnitRange{Int})
    nd = ndims(ds)
    n = size(ds, nd)
    (first(range) >= 1 && last(range) <= n) || error(
        "computed slice $range is out of bounds for $(HDF5.name(ds)) with $n rows; " *
            "the file's offset information is inconsistent"
    )
    isempty(range) && return nd == 1 ? read(ds)[1:0] : read(ds)[:, 1:0]
    nd == 1 && return ds[range]
    nd == 2 && return ds[:, range]
    error("unexpected rank $nd of data array $(HDF5.name(ds))")
end

# Per-step FieldData: offsets from FieldDataOffsets, tuple counts and ncomp
# from FieldDataSizes; spec defaults when either is absent.
function read_step_field_array(r::VTKHDFReader, name::String, ds::HDF5.Dataset, i::Int)
    is_string_dataset(ds) &&
        error("string field data is not supported in temporal files")
    si = steps_info(r)
    sizes = steps_dataset(si, "FieldDataSizes/" * name)
    if sizes === nothing
        ntuples = 1
        offset = steps_entry(si, "FieldDataOffsets/" * name, i, i - 1)
    else
        (sizes isa Matrix && size(sizes, 1) == 2) ||
            error("Steps/FieldDataSizes/$name must be (2, NSteps): (ncomponents, ntuples) per step")
        ntuples = Int(sizes[2, i])
        # default offset: cumulative tuple count of the preceding steps
        default = Int(sum(view(sizes, 2, 1:(i - 1)); init = Int64(0)))
        offset = steps_entry(si, "FieldDataOffsets/" * name, i, default)
    end
    A = read_tuple_rows(ds, (offset + 1):(offset + ntuples))
    if sizes !== nothing && A isa Matrix
        # the dataset holds the maximum component count over all steps; this
        # step may use fewer components
        ncomp = Int(sizes[1, i])
        1 <= ncomp <= size(A, 1) || error(
            "Steps/FieldDataSizes/$name declares $ncomp components for step $i, " *
                "but the dataset has $(size(A, 1))"
        )
        ncomp == size(A, 1) || (A = ncomp == 1 ? vec(A[1:1, :]) : A[1:ncomp, :])
    end
    return A
end

# ---- active attributes ----

"""
    active_attributes(r, loc) -> Dict{Symbol, String}

The active-attribute roles (`:Scalars`, `:Vectors`, ...) defined for the
data-array group at `loc`, mapping role to array name. Both spec encodings
are understood: group attributes (`PointData` attribute `"Vectors" => "u"`)
take precedence; per-array `Attribute` attributes (spec 2.6, matched
case-insensitively) fill in roles the group does not define.
"""
function active_attributes(r::Union{VTKHDFReader, VTKHDFTimeStep}, loc::AbstractFieldData)
    rr = r isa VTKHDFTimeStep ? r.reader : r
    check_open(rr)
    out = Dict{Symbol, String}()
    haskey(rr.root, read_group_name(loc)) || return out
    grp = rr.root[read_group_name(loc)]::HDF5.Group
    names = collect(String, keys(grp))
    gattrs = HDF5.attrs(grp)
    for role in ATTRIBUTE_KINDS
        haskey(gattrs, String(role)) || continue
        v = gattrs[String(role)]
        v isa AbstractString || continue
        out[role] = rstrip(String(v), '\0')
    end
    # per-array Attribute values fill roles the group attributes do not define
    for name in sort(names)
        ds = grp[name]
        ds isa HDF5.Dataset || continue
        haskey(HDF5.attrs(ds), "Attribute") || continue
        v = HDF5.attrs(ds)["Attribute"]
        v isa AbstractString || continue
        role = findfirst(k -> lowercase(String(k)) == lowercase(rstrip(String(v), '\0')), ATTRIBUTE_KINDS)
        role === nothing && continue  # unknown role names are ignored
        rolesym = ATTRIBUTE_KINDS[role]
        if haskey(out, rolesym)
            out[rolesym] == name || @warn(
                "conflicting active-attribute encodings for $rolesym in $(read_group_name(loc)): " *
                    "group attribute names $(repr(out[rolesym])), array $(repr(name)) claims the " *
                    "role too; using the group attribute"
            )
        else
            out[rolesym] = name
        end
    end
    return out
end

"""
    data_attributes(r, name, loc) -> Vector{Symbol}

All active-attribute roles resolved to the array `name` at `loc`
(see [`active_attributes`](@ref)).
"""
function data_attributes(r::Union{VTKHDFReader, VTKHDFTimeStep}, name::AbstractString, loc::AbstractFieldData)
    active = active_attributes(r, loc)
    return sort!([role for (role, n) in active if n == String(name)])
end

# ---- shared helpers for kind implementations ----

# 1-based prefix sums of a per-partition count dataset: prefix[p] is the
# total before partition p, prefix[end] the grand total.
function count_prefix(root::HDF5.Group, name::String)
    haskey(root, name) || error("missing dataset $name in $(HDF5.name(root))")
    counts = read(get_dataset(root, name))
    all(>=(0), counts) || error("dataset $name contains negative counts")
    prefix = Vector{Int}(undef, length(counts) + 1)
    prefix[1] = 0
    for (i, c) in enumerate(counts)
        prefix[i + 1] = prefix[i] + Int(c)
    end
    return prefix
end

npartitions_prefix(prefix::Vector{Int}) = length(prefix) - 1

function require_dataset(root::HDF5.Group, name::String)
    haskey(root, name) ||
        error("missing dataset $name in $(HDF5.name(root)); not a valid VTKHDF layout")
    return get_dataset(root, name)
end

function read_attr_tuple(root, name, ::Type{T}, ::Val{N}) where {T, N}
    haskey(HDF5.attrs(root), name) || return nothing
    v = HDF5.read_attribute(root, name)
    length(v) == N || error("attribute $name has $(length(v)) entries, expected $N")
    return ntuple(i -> T(v[i]), Val(N))::NTuple{N, T}
end

# Kinds define: reader_type_string, supports_temporal, data_locations,
# and geometry accessors. Defaults:
supports_temporal(::ReaderKind) = true

"""
    grid_info(r) -> NamedTuple

The metadata attributes stored at the dataset root, per dataset type:

- ImageData: `(; dims, origin, spacing, direction, whole_extent)`
- RectilinearGrid, StructuredGrid: `(; dims, whole_extent)`
- OverlappingAMR: `(; origin, grid_description)`
- HyperTreeGrid: `(; dimensions, branch_factor, transposed_root_indexing)`

UnstructuredGrid and PolyData store no such attributes and have no
`grid_info` method.
"""
function grid_info end

"""
    read_points(r_or_step)

Point coordinates: a `3 × N` matrix for UnstructuredGrid/PolyData, a
`(3, ni, nj, nk)` array for StructuredGrid.
"""
function read_points end

"""
    read_cells(r_or_step)

The cells of an UnstructuredGrid (a vector of `MeshCell`/`VTKPolyhedron`
with 1-based connectivity into [`read_points`](@ref)), or of a PolyData
(a NamedTuple `(; vertices, lines, polygons, strips)` of `MeshCell`
vectors). Multi-partition files are concatenated with point ids rebased;
see [`partition_ranges`](@ref) for the partition structure.
"""
function read_cells end

"""
    read_coordinates(r_or_step) -> (x, y, z)

The coordinate vectors of a RectilinearGrid.
"""
function read_coordinates end

"""
    npoints(r_or_step) -> Int

Total number of points (UnstructuredGrid, PolyData, ImageData,
RectilinearGrid, StructuredGrid; also [`amr_level`](@ref) handles).
"""
function npoints end

"""
    ncells(r_or_step) -> Int

Total number of cells (UnstructuredGrid, PolyData, ImageData,
RectilinearGrid, StructuredGrid, HyperTreeGrid; also [`amr_level`](@ref)
handles).
"""
function ncells end

"""
    npartitions(r_or_step) -> Int

Number of geometry partitions (UnstructuredGrid, PolyData).
"""
function npartitions end

"""
    partition_ranges(r_or_step) -> NamedTuple

The index ranges of each partition into the concatenated point/cell data
arrays. For UnstructuredGrid: `(; points, cells)`, each a vector with one
range per partition. For PolyData additionally `cells_by_category`: per
partition, a `(; vertices, lines, polygons, strips)` NamedTuple of ranges
into the cell-data arrays (on disk, PolyData cell data is partition-major
with the categories in that order inside each partition). For an
[`amr_level`](@ref) handle: `(; points, cells)` with one range per AMR box.
"""
function partition_ranges end

function check_static_geometry(r::VTKHDFReader, what::String)
    check_open(r)
    is_temporal(r) &&
        error("temporal file: read $what through read_timestep, e.g. $what(read_timestep(r, i))")
    return nothing
end
