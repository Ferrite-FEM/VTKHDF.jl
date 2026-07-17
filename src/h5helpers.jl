# Low-level HDF5 helpers shared by all writers.

const CHUNK_TARGET_BYTES = 1 << 20  # ~1 MiB chunks

function check_name(name::AbstractString, what::String = "array name")
    isempty(name) && throw(ArgumentError("$what must not be empty"))
    isascii(name) || throw(ArgumentError("$what must be ASCII, got $(repr(name))"))
    if occursin('/', name) || occursin('.', name)
        throw(ArgumentError("VTKHDF does not allow '/' or '.' in names, got $(repr(name))"))
    end
    occursin('\0', name) && throw(ArgumentError("$what must not contain NUL bytes"))
    return name
end

# Fixed-length ASCII scalar string attribute (VTK's readers require this
# encoding for e.g. the "Type" attribute; variable-length UTF-8 breaks them).
function write_ascii_attribute(parent, name::AbstractString, value::String)
    isascii(value) || throw(ArgumentError("attribute $name must be ASCII, got $(repr(value))"))
    occursin('\0', value) && throw(ArgumentError("attribute $name must not contain NUL bytes"))
    dtype_id = HDF5.API.h5t_copy(HDF5.API.H5T_C_S1)
    HDF5.API.h5t_set_size(dtype_id, ncodeunits(value))
    HDF5.API.h5t_set_cset(dtype_id, HDF5.API.H5T_CSET_ASCII)
    HDF5.API.h5t_set_strpad(dtype_id, HDF5.API.H5T_STR_NULLPAD)
    dtype = HDF5.Datatype(dtype_id)
    dspace = HDF5.Dataspace(HDF5.API.h5s_create(HDF5.API.H5S_SCALAR))
    haskey(HDF5.attrs(parent), name) && HDF5.delete_attribute(parent, name)
    attr = HDF5.create_attribute(parent, name, dtype, dspace)
    try
        GC.@preserve value HDF5.API.h5a_write(attr.id, dtype_id, pointer(value))
    finally
        close(attr)
        close(dspace)
        close(dtype)
    end
    return nothing
end

# Variable-length ASCII string dataset (VTK cannot read HDF5.jl's default
# variable-length UTF-8 strings).
function ascii_string_datatype()
    dtype_id = HDF5.API.h5t_copy(HDF5.API.H5T_C_S1)
    HDF5.API.h5t_set_size(dtype_id, HDF5.API.H5T_VARIABLE)
    HDF5.API.h5t_set_cset(dtype_id, HDF5.API.H5T_CSET_ASCII)
    return HDF5.Datatype(dtype_id)
end

function write_string_dataset(parent, name::AbstractString, strs::AbstractVector{String})
    for s in strs
        isascii(s) || throw(ArgumentError("string data must be ASCII, got $(repr(s))"))
    end
    dtype = ascii_string_datatype()
    dspace = HDF5.dataspace((length(strs),))
    ds = HDF5.create_dataset(parent, name, dtype, dspace)
    try
        ptrs = Base.unsafe_convert.(Cstring, strs)
        GC.@preserve strs HDF5.API.h5d_write(
            ds.id, dtype.id, HDF5.API.H5S_ALL, HDF5.API.H5S_ALL, HDF5.API.H5P_DEFAULT, ptrs,
        )
    finally
        close(ds)
        close(dspace)
        close(dtype)
    end
    return nothing
end

function create_soft_link(parent::HDF5.Group, name::AbstractString, target::AbstractString)
    HDF5.API.h5l_create_soft(target, parent.id, name, HDF5.API.H5P_DEFAULT, HDF5.API.H5P_DEFAULT)
    return nothing
end

# Chunk dimensions for a dataset with Julia dims (rowdims..., append-dim).
# Aims for CHUNK_TARGET_BYTES, tiling the row dims when a single row is
# already larger than the target (e.g. big ImageData time slices).
function pick_chunk(::Type{T}, rowdims::Dims, chunk_size::Int) where {T}
    rowbytes = max(sizeof(T), 1) * max(prod(rowdims), 1)
    if rowbytes >= CHUNK_TARGET_BYTES
        # tile row dims greedily (first Julia dims kept whole while they fit)
        remaining = max(CHUNK_TARGET_BYTES ÷ max(sizeof(T), 1), 1)
        tile = Int[]
        for n in rowdims
            if remaining >= n
                push!(tile, n)
                remaining = max(remaining ÷ max(n, 1), 1)
            else
                push!(tile, max(remaining, 1))
                remaining = 1
            end
        end
        return (tile..., 1)
    end
    nrows = chunk_size > 0 ? chunk_size : clamp(CHUNK_TARGET_BYTES ÷ rowbytes, 1, 1 << 16)
    return (rowdims..., Int(nrows))
end

compression_kwargs(vtk, ::Type{T}) where {T} =
    (vtk.compress > 0 && isbitstype(T)) ? (shuffle = true, deflate = vtk.compress) : NamedTuple()

# Appendable dataset: Julia dims (rowdims..., N) with unlimited last dimension
# (on disk: first HDF dimension unlimited, matching VTK's own files).
function create_appendable(vtk, parent, name::AbstractString, ::Type{T}, rowdims::Dims) where {T}
    dims = (rowdims..., 0)
    maxdims = (rowdims..., -1)
    chunk = pick_chunk(T, rowdims, vtk.chunk_size)
    return HDF5.create_dataset(
        parent, name, T, (dims, maxdims);
        chunk, compression_kwargs(vtk, T)...
    )
end

# Append `n` rows to an appendable dataset; `data`'s last dimension is the row
# dimension (a Vector for scalar rows).
function append_rows(ds::HDF5.Dataset, data::AbstractArray)
    nd = ndims(ds)
    ndims(data) == nd || throw(ArgumentError("expected $(nd)-dimensional data for $(HDF5.name(ds))"))
    old = size(ds)
    n = size(data, nd)
    if nd > 1 && size(data)[1:(nd - 1)] != old[1:(nd - 1)]
        throw(ArgumentError("row dimensions $(size(data)[1:(nd - 1)]) do not match dataset $(HDF5.name(ds)) $(old[1:(nd - 1)])"))
    end
    HDF5.set_extent_dims(ds, (old[1:(nd - 1)]..., old[nd] + n))
    if n > 0
        idx = ntuple(i -> i == nd ? ((old[nd] + 1):(old[nd] + n)) : Colon(), nd)
        ds[idx...] = data
    end
    return old[nd]  # previous number of rows (the offset of the appended data)
end

append_rows(ds::HDF5.Dataset, data::Number) = append_rows(ds, [data])

n_rows(ds::HDF5.Dataset) = size(ds, ndims(ds))

# Get-or-create an appendable dataset under `parent`.
function appendable(vtk, parent, name::AbstractString, ::Type{T}, rowdims::Dims) where {T}
    return haskey(parent, name) ? parent[name]::HDF5.Dataset :
        create_appendable(vtk, parent, name, T, rowdims)
end

get_or_create_group(parent, name::AbstractString; track_order::Bool = false) =
    haskey(parent, name) ? parent[name]::HDF5.Group :
    HDF5.create_group(parent, name; track_order)

get_dataset(parent, name::AbstractString) = parent[name]::HDF5.Dataset

function write_version_attribute(root, version::NTuple{2, Int})
    HDF5.attrs(root)["Version"] = Int64[version...]
    return nothing
end
