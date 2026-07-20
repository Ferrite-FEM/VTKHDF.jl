# Public constructors and dispatch.

const GRID_KWARGS_DOC = """
Common keyword arguments:

- `temporal::Bool = false`: open the file for time-dependent writing with
  [`write_timestep`](@ref).
- `compress = false`: `true` (gzip level 6) or an integer level `0-9`.
- `chunk_size = 0`: override the automatic HDF5 chunk length along the
  append dimension.

All forms support a leading function for the do-block idiom, which closes the
file afterwards. When `filename` has no extension, `.vtkhdf` is appended.
Instead of a filename, `(collection, blockname)` can be passed to create the
dataset as a block of a composite file (see [`vtkhdf_collection`](@ref)).
"""

"""
    vtkhdf_grid(filename, points, cells; kwargs...)
    vtkhdf_grid(filename, x, [y, [z]]; kwargs...)
    vtkhdf_grid(filename, xyz::AbstractArray{T,4}; kwargs...)
    vtkhdf_grid(dataset_type, filename, args...; kwargs...)

Create a VTKHDF file for one of the grid dataset types. The type is inferred
from the arguments (mirroring WriteVTK.jl):

- **UnstructuredGrid**: `points` (a `3×N` matrix — smaller first dimensions
  are zero-padded — a vector of point-like objects, or a tuple of coordinate
  vectors) and `cells`, a vector of `MeshCell`/`VTKPolyhedron`.
- **PolyData**: `points` plus one or more vectors of `MeshCell`s with
  `PolyData.*` cell types (e.g. `PolyData.Polys()`); each vector must contain
  a single category. On disk (and for cell data) cells are ordered Vertices,
  Lines, Polygons, Strips.
- **ImageData**: 1–3 `AbstractRange`s (uniform spacing), or
  `vtkhdf_grid(VTKImageData(), filename, dims; origin, spacing, direction,
  whole_extent)`.
- **RectilinearGrid**: 1–3 coordinate `AbstractVector`s, at least one of
  which is not a range (or the explicit `VTKRectilinearGrid()` tag).
- **StructuredGrid**: a `(3, ni, nj, nk)` coordinate array, or three
  `(ni, nj, nk)` arrays `x, y, z`.

Data arrays are then written with `vtk["name"] = data` (the location is
inferred from the size, or passed explicitly:
`vtk["name", VTKCellData()] = data`). Vector/tensor data uses the component-
first convention (`(3, N)` for vectors, like WriteVTK.jl), or a vector of
static vectors / tuples. The `attribute` keyword marks active attributes:
`vtk["u", VTKPointData(), attribute = :Vectors] = u`.

$GRID_KWARGS_DOC
"""
function vtkhdf_grid end

const Dest = Union{AbstractString, BlockDest}

# -- UnstructuredGrid --
vtkhdf_grid(dest::Dest, points, cells::UnstructuredCells; kwargs...) =
    register_and_init(dest, init_unstructured, points, cells; kwargs...)
vtkhdf_grid(::VTKUnstructuredGrid, dest::Dest; kwargs...) =
    register_and_init(dest, init_unstructured, nothing, nothing; kwargs...)
vtkhdf_grid(::VTKUnstructuredGrid, dest::Dest, points, cells::UnstructuredCells; kwargs...) =
    register_and_init(dest, init_unstructured, points, cells; kwargs...)

# -- PolyData --
vtkhdf_grid(
    dest::Dest, points, cells1::AbstractVector{<:PolyDataCell},
    cells::Vararg{AbstractVector{<:PolyDataCell}}; kwargs...
) =
    register_and_init(dest, init_polydata, points, (cells1, cells...); kwargs...)
vtkhdf_grid(::VTKPolyData, dest::Dest; kwargs...) =
    register_and_init(dest, init_polydata, nothing, (); kwargs...)
vtkhdf_grid(
    ::VTKPolyData, dest::Dest, points, cells1::AbstractVector{<:PolyDataCell},
    cells::Vararg{AbstractVector{<:PolyDataCell}}; kwargs...
) =
    register_and_init(dest, init_polydata, points, (cells1, cells...); kwargs...)

# -- abstractly-typed cell vectors (e.g. `MeshCell[...]`): classify at runtime --
const LooseCells = Union{
    AbstractVector{<:Union{MeshCell, VTKPolyhedron}},
    AbstractVector{Any},
}
function vtkhdf_grid(dest::Dest, points, cells::LooseCells; kwargs...)
    if !isempty(cells) && first(cells) isa MeshCell{<:PolyData.CellType}
        typed = PolyDataCell[c for c in cells]
        return register_and_init(dest, init_polydata, points, (typed,); kwargs...)
    end
    typed = Union{MeshCell{VTKCellType}, VTKPolyhedron}[c for c in cells]
    return register_and_init(dest, init_unstructured, points, typed; kwargs...)
end

# -- ImageData / RectilinearGrid from coordinate vectors --
function vtkhdf_grid(dest::Dest, coords::Vararg{AbstractVector{<:Real}, N}; kwargs...) where {N}
    1 <= N <= 3 || throw(ArgumentError("expected 1 to 3 coordinate vectors"))
    if all(c -> c isa AbstractRange, coords)
        dims3 = ntuple(i -> i <= N ? length(coords[i]) : 1, 3)
        origin = ntuple(i -> i <= N ? Float64(first(coords[i])) : 0.0, 3)
        spacing = ntuple(i -> i <= N && length(coords[i]) > 1 ? Float64(step(coords[i])) : 1.0, 3)
        return register_and_init(dest, init_image, dims3; origin, spacing, kwargs...)
    end
    x = coords[1]
    y = N >= 2 ? coords[2] : [0.0]
    z = N >= 3 ? coords[3] : [0.0]
    return register_and_init(dest, init_rectilinear, x, y, z; kwargs...)
end

function vtkhdf_grid(::VTKImageData, dest::Dest, dims::Tuple{Vararg{Integer}}; kwargs...)
    1 <= length(dims) <= 3 ||
        throw(ArgumentError("dimensions must have 1 to 3 entries, got $(length(dims))"))
    return register_and_init(dest, init_image, ntuple(i -> i <= length(dims) ? Int(dims[i]) : 1, 3); kwargs...)
end
vtkhdf_grid(
    ::VTKRectilinearGrid, dest::Dest, x::AbstractVector{<:Real},
    y::AbstractVector{<:Real} = [0.0], z::AbstractVector{<:Real} = [0.0]; kwargs...
) =
    register_and_init(dest, init_rectilinear, x, y, z; kwargs...)

# -- StructuredGrid --
function vtkhdf_grid(dest::Dest, xyz::AbstractArray{T, 4}; kwargs...) where {T <: Real}
    dims, pts = structured_points(xyz)
    return register_and_init(dest, init_structured, dims, pts; kwargs...)
end
function vtkhdf_grid(
        dest::Dest, x::AbstractArray{T, 3}, y::AbstractArray{T, 3},
        z::AbstractArray{T, 3}; kwargs...
    ) where {T <: Real}
    dims, pts = structured_points(x, y, z)
    return register_and_init(dest, init_structured, dims, pts; kwargs...)
end
function vtkhdf_grid(::VTKStructuredGrid, dest::Dest, xyz::AbstractArray{T, 4}; kwargs...) where {T <: Real}
    dims, pts = structured_points(xyz)
    return register_and_init(dest, init_structured, dims, pts; kwargs...)
end

register_and_init(dest, init, args...; kwargs...) =
    register_block(dest, init(dest, args...; kwargs...))

# Composite blocks: roll back the block group and its index when construction
# fails after the group was created (invalid options, bad geometry, ...), so a
# caller that catches the error is not left with an orphan block.
function register_and_init(dest::BlockDest, init, args...; kwargs...)
    col = dest.col
    col.isopen || error("collection is closed")
    check_name(dest.name, "block name")
    preexisting = haskey(col.root, dest.name)
    saved_index = col.next_index
    try
        return register_block(dest, init(dest, args...; kwargs...))
    catch
        if !preexisting && haskey(col.root, dest.name)
            HDF5.delete_object(col.root, dest.name)
            col.next_index = saved_index
        end
        rethrow()
    end
end

# tag-second form, so `vtkhdf_grid(col, "name", VTKImageData(), ...)` works
vtkhdf_grid(dest::Dest, tag::AbstractVTKDataset, args...; kwargs...) =
    vtkhdf_grid(tag, dest, args...; kwargs...)
# disambiguation against the (dest, points, cells...) methods; a dataset tag
# in the points slot is always a usage error
_tag_points_error(tag) = throw(
    ArgumentError(
        "dataset tag $(nameof(typeof(tag))) cannot be used as points; pass the tag before the filename or right after (collection, name)"
    )
)
vtkhdf_grid(dest::Dest, tag::AbstractVTKDataset, cells::UnstructuredCells; kwargs...) =
    _tag_points_error(tag)
vtkhdf_grid(dest::Dest, tag::AbstractVTKDataset, cells::LooseCells; kwargs...) =
    _tag_points_error(tag)
vtkhdf_grid(
    dest::Dest, tag::AbstractVTKDataset, cells1::AbstractVector{<:PolyDataCell},
    cells::Vararg{AbstractVector{<:PolyDataCell}}; kwargs...
) =
    _tag_points_error(tag)

# -- composite block forms: (collection, name, args...) --
for fn in (:vtkhdf_grid, :vtkhdf_table, :vtkhdf_amr, :vtkhdf_htg)
    @eval $fn(col::VTKHDFCollection, name::AbstractString, args...; kwargs...) =
        $fn(BlockDest(col, name), args...; kwargs...)
end
vtkhdf_table(dest::BlockDest; kwargs...) = register_and_init(dest, init_table; kwargs...)
vtkhdf_amr(dest::BlockDest; kwargs...) = register_and_init(dest, init_amr; kwargs...)
vtkhdf_htg(dest::BlockDest; kwargs...) = register_and_init(dest, init_htg; kwargs...)

# -- do-block forms --
for fn in (
        :vtkhdf_grid, :vtkhdf_table, :vtkhdf_amr, :vtkhdf_htg,
        :vtkhdf_collection, :vtkhdf_multiblock,
    )
    @eval function $fn(f::Function, args...; kwargs...)
        obj = $fn(args...; kwargs...)
        try
            f(obj)
        finally
            close(obj)
        end
        return obj
    end
end
