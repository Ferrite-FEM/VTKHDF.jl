# StructuredGrid writer (spec 2.7): like ImageData but with explicit point
# positions in a Points dataset of HDF shape (nz, ny, nx, 3) — a Julia
# (3, nx, ny, nz) array — with a prepended time dimension for temporal files.

mutable struct StructuredState <: StructuredKind
    pdims::NTuple{3, Int}
    point_slab::Int       # current read offset into the time dimension of Points
    new_point_slab::Union{Nothing, Int}
end

function structured_points(xyz::AbstractArray{T, 4}) where {T <: Real}
    size(xyz, 1) == 3 || throw(ArgumentError("points array must be (3, ni, nj, nk), got $(size(xyz))"))
    return (size(xyz)[2:4], xyz)
end

function structured_points(x::AbstractArray{T, 3}, y::AbstractArray{T, 3}, z::AbstractArray{T, 3}) where {T <: Real}
    size(x) == size(y) == size(z) || throw(ArgumentError("coordinate arrays must have equal size"))
    xyz = Array{T, 4}(undef, 3, size(x)...)
    xyz[1, :, :, :] .= x
    xyz[2, :, :, :] .= y
    xyz[3, :, :, :] .= z
    return (size(x), xyz)
end

function init_structured(
        dest, dims::NTuple{3, Int}, xyz::AbstractArray{<:Real, 4};
        whole_extent = nothing, temporal = false, compress = false, chunk_size = 0
    )
    file, root = open_dest(dest)
    write_ascii_attribute(root, "Type", "StructuredGrid")
    write_version_attribute(root, (2, 7))
    ext = whole_extent === nothing ?
        (0, dims[1] - 1, 0, dims[2] - 1, 0, dims[3] - 1) : Tuple(whole_extent)
    HDF5.attrs(root)["WholeExtent"] = Int64[ext...]
    vtk = make_vtkfile(
        file, root, StructuredState(dims, 0, nothing);
        temporal, compress, chunk_size, version = (2, 7)
    )
    T = eltype(xyz)
    if temporal
        ds = create_appendable(vtk, root, "Points", T, (3, dims...))
        append_rows(ds, reshape(xyz, (3, dims..., 1)))
    else
        HDF5.write_dataset(root, "Points", Array(xyz))
    end
    return vtk
end

# ---- temporal hooks ----

function kind_begin_step!(vtk, kind::StructuredState)
    kind.new_point_slab = nothing
    return nothing
end

function supply_step_geometry!(vtk, kind::StructuredState; points = nothing)
    points === nothing && return nothing
    xyz = points isa AbstractArray{<:Real, 4} ? points : structured_points(points...)[2]
    dims, xyz = structured_points(xyz)
    dims == kind.pdims || throw(ArgumentError("points must keep dimensions $(kind.pdims), got $dims"))
    kind.new_point_slab = append_rows(get_dataset(vtk.root, "Points"), reshape(xyz, (3, dims..., 1)))
    return nothing
end

function finish_step_geometry!(vtk, kind::StructuredState)
    if kind.new_point_slab !== nothing
        kind.point_slab = kind.new_point_slab
    end
    return nothing
end

function append_step_offsets!(vtk, kind::StructuredState, sg, geom)
    append_rows(appendable(vtk, sg, "PointOffsets", Int64, ()), Int64(kind.point_slab))
    return nothing
end
