# ImageData writer, plus the slab-style data machinery shared with
# RectilinearGrid and StructuredGrid.
#
# On disk, arrays of these kinds have HDF shape (nz, ny, nx[, ncomp]) which is
# exactly a Julia array ([ncomp,] nx, ny, nz) written unpermuted. Temporal
# files prepend a time dimension on disk, i.e. append along the last Julia
# dimension.

abstract type StructuredKind <: DatasetKind end

point_dims(kind::StructuredKind) = kind.pdims
cell_dims(kind::StructuredKind) = max.(point_dims(kind) .- 1, 1)

mutable struct ImageState <: StructuredKind
    pdims::NTuple{3, Int}
end

function init_image(
        dest, dims::NTuple{3, Integer};
        origin = (0.0, 0.0, 0.0), spacing = (1.0, 1.0, 1.0),
        direction = (1.0, 0, 0, 0, 1, 0, 0, 0, 1),
        whole_extent = nothing,
        temporal = false, compress = false, chunk_size = 0
    )
    all(>=(1), dims) || throw(ArgumentError("image dimensions must be positive, got $dims"))
    length(origin) == 3 || throw(ArgumentError("origin must have 3 entries"))
    length(spacing) == 3 || throw(ArgumentError("spacing must have 3 entries"))
    length(direction) == 9 || throw(ArgumentError("direction must have 9 entries (row-major 3×3)"))
    file, root = open_dest(dest)
    write_ascii_attribute(root, "Type", "ImageData")
    write_version_attribute(root, (2, 0))
    ext = whole_extent === nothing ?
        (0, dims[1] - 1, 0, dims[2] - 1, 0, dims[3] - 1) : Tuple(whole_extent)
    length(ext) == 6 || throw(ArgumentError("whole_extent must have 6 entries"))
    HDF5.attrs(root)["WholeExtent"] = Int64[ext...]
    HDF5.attrs(root)["Origin"] = Float64[origin...]
    HDF5.attrs(root)["Spacing"] = Float64[spacing...]
    HDF5.attrs(root)["Direction"] = Float64[direction...]
    return make_vtkfile(
        file, root, ImageState(Tuple(Int.(dims)));
        temporal, compress, chunk_size, version = (2, 0)
    )
end

# ---- shared slab data path ----

strip_trailing_ones(dims::Tuple) = Tuple(dims[1:something(findlast(!=(1), dims), 0)])

# Match an array shape against spatial dims; returns ncomp or nothing.
# `exact = true` restricts to component-free (scalar) matches.
function spatial_ncomp(ashape::Tuple, dims::NTuple{3, Int}; exact::Bool = false)
    sd = strip_trailing_ones(dims)
    strip_trailing_ones(ashape) == sd && return 1
    if !exact && length(ashape) >= 2 && strip_trailing_ones(ashape[2:end]) == sd
        return ashape[1]
    end
    return nothing
end

function resolve_location(vtk, kind::StructuredKind, data)
    data isa AbstractArray || throw(ArgumentError("expected an array"))
    # an exact (component-free) spatial match takes precedence over
    # interpreting the first dimension as components
    for exact in (true, false)
        p = spatial_ncomp(size(data), point_dims(kind); exact)
        c = spatial_ncomp(size(data), cell_dims(kind); exact)
        if p !== nothing && c !== nothing
            error(
                "data shape $(size(data)) matches both the point dimensions $(point_dims(kind)) " *
                    "and the cell dimensions $(cell_dims(kind)); pass VTKPointData() or VTKCellData() explicitly"
            )
        end
        p !== nothing && return VTKPointData()
        c !== nothing && return VTKCellData()
    end
    error(
        "data shape $(size(data)) matches neither the point dimensions $(point_dims(kind)) " *
            "nor the cell dimensions $(cell_dims(kind)); for field data pass VTKFieldData()"
    )
end

function write_array!(
        vtk::VTKHDFFile, kind::StructuredKind,
        loc::Union{VTKPointData, VTKCellData}, name::AbstractString, data::AbstractArray{T}
    ) where {T <: Real}
    dims = loc isa VTKPointData ? point_dims(kind) : cell_dims(kind)
    ncomp = spatial_ncomp(size(data), dims)
    ncomp === nothing && error(
        "data shape $(size(data)) does not match the $(
            loc isa VTKPointData ?
                "point" : "cell"
        ) dimensions $dims"
    )
    rowdims = ncomp == 1 ? dims : (ncomp, dims...)
    A = reshape(data, rowdims)
    groupname = location_group(loc)
    grp = get_or_create_group(vtk.root, groupname)
    key = groupname * "/" * name
    if vtk.temporal
        if haskey(grp, name)
            vtk.schema !== nothing && !(key in vtk.schema) &&
                error("array $key was not part of the first time step; the array schema is fixed by the first step")
            ds = get_dataset(grp, name)
        else
            vtk.schema !== nothing &&
                error("array $key was not part of the first time step; the array schema is fixed by the first step")
            ds = create_appendable(vtk, grp, name, T, rowdims)
        end
        append_rows(ds, reshape(A, (rowdims..., 1)))
        vtk.data_rows[key] = get(vtk.data_rows, key, 0) + 1
    else
        haskey(grp, name) && error("array $key already written")
        if vtk.compress > 0
            chunk = pick_chunk(T, rowdims, vtk.chunk_size)[1:length(rowdims)]
            ds = HDF5.create_dataset(
                grp, name, T, rowdims;
                chunk = Tuple(max.(chunk, 1)), compression_kwargs(vtk, T)...
            )
            HDF5.write_dataset(ds, HDF5.datatype(T), Array(A))
        else
            HDF5.write_dataset(grp, name, Array(A))
        end
        vtk.data_rows[key] = 1
    end
    return nothing
end

write_array!(vtk::VTKHDFFile, kind::StructuredKind, loc::VTKFieldData, name::AbstractString, data) =
    append_tuple_data!(vtk, "FieldData", name, data)

# steps: only Values (+ field data offsets) are needed; data arrays are
# indexed by the time dimension directly.
step_expected_totals(vtk, kind::StructuredKind, geom) = Dict(
    "PointData" => vtk.nsteps + 1,
    "CellData" => vtk.nsteps + 1,
)
