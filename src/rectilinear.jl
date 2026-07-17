# RectilinearGrid writer (spec 2.7): like ImageData plus explicit
# X/Y/ZCoordinates datasets.

mutable struct RectilinearState <: StructuredKind
    pdims::NTuple{3, Int}
    coord_offsets::Vector{Int}   # current read offsets per axis
    new_coord_offsets::Vector{Union{Nothing, Int}}  # set when a step updates an axis
end

const COORD_NAMES = ("XCoordinates", "YCoordinates", "ZCoordinates")

function init_rectilinear(
        dest, x::AbstractVector{<:Real}, y::AbstractVector{<:Real},
        z::AbstractVector{<:Real};
        whole_extent = nothing, temporal = false, compress = false, chunk_size = 0
    )
    dims = (length(x), length(y), length(z))
    all(>=(1), dims) || throw(ArgumentError("coordinate vectors must be non-empty"))
    ext = check_whole_extent(whole_extent, dims)
    file, root = open_dest(dest)
    write_ascii_attribute(root, "Type", "RectilinearGrid")
    write_version_attribute(root, (2, 7))
    HDF5.attrs(root)["WholeExtent"] = Int64[ext...]
    vtk = make_vtkfile(
        file, root, RectilinearState(dims, [0, 0, 0], [nothing, nothing, nothing]);
        temporal, compress, chunk_size, version = (2, 7)
    )
    for (name, coords) in zip(COORD_NAMES, (x, y, z))
        c = Vector{Float64}(coords)
        if temporal
            append_rows(create_appendable(vtk, root, name, Float64, ()), c)
        else
            HDF5.write_dataset(root, name, c)
        end
    end
    return vtk
end

# ---- temporal hooks ----

function kind_begin_step!(vtk, kind::RectilinearState)
    fill!(kind.new_coord_offsets, nothing)
    return nothing
end

function supply_step_geometry!(
        vtk, kind::RectilinearState;
        x = nothing, y = nothing, z = nothing
    )
    for (i, coords) in enumerate((x, y, z))
        coords === nothing && continue
        length(coords) == kind.pdims[i] || throw(
            ArgumentError(
                "$(COORD_NAMES[i]) must keep length $(kind.pdims[i]), got $(length(coords))"
            )
        )
        ds = get_dataset(vtk.root, COORD_NAMES[i])
        kind.new_coord_offsets[i] = append_rows(ds, Vector{Float64}(coords))
    end
    return nothing
end

function finish_step_geometry!(vtk, kind::RectilinearState)
    for i in 1:3
        if kind.new_coord_offsets[i] !== nothing
            kind.coord_offsets[i] = kind.new_coord_offsets[i]
        end
    end
    return nothing
end

function append_step_offsets!(vtk, kind::RectilinearState, sg, geom)
    for i in 1:3
        append_rows(
            appendable(vtk, sg, COORD_NAMES[i] * "Offsets", Int64, ()),
            Int64(kind.coord_offsets[i])
        )
    end
    return nothing
end
