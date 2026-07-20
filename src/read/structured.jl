# ImageData / RectilinearGrid / StructuredGrid readers.
#
# Point/cell arrays of these kinds are stored whole, HDF shape
# (nz, ny, nx[, ncomp]) == Julia ([ncomp,] nx, ny, nz); temporal files append
# a trailing Julia time dimension, so step i is a plain index there.

abstract type ReadStructuredKind <: ReaderKind end

struct ReadImage <: ReadStructuredKind
    pdims::NTuple{3, Int}
    origin::NTuple{3, Float64}
    spacing::NTuple{3, Float64}
    direction::NTuple{9, Float64}
    whole_extent::NTuple{6, Int}
end

struct ReadRectilinear <: ReadStructuredKind
    pdims::NTuple{3, Int}
    whole_extent::NTuple{6, Int}
end

struct ReadStructured <: ReadStructuredKind
    pdims::NTuple{3, Int}
    whole_extent::NTuple{6, Int}
end

reader_type_string(::ReadImage) = "ImageData"
reader_type_string(::ReadRectilinear) = "RectilinearGrid"
reader_type_string(::ReadStructured) = "StructuredGrid"

function read_whole_extent(root)
    ext = read_attr_tuple(root, "WholeExtent", Int, Val(6))
    ext === nothing && error("$(HDF5.name(root)) has no WholeExtent attribute")
    dims = ntuple(i -> ext[2i] - ext[2i - 1] + 1, 3)
    all(>=(1), dims) || error("invalid WholeExtent $ext")
    return ext, dims
end

function read_image_kind(root::HDF5.Group)
    ext, dims = read_whole_extent(root)
    origin = something(read_attr_tuple(root, "Origin", Float64, Val(3)), (0.0, 0.0, 0.0))
    spacing = something(read_attr_tuple(root, "Spacing", Float64, Val(3)), (1.0, 1.0, 1.0))
    direction = something(
        read_attr_tuple(root, "Direction", Float64, Val(9)),
        (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
    )
    return ReadImage(dims, origin, spacing, direction, ext)
end

function read_rectilinear_kind(root::HDF5.Group)
    ext, dims = read_whole_extent(root)
    return ReadRectilinear(dims, ext)
end

function read_structured_kind(root::HDF5.Group)
    ext, dims = read_whole_extent(root)
    return ReadStructured(dims, ext)
end

point_dims(k::ReadStructuredKind) = k.pdims
cell_dims(k::ReadStructuredKind) = max.(k.pdims .- 1, 1)

grid_info(r::VTKHDFReader{ReadImage}) = (
    dims = r.kind.pdims, origin = r.kind.origin, spacing = r.kind.spacing,
    direction = r.kind.direction, whole_extent = r.kind.whole_extent,
)
grid_info(r::VTKHDFReader{ReadRectilinear}) =
    (dims = r.kind.pdims, whole_extent = r.kind.whole_extent)
grid_info(r::VTKHDFReader{ReadStructured}) =
    (dims = r.kind.pdims, whole_extent = r.kind.whole_extent)

npoints(r::VTKHDFReader{<:ReadStructuredKind}) = (check_open(r); prod(point_dims(r.kind)))
ncells(r::VTKHDFReader{<:ReadStructuredKind}) = (check_open(r); prod(cell_dims(r.kind)))
npoints(s::VTKHDFTimeStep{<:VTKHDFReader{<:ReadStructuredKind}}) = npoints(s.reader)
ncells(s::VTKHDFTimeStep{<:VTKHDFReader{<:ReadStructuredKind}}) = ncells(s.reader)

# Temporal point/cell arrays must have one trailing slab per step.
function validate_reader(kind::ReadStructuredKind, root::HDF5.Group, steps::StepsInfo)
    for groupname in ("PointData", "CellData")
        haskey(root, groupname) || continue
        grp = root[groupname]::HDF5.Group
        for name in keys(grp)
            ds = grp[name]
            ds isa HDF5.Dataset || continue
            n = size(ds, ndims(ds))
            n == steps.nsteps || error(
                "temporal array $groupname/$name has $n time slabs, expected NSteps = $(steps.nsteps)"
            )
        end
    end
    return nothing
end

# ---- data arrays ----

function read_static_array(kind::ReadStructuredKind, loc::Union{VTKPointData, VTKCellData}, ds::HDF5.Dataset)
    A = read(ds)
    check_spatial_shape(kind, loc, size(A), HDF5.name(ds))
    return A
end

function read_step_array(
        r::VTKHDFReader{<:ReadStructuredKind}, kind::ReadStructuredKind,
        loc::Union{VTKPointData, VTKCellData}, name::String, ds::HDF5.Dataset, i::Int
    )
    nd = ndims(ds)
    A = ds[ntuple(_ -> Colon(), nd - 1)..., i]
    check_spatial_shape(kind, loc, size(A), HDF5.name(ds))
    return A
end

function check_spatial_shape(kind::ReadStructuredKind, loc, ashape::Tuple, name::String)
    dims = loc isa VTKPointData ? point_dims(kind) : cell_dims(kind)
    spatial_ncomp(ashape, dims) === nothing && error(
        "array $name has shape $ashape, which does not match the " *
            (loc isa VTKPointData ? "point" : "cell") * " dimensions $dims"
    )
    return nothing
end

# ---- geometry ----

function read_coordinates(r::VTKHDFReader{ReadRectilinear})
    check_static_geometry(r, "read_coordinates")
    return read_coordinate_slabs(r, (0, 0, 0))
end

function read_coordinates(s::VTKHDFTimeStep{<:VTKHDFReader{ReadRectilinear}})
    r = s.reader
    offsets = ntuple(d -> steps_entry(steps_info(r), COORD_NAMES[d] * "Offsets", s.index, 0), 3)
    return read_coordinate_slabs(r, offsets)
end

function read_coordinate_slabs(r::VTKHDFReader{ReadRectilinear}, offsets::NTuple{3, Int})
    return ntuple(3) do d
        ds = require_dataset(r.root, COORD_NAMES[d])
        n = r.kind.pdims[d]
        read_tuple_rows(ds, (offsets[d] + 1):(offsets[d] + n))
    end
end

function read_points(r::VTKHDFReader{ReadStructured})
    check_static_geometry(r, "read_points")
    ds = require_dataset(r.root, "Points")
    ndims(ds) == 4 || error("static StructuredGrid Points must be 4-dimensional")
    return read(ds)
end

function read_points(s::VTKHDFTimeStep{<:VTKHDFReader{ReadStructured}})
    r = s.reader
    ds = require_dataset(r.root, "Points")
    ndims(ds) == 5 || error("temporal StructuredGrid Points must be 5-dimensional")
    slab = steps_entry(steps_info(r), "PointOffsets", s.index, 0)
    0 <= slab < size(ds, 5) ||
        error("Steps/PointOffsets entry $slab is out of bounds for Points with $(size(ds, 5)) slabs")
    return ds[:, :, :, :, slab + 1]
end
