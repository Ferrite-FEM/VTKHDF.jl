# HyperTreeGrid writer (spec 2.4, static).

mutable struct HTGState <: DatasetKind
    dims::NTuple{3, Int}
    total_cells::Int
    npieces::Int
    any_mask::Bool
end

"""
    vtkhdf_htg(filename; dimensions, branch_factor = 2, kwargs...)

Create a VTKHDF `HyperTreeGrid` file. The grid of trees has `dimensions`
(number of coordinate points per direction; `(n₁-1)·(n₂-1)·(n₃-1)` trees).
Pieces (partitions) are appended with [`add_piece`](@ref).

Optional keywords: `transposed_root_indexing::Bool`,
`interface_normals_name`, `interface_intercepts_name`.

Temporal HyperTreeGrid writing is not supported.
"""
function vtkhdf_htg(filename::AbstractString; kwargs...)
    return init_htg(filename; kwargs...)
end

function init_htg(
        dest; dimensions, branch_factor::Integer = 2,
        transposed_root_indexing::Bool = false,
        interface_normals_name = nothing, interface_intercepts_name = nothing,
        compress = false, chunk_size = 0, temporal = false
    )
    temporal && throw(ArgumentError("temporal HyperTreeGrid writing is not supported"))
    length(dimensions) == 3 || throw(ArgumentError("dimensions must have 3 entries"))
    branch_factor in (2, 3) || throw(ArgumentError("branch_factor must be 2 or 3"))
    file, root = open_dest(dest)
    write_ascii_attribute(root, "Type", "HyperTreeGrid")
    write_version_attribute(root, (2, 4))
    HDF5.attrs(root)["Dimensions"] = Int64[dimensions...]
    HDF5.attrs(root)["BranchFactor"] = Int64(branch_factor)
    HDF5.attrs(root)["TransposedRootIndexing"] = Int64(transposed_root_indexing)
    interface_normals_name === nothing ||
        write_ascii_attribute(root, "InterfaceNormalsName", String(interface_normals_name))
    interface_intercepts_name === nothing ||
        write_ascii_attribute(root, "InterfaceInterceptsName", String(interface_intercepts_name))
    return make_vtkfile(
        file, root, HTGState(Tuple(Int.(dimensions)), 0, 0, false);
        temporal = false, compress, chunk_size, version = (2, 4)
    )
end

# MSB-first bit packing (vtkBitArray convention).
function pack_bits(bits::AbstractVector{Bool})
    bytes = zeros(UInt8, cld(length(bits), 8))
    for (i, b) in enumerate(bits)
        if b
            bytes[(i - 1) >> 3 + 1] |= 0x80 >> ((i - 1) & 7)
        end
    end
    return bytes
end

"""
    add_piece(
        htg; descriptors, depth_per_tree, tree_ids,
        number_of_cells_per_tree_depth, xcoordinates, ycoordinates,
        zcoordinates, mask = nothing, celldata = ()
    )

Append one piece (partition) to a HyperTreeGrid file:

- `descriptors::AbstractVector{Bool}`: refinement bits, level by level, for all
  trees of the piece (deepest level excluded); bit-packed MSB-first on disk.
- `depth_per_tree`, `tree_ids`: depth and id of each tree in the piece.
- `number_of_cells_per_tree_depth`: cells per depth, tree by tree
  (`sum(depth_per_tree)` entries); its total is the piece's cell count.
- `x/y/zcoordinates`: tree coordinates for this piece (lengths matching
  `dimensions`).
- `mask`: optional `AbstractVector{Bool}` of length equal to the cell count.
- `celldata`: iterable of `name => array` pairs (HyperTreeGrids have no point
  data).
"""
function add_piece(
        vtk::VTKHDFFile{HTGState};
        descriptors::AbstractVector{Bool},
        depth_per_tree::AbstractVector{<:Integer},
        tree_ids::AbstractVector{<:Integer},
        number_of_cells_per_tree_depth::AbstractVector{<:Integer},
        xcoordinates::AbstractVector{<:Real},
        ycoordinates::AbstractVector{<:Real},
        zcoordinates::AbstractVector{<:Real},
        mask::Union{Nothing, AbstractVector{Bool}} = nothing,
        celldata = ()
    )
    vtk.isopen || error("file is closed")
    kind = vtk.kind
    length(depth_per_tree) == length(tree_ids) ||
        throw(ArgumentError("depth_per_tree and tree_ids must have equal length"))
    length(number_of_cells_per_tree_depth) == sum(depth_per_tree; init = 0) ||
        throw(ArgumentError("number_of_cells_per_tree_depth must have sum(depth_per_tree) entries"))
    ncells = sum(number_of_cells_per_tree_depth; init = 0)
    if mask !== nothing
        length(mask) == ncells ||
            throw(ArgumentError("mask must have one entry per cell ($ncells), got $(length(mask))"))
        kind.npieces > 0 && !kind.any_mask &&
            throw(ArgumentError("either all pieces or no piece must define a mask"))
        kind.any_mask = true
    elseif kind.any_mask
        throw(ArgumentError("either all pieces or no piece must define a mask"))
    end
    for (coords, n) in zip((xcoordinates, ycoordinates, zcoordinates), kind.dims)
        length(coords) == n ||
            throw(ArgumentError("coordinate vector length $(length(coords)) does not match dimensions $(kind.dims)"))
    end
    root = vtk.root
    append_rows(appendable(vtk, root, "XCoordinates", Float64, ()), Vector{Float64}(xcoordinates))
    append_rows(appendable(vtk, root, "YCoordinates", Float64, ()), Vector{Float64}(ycoordinates))
    append_rows(appendable(vtk, root, "ZCoordinates", Float64, ()), Vector{Float64}(zcoordinates))
    # Descriptors and Mask are bit-packed; every piece starts on a fresh byte.
    append_rows(appendable(vtk, root, "Descriptors", UInt8, ()), pack_bits(descriptors))
    append_rows(appendable(vtk, root, "DescriptorsSize", Int64, ()), Int64(length(descriptors)))
    append_rows(appendable(vtk, root, "TreeIds", Int64, ()), Vector{Int64}(tree_ids))
    append_rows(appendable(vtk, root, "DepthPerTree", Int64, ()), Vector{Int64}(depth_per_tree))
    append_rows(
        appendable(vtk, root, "NumberOfCellsPerTreeDepth", Int64, ()),
        Vector{Int64}(number_of_cells_per_tree_depth)
    )
    append_rows(appendable(vtk, root, "NumberOfTrees", Int64, ()), Int64(length(tree_ids)))
    append_rows(appendable(vtk, root, "NumberOfDepths", Int64, ()), Int64(length(number_of_cells_per_tree_depth)))
    append_rows(appendable(vtk, root, "NumberOfCells", Int64, ()), Int64(ncells))
    mask === nothing ||
        append_rows(appendable(vtk, root, "Mask", UInt8, ()), pack_bits(mask))
    for (name, data) in celldata
        n = tuple_count(data)
        n == ncells || error("piece cell data $name has $n tuples, expected $ncells")
        append_tuple_data!(vtk, "CellData", String(name), data)
    end
    kind.total_cells += ncells
    kind.npieces += 1
    return vtk
end

resolve_location(vtk, kind::HTGState, data) = VTKCellData()

write_array!(vtk::VTKHDFFile, kind::HTGState, loc::VTKCellData, name::AbstractString, data) =
    append_tuple_data!(vtk, "CellData", name, data)

write_array!(vtk::VTKHDFFile, kind::HTGState, loc::VTKPointData, name::AbstractString, data) =
    throw(ArgumentError("HyperTreeGrids cannot store point data"))

write_array!(vtk::VTKHDFFile, kind::HTGState, loc::VTKFieldData, name::AbstractString, data) =
    append_tuple_data!(vtk, "FieldData", name, data)

expected_totals(vtk, kind::HTGState) = Dict("CellData" => kind.total_cells)
