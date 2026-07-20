# HyperTreeGrid writer (spec 2.4, static).

mutable struct HTGState <: DatasetKind
    dims::NTuple{3, Int}
    branch_factor::Int
    total_cells::Int
    npieces::Int
    any_mask::Bool
end

"""
    vtkhdf_htg(filename; dimensions, branch_factor = 2, kwargs...)

Create a VTKHDF `HyperTreeGrid` file. The grid of trees has `dimensions`
(number of coordinate points per direction); the number of trees is the
product of `max(nᵢ - 1, 1)` over the directions (degenerate directions with a
single coordinate do not contribute). Pieces (partitions) are appended with
[`add_piece`](@ref).

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
        file, root, HTGState(Tuple(Int.(dimensions)), Int(branch_factor), 0, 0, false);
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

- `descriptors::AbstractVector{Bool}`: refinement bits for each tree in turn,
  level by level within the tree (each tree's deepest level carries no bits);
  bit-packed MSB-first on disk. Each refined cell must have exactly
  `branch_factor^dim` cells on the next level.
- `depth_per_tree`, `tree_ids`: depth and id of each tree in the piece.
- `number_of_cells_per_tree_depth`: cells per depth, tree by tree
  (`sum(depth_per_tree)` entries); its total is the piece's cell count.
- `x/y/zcoordinates`: tree coordinates for this piece (lengths matching
  `dimensions`).
- `mask`: optional `AbstractVector{Bool}` of length equal to the cell count,
  in descriptor cell order; refined cells cannot be masked.
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
    vtk.failed && error("a previous write to this file failed; the file is incomplete")
    kind = vtk.kind
    length(depth_per_tree) == length(tree_ids) ||
        throw(ArgumentError("depth_per_tree and tree_ids must have equal length"))
    all(>=(1), depth_per_tree) || throw(ArgumentError("tree depths must be at least 1"))
    length(number_of_cells_per_tree_depth) == sum(depth_per_tree; init = 0) ||
        throw(ArgumentError("number_of_cells_per_tree_depth must have sum(depth_per_tree) entries"))
    all(>=(1), number_of_cells_per_tree_depth) ||
        throw(ArgumentError("cell counts per tree depth must be at least 1"))
    ntrees_total = prod(max(d - 1, 1) for d in kind.dims)
    allunique(tree_ids) || throw(ArgumentError("tree_ids must be unique within a piece"))
    all(id -> 0 <= id < ntrees_total, tree_ids) ||
        throw(ArgumentError("tree_ids must be in 0:$(ntrees_total - 1)"))
    # per-tree structure: the root level has one cell, and only non-deepest
    # levels carry descriptor bits
    edim = count(>(1), kind.dims)
    fanout = kind.branch_factor^edim
    nbits = 0
    offset = 0
    for depth in depth_per_tree
        counts = view(number_of_cells_per_tree_depth, (offset + 1):(offset + depth))
        counts[1] == 1 || throw(ArgumentError("the root level of each tree must have exactly 1 cell"))
        nbits += sum(counts) - counts[end]  # deepest level carries no bits
        offset += depth
    end
    length(descriptors) == nbits || throw(
        ArgumentError(
            "descriptors must have exactly one bit per cell of every non-deepest tree level ($nbits), got $(length(descriptors))"
        )
    )
    ncells = sum(number_of_cells_per_tree_depth; init = 0)
    if mask !== nothing
        length(mask) == ncells ||
            throw(ArgumentError("mask must have one entry per cell ($ncells), got $(length(mask))"))
        kind.npieces > 0 && !kind.any_mask &&
            throw(ArgumentError("either all pieces or no piece must define a mask"))
    elseif kind.any_mask
        throw(ArgumentError("either all pieces or no piece must define a mask"))
    end
    # the descriptor bits must produce exactly the next level's cells, and a
    # refined cell cannot be masked (both required by the VTKHDF spec)
    bit = 0
    cell = 0
    offset = 0
    for depth in depth_per_tree
        counts = view(number_of_cells_per_tree_depth, (offset + 1):(offset + depth))
        for d in 1:(depth - 1)
            nrefined = 0
            for j in 1:counts[d]
                descriptors[bit + j] || continue
                nrefined += 1
                mask !== nothing && mask[cell + j] &&
                    throw(ArgumentError("refined cells cannot be masked (cell $(cell + j) of the piece)"))
            end
            counts[d + 1] == nrefined * fanout || throw(
                ArgumentError(
                    "invalid refinement: depth $(d + 1) of a tree has $(counts[d + 1]) cells, but its " *
                        "descriptors refine $nrefined cells at depth $d (branch factor " *
                        "$(kind.branch_factor), $edim dimensions => $(nrefined * fanout) cells)"
                )
            )
            bit += counts[d]
            cell += counts[d]
        end
        cell += counts[depth]  # deepest-level cells carry no descriptor bits
        offset += depth
    end
    for (coords, n) in zip((xcoordinates, ycoordinates, zcoordinates), kind.dims)
        length(coords) == n ||
            throw(ArgumentError("coordinate vector length $(length(coords)) does not match dimensions $(kind.dims)"))
    end
    for (name, data) in celldata  # preflight before mutating the file
        check_name(String(name))
        n = tuple_count(data)
        n == ncells || error("piece cell data $name has $n tuples, expected $ncells")
        check_tuple_data(vtk, "CellData", String(name), data)
    end
    root = vtk.root
    mutating(vtk) do
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
            append_tuple_data!(vtk, "CellData", String(name), data)
        end
        kind.total_cells += ncells
        kind.npieces += 1
        mask === nothing || (kind.any_mask = true)
    end
    return vtk
end

# an HTG closed without pieces still gets its (empty) datasets
function finalize_kind!(vtk, kind::HTGState)
    if kind.npieces == 0
        root = vtk.root
        for name in ("XCoordinates", "YCoordinates", "ZCoordinates")
            appendable(vtk, root, name, Float64, ())
        end
        appendable(vtk, root, "Descriptors", UInt8, ())
        for name in (
                "DescriptorsSize", "TreeIds", "DepthPerTree",
                "NumberOfCellsPerTreeDepth", "NumberOfTrees", "NumberOfDepths",
                "NumberOfCells",
            )
            appendable(vtk, root, name, Int64, ())
        end
    end
    return nothing
end

resolve_location(vtk, kind::HTGState, data) = VTKCellData()

write_array!(vtk::VTKHDFFile, kind::HTGState, loc::VTKCellData, name::AbstractString, data) =
    append_tuple_data!(vtk, "CellData", name, data)

write_array!(vtk::VTKHDFFile, kind::HTGState, loc::VTKPointData, name::AbstractString, data) =
    throw(ArgumentError("HyperTreeGrids cannot store point data"))

write_array!(vtk::VTKHDFFile, kind::HTGState, loc::VTKFieldData, name::AbstractString, data) =
    append_tuple_data!(vtk, "FieldData", name, data)

expected_totals(vtk, kind::HTGState) = Dict("CellData" => kind.total_cells)
