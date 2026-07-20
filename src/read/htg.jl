# HyperTreeGrid reader (static; temporal HTG is rejected at open).
#
# All per-piece quantities are concatenations over the pieces; the
# bit-packed Descriptors/Mask blocks start on a fresh byte per piece.

struct ReadHTG <: ReaderKind
    dims::NTuple{3, Int}
    branch_factor::Int
    transposed_root_indexing::Bool
    trees_prefix::Vector{Int}     # NumberOfTrees
    depths_prefix::Vector{Int}    # NumberOfDepths
    descbits_prefix::Vector{Int}  # DescriptorsSize
    cells_prefix::Vector{Int}     # NumberOfCells
end

reader_type_string(::ReadHTG) = "HyperTreeGrid"
supports_temporal(::ReadHTG) = false
data_locations(::ReadHTG) = (VTKCellData(), VTKFieldData())

function read_htg_kind(root::HDF5.Group)
    dims = read_attr_tuple(root, "Dimensions", Int, Val(3))
    dims === nothing && error("$(HDF5.name(root)) has no Dimensions attribute")
    bf = haskey(HDF5.attrs(root), "BranchFactor") ?
        Int(HDF5.read_attribute(root, "BranchFactor")) : 2
    transposed = haskey(HDF5.attrs(root), "TransposedRootIndexing") &&
        HDF5.read_attribute(root, "TransposedRootIndexing") != 0
    tp = count_prefix(root, "NumberOfTrees")
    dp = count_prefix(root, "NumberOfDepths")
    bp = count_prefix(root, "DescriptorsSize")
    cp = count_prefix(root, "NumberOfCells")
    npieces = npartitions_prefix(tp)
    (npartitions_prefix(dp) == npieces && npartitions_prefix(bp) == npieces && npartitions_prefix(cp) == npieces) ||
        error("inconsistent piece counts in the HyperTreeGrid datasets")
    return ReadHTG(dims, bf, transposed, tp, dp, bp, cp)
end

grid_info(r::VTKHDFReader{ReadHTG}) = (
    dimensions = r.kind.dims, branch_factor = r.kind.branch_factor,
    transposed_root_indexing = r.kind.transposed_root_indexing,
)

ncells(r::VTKHDFReader{ReadHTG}) = (check_open(r); r.kind.cells_prefix[end])

"""
    npieces(r) -> Int

Number of pieces (partitions) of a HyperTreeGrid file.
"""
npieces(r::VTKHDFReader{ReadHTG}) = (check_open(r); npartitions_prefix(r.kind.trees_prefix))

# MSB-first bit unpacking (inverse of pack_bits).
function unpack_bits(bytes::AbstractVector{UInt8}, nbits::Int)
    length(bytes) >= cld(nbits, 8) ||
        error("bit-packed dataset too short: $(length(bytes)) bytes for $nbits bits")
    bits = Vector{Bool}(undef, nbits)
    for i in 1:nbits
        bits[i] = (bytes[(i - 1) >> 3 + 1] >> (7 - ((i - 1) & 7))) & 0x01 == 0x01
    end
    return bits
end

# Byte range of piece p in a per-piece bit-packed dataset (each piece starts
# on a fresh byte).
function piece_byte_range(bits_prefix::Vector{Int}, p::Int)
    start = sum(i -> cld(bits_prefix[i + 1] - bits_prefix[i], 8), 1:(p - 1); init = 0)
    nbytes = cld(bits_prefix[p + 1] - bits_prefix[p], 8)
    return (start + 1):(start + nbytes)
end

"""
    htg_piece(r, p) -> NamedTuple

Piece `p` (1-based) of a HyperTreeGrid file, as the exact inverse of
[`add_piece`](@ref): a NamedTuple with `descriptors::Vector{Bool}`,
`tree_ids`, `depth_per_tree`, `number_of_cells_per_tree_depth`,
`xcoordinates`/`ycoordinates`/`zcoordinates`,
`mask::Union{Nothing, Vector{Bool}}` and
`celldata::Dict{String, Any}` (this piece's slice of every cell-data array).
"""
function htg_piece(r::VTKHDFReader{ReadHTG}, p::Integer)
    check_open(r)
    k = r.kind
    np = npieces(r)
    1 <= p <= np || error("piece $p does not exist; the file has $np pieces")
    p = Int(p)
    tree_ids = Int.(read_tuple_rows(require_dataset(r.root, "TreeIds"), (k.trees_prefix[p] + 1):k.trees_prefix[p + 1]))
    depth_per_tree = Int.(read_tuple_rows(require_dataset(r.root, "DepthPerTree"), (k.trees_prefix[p] + 1):k.trees_prefix[p + 1]))
    ncptd = Int.(
        read_tuple_rows(
            require_dataset(r.root, "NumberOfCellsPerTreeDepth"),
            (k.depths_prefix[p] + 1):k.depths_prefix[p + 1]
        )
    )
    nbits = k.descbits_prefix[p + 1] - k.descbits_prefix[p]
    desc_bytes = read_tuple_rows(require_dataset(r.root, "Descriptors"), piece_byte_range(k.descbits_prefix, p))
    descriptors = unpack_bits(desc_bytes, nbits)
    ncells_p = k.cells_prefix[p + 1] - k.cells_prefix[p]
    mask = nothing
    if haskey(r.root, "Mask")
        mask_bytes = read_tuple_rows(get_dataset(r.root, "Mask"), piece_byte_range(k.cells_prefix, p))
        mask = unpack_bits(mask_bytes, ncells_p)
    end
    coords = ntuple(3) do d
        ds = require_dataset(r.root, COORD_NAMES[d])
        n = k.dims[d]
        read_tuple_rows(ds, ((p - 1) * n + 1):(p * n))
    end
    celldata = Dict{String, Any}()
    for name in keys(r, VTKCellData())
        ds = get_dataset(r.root["CellData"]::HDF5.Group, name)
        celldata[name] = read_tuple_rows(ds, (k.cells_prefix[p] + 1):k.cells_prefix[p + 1])
    end
    return (
        descriptors = descriptors, tree_ids = tree_ids, depth_per_tree = depth_per_tree,
        number_of_cells_per_tree_depth = ncptd,
        xcoordinates = coords[1], ycoordinates = coords[2], zcoordinates = coords[3],
        mask = mask, celldata = celldata,
    )
end
