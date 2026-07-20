# Composite dataset reader: PartitionedDataSetCollection / MultiBlockDataSet.
#
# Blocks share the collection's file handle: block readers hold a reference
# to the collection and are invalidated when it closes; closing a block
# reader never closes the file.

mutable struct VTKHDFCollectionReader
    file::HDF5.File
    root::HDF5.Group
    multiblock::Bool
    version::NTuple{2, Int}
    block_names::Vector{String}   # PDC: ordered by Index; MBD: creation order
    blocks::Dict{String, Any}     # opened block readers, by name
    isopen::Bool
end

function open_collection(file::HDF5.File, root::HDF5.Group, type::String)
    multiblock = type == "MultiBlockDataSet"
    version = read_version_attribute(root)
    names = [n for n in links_in_creation_order(root) if n != "Assembly"]
    if !multiblock
        # PDC blocks are ordered by their globally unique Index attribute
        indices = map(names) do n
            grp = root[n]
            grp isa HDF5.Group && haskey(HDF5.attrs(grp), "Index") ||
                error("PartitionedDataSetCollection block $n has no Index attribute")
            Int(HDF5.read_attribute(grp, "Index"))
        end
        allunique(indices) ||
            error("PartitionedDataSetCollection block Index attributes are not unique")
        names = names[sortperm(indices)]
    end
    return VTKHDFCollectionReader(
        file, root, multiblock, version, names, Dict{String, Any}(), true
    )
end

# Placeholder handle for a type-less (empty) block.
struct VTKHDFEmptyBlockReader
    parent::VTKHDFCollectionReader
    name::String
end

dataset_type(b::VTKHDFEmptyBlockReader) = (check_open(b.parent); nothing)
Base.show(io::IO, b::VTKHDFEmptyBlockReader) =
    print(io, "VTKHDFEmptyBlockReader(", repr(b.name), ")")

dataset_type(col::VTKHDFCollectionReader) =
    (check_open(col); col.multiblock ? "MultiBlockDataSet" : "PartitionedDataSetCollection")
file_version(col::VTKHDFCollectionReader) = (check_open(col); col.version)

function check_open(col::VTKHDFCollectionReader)
    col.isopen || error("collection reader is closed")
    return nothing
end

function Base.close(col::VTKHDFCollectionReader)
    col.isopen || return nothing
    col.isopen = false
    for (_, blk) in col.blocks
        blk isa VTKHDFReader && (blk.isopen = false)
    end
    close(col.root)
    close(col.file)
    return nothing
end

Base.isopen(col::VTKHDFCollectionReader) = col.isopen

function Base.show(io::IO, col::VTKHDFCollectionReader)
    return print(
        io, "VTKHDFCollectionReader (", dataset_type(col), ", ",
        length(col.block_names), " blocks, ", col.isopen ? "open" : "closed", ")"
    )
end

"""
    keys(col) -> Vector{String}

The block names of a composite file, ordered by the `Index` attribute for
`PartitionedDataSetCollection` and by creation order for
`MultiBlockDataSet`. Blocks are opened by indexing: `col["mesh"]`.
"""
Base.keys(col::VTKHDFCollectionReader) = (check_open(col); copy(col.block_names))
Base.haskey(col::VTKHDFCollectionReader, name::AbstractString) =
    (check_open(col); String(name) in col.block_names)

"""
    col[name] -> block reader

Open block `name` of a composite file. Typed blocks return a reader with
the full reading API; type-less placeholder blocks return a handle whose
`dataset_type` is `nothing`. Block readers share the collection's file
handle and are invalidated when the collection is closed.
"""
function Base.getindex(col::VTKHDFCollectionReader, name::AbstractString)
    check_open(col)
    name = String(name)
    haskey(col, name) ||
        error("no block named $(repr(name)); available blocks: $(join(col.block_names, ", "))")
    return get!(col.blocks, name) do
        grp = col.root[name]::HDF5.Group
        if haskey(HDF5.attrs(grp), "Type")
            open_reader(nothing, grp, col)
        else
            VTKHDFEmptyBlockReader(col, name)
        end
    end
end

# ---- temporal metadata (per spec, blocks carry their own Steps groups and
# must agree on the time values) ----

typed_blocks(col::VTKHDFCollectionReader) =
    (b for b in (col[name] for name in col.block_names) if b isa VTKHDFReader)

is_temporal(col::VTKHDFCollectionReader) =
    (check_open(col); any(is_temporal, typed_blocks(col)))

function collection_time_values(col::VTKHDFCollectionReader)
    ref = nothing
    refname = ""
    for name in col.block_names
        blk = col[name]
        blk isa VTKHDFReader && is_temporal(blk) || continue
        vals = blk.steps.values
        if ref === nothing
            ref = vals
            refname = name
        elseif vals != ref
            error("temporal blocks $refname and $name have different time values")
        end
    end
    ref === nothing && error("not a temporal file: no block has time steps")
    return ref
end

nsteps(col::VTKHDFCollectionReader) = (check_open(col); length(collection_time_values(col)))
time_values(col::VTKHDFCollectionReader) = (check_open(col); copy(collection_time_values(col)))

# ---- assembly ----

"""
    read_assembly(col) -> AssemblyTreeNode

The Assembly hierarchy of a composite file as a tree. Each node has `name`,
`children::Vector{AssemblyTreeNode}` (nested assembly nodes) and
`blocks::Vector{String}` (names of the blocks referenced from this node),
both in creation order. The root node is named `"Assembly"`.
"""
function read_assembly(col::VTKHDFCollectionReader)
    check_open(col)
    haskey(col.root, "Assembly") ||
        return AssemblyTreeNode("Assembly", AssemblyTreeNode[], String[])
    block_addrs = Dict{Tuple{UInt, UInt64}, String}()
    for name in col.block_names
        obj = col.root[name]
        obj isa HDF5.Group && (block_addrs[object_address(obj)] = name)
    end
    visited = Set{Tuple{UInt, UInt64}}()
    return read_assembly_node(col, col.root["Assembly"]::HDF5.Group, "Assembly", block_addrs, visited)
end

# File-unique address of an HDF5 object, used to recognize which block a
# soft link resolves to (HDF5.name reports the access path, not the target).
function object_address(obj)
    info = HDF5.API.h5o_get_info1(obj)
    return (UInt(info.fileno), UInt64(info.addr))
end

struct AssemblyTreeNode
    name::String
    children::Vector{AssemblyTreeNode}
    blocks::Vector{String}
end

function Base.show(io::IO, node::AssemblyTreeNode)
    return print(
        io, "AssemblyTreeNode(", repr(node.name), ", ", length(node.children),
        " children, ", length(node.blocks), " blocks)"
    )
end

function read_assembly_node(
        col::VTKHDFCollectionReader, group::HDF5.Group, name::String,
        block_addrs::Dict, visited::Set
    )
    addr = object_address(group)
    addr in visited &&
        error("assembly hierarchy contains a cycle or repeated node at $(HDF5.name(group))")
    push!(visited, addr)
    children = AssemblyTreeNode[]
    blocks = String[]
    path = HDF5.name(group)
    for link in links_in_creation_order(group)
        obj = try
            group[link]
        catch
            error("assembly link $path/$link cannot be resolved (dangling soft link)")
        end
        obj isa HDF5.Group ||
            error("assembly entry $path/$link is not a group")
        # Classification is by the resolved object's identity (HDF5.jl has no
        # working link-introspection API against libhdf5 1.14): a link
        # resolving to a top-level block group is a block reference (the
        # writer stores those as soft links); any other group is treated as a
        # nested assembly node. A soft link pointing into the Assembly itself
        # makes some node resolve twice and is rejected by the `visited` set.
        blockname = get(block_addrs, object_address(obj), nothing)
        if blockname !== nothing
            push!(blocks, blockname)
        else
            push!(children, read_assembly_node(col, obj, link, block_addrs, visited))
        end
    end
    return AssemblyTreeNode(name, children, blocks)
end

# ---- creation-order link iteration ----

# Links of a group in creation order, falling back to name order for groups
# without a creation-order index (files from other writers).
function links_in_creation_order(group::HDF5.Group)
    names = String[]
    try
        n = length(keys(group))
        for i in 0:(n - 1)
            name = HDF5.API.h5l_get_name_by_idx(
                group, ".", HDF5.API.H5_INDEX_CRT_ORDER, HDF5.API.H5_ITER_INC,
                i, HDF5.API.H5P_DEFAULT
            )
            push!(names, name)
        end
    catch
        return collect(String, keys(group))
    end
    return names
end
