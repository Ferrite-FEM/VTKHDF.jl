# Composite datasets: PartitionedDataSetCollection and MultiBlockDataSet.
#
# Blocks are simple datasets stored as subgroups of /VTKHDF (each with its own
# Type/Version and, for PDC, an Index attribute). The block hierarchy is
# described in the Assembly group with soft links; all involved groups track
# link creation order as required by the spec.

mutable struct VTKHDFCollection
    file::HDF5.File
    root::HDF5.Group
    multiblock::Bool
    blocks::Vector{VTKHDFFile}
    next_index::Int
    isopen::Bool
end

struct AssemblyNode
    col::VTKHDFCollection
    group::HDF5.Group
end

"""
    vtkhdf_collection(filename) -> collection
    vtkhdf_multiblock(filename) -> collection

Create a composite VTKHDF file of type `PartitionedDataSetCollection`
(respectively `MultiBlockDataSet`). Blocks are created with the usual
constructors taking `(collection, name)` instead of a filename, e.g.
`vtkhdf_grid(col, "mesh", points, cells)`, and behave like standalone files
(including temporal writing; all temporal blocks must use the same time
values). The block hierarchy is defined with [`add_node`](@ref) and
[`add_block_ref`](@ref); [`add_empty_block`](@ref) creates a type-less block.

Supports the do-block form. Closing the collection closes all blocks.
"""
function vtkhdf_collection(filename::AbstractString; multiblock::Bool = false)
    file = h5open(add_extension(filename), "w")
    root = HDF5.create_group(file, "VTKHDF"; track_order = true)
    write_ascii_attribute(
        root, "Type",
        multiblock ? "MultiBlockDataSet" : "PartitionedDataSetCollection"
    )
    write_version_attribute(root, (2, 1))
    # always present, even when empty: vtkHDFReader crashes without it
    HDF5.create_group(root, "Assembly"; track_order = true)
    return VTKHDFCollection(file, root, multiblock, VTKHDFFile[], 0, true)
end

"""
    vtkhdf_multiblock(filename)

Shorthand for `vtkhdf_collection(filename; multiblock = true)`; see
[`vtkhdf_collection`](@ref).
"""
vtkhdf_multiblock(filename::AbstractString) = vtkhdf_collection(filename; multiblock = true)

assembly_group(col::VTKHDFCollection) =
    get_or_create_group(col.root, "Assembly"; track_order = true)

# Destination handle used by the init_* constructors: either a filename or a
# block inside a collection.
struct BlockDest
    col::VTKHDFCollection
    name::String
end

open_dest(filename::AbstractString) = open_vtkhdf(filename)

function open_dest(d::BlockDest)
    d.col.isopen || error("collection is closed")
    check_name(d.name, "block name")
    haskey(d.col.root, d.name) && throw(ArgumentError("block $(d.name) already exists"))
    grp = HDF5.create_group(d.col.root, d.name; track_order = true)
    if !d.col.multiblock
        HDF5.attrs(grp)["Index"] = Int64(d.col.next_index)
        d.col.next_index += 1
    end
    return nothing, grp
end

register_block(dest, vtk) = vtk
function register_block(dest::BlockDest, vtk::VTKHDFFile)
    push!(dest.col.blocks, vtk)
    return vtk
end

"""
    add_empty_block(col, name)

Add a block with no dataset type to a composite file (a placeholder allowed by
the spec).
"""
function add_empty_block(col::VTKHDFCollection, name::AbstractString)
    col.isopen || error("collection is closed")
    check_name(name, "block name")
    haskey(col.root, name) && throw(ArgumentError("block $name already exists"))
    HDF5.create_group(col.root, name; track_order = true)
    return col
end

"""
    add_node(col_or_node, name) -> node

Create a group node in the Assembly hierarchy of a composite file (under the
root assembly, or nested under another node).
"""
function add_node(col::VTKHDFCollection, name::AbstractString)
    check_name(name, "assembly node name")
    return AssemblyNode(col, get_or_create_group(assembly_group(col), name; track_order = true))
end
function add_node(node::AssemblyNode, name::AbstractString)
    check_name(name, "assembly node name")
    return AssemblyNode(node.col, get_or_create_group(node.group, name; track_order = true))
end

"""
    add_block_ref(col_or_node, block)

Reference a block from the Assembly hierarchy (a soft link named after the
block). A block may be referenced from several nodes.
"""
function add_block_ref(node::AssemblyNode, blk::VTKHDFFile)
    target = HDF5.name(blk.root)
    create_soft_link(node.group, basename(target), target)
    return node
end
function add_block_ref(col::VTKHDFCollection, blk::VTKHDFFile)
    target = HDF5.name(blk.root)
    create_soft_link(assembly_group(col), basename(target), target)
    return col
end

function Base.close(col::VTKHDFCollection)
    col.isopen || return nothing
    # temporal blocks must agree on their time values
    ref = nothing
    for blk in col.blocks
        blk.temporal || continue
        if ref === nothing
            ref = blk.step_values
        elseif blk.step_values != ref
            error("temporal blocks have different time steps: $(blk.step_values) vs $ref")
        end
    end
    for blk in col.blocks
        close(blk)
    end
    col.isopen = false
    close(col.root)
    close(col.file)
    return nothing
end

function Base.show(io::IO, col::VTKHDFCollection)
    return print(
        io, "VTKHDFCollection (", col.multiblock ? "MultiBlockDataSet" : "PartitionedDataSetCollection",
        ", ", length(col.blocks), " blocks, ", col.isopen ? "open" : "closed", ")"
    )
end
