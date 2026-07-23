# # PartitionedDataSetCollection: blocks and assembly
#
# Port of the `composite.hdf` example from the VTKHDF specification: a
# collection holding a PolyData block and an UnstructuredGrid block, plus an
# Assembly hierarchy grouping them. Blocks are written with the same API as
# standalone files; assembly references are HDF5 soft links, and a block may
# appear under several nodes. `vtkhdf_multiblock` writes a
# MultiBlockDataSet the same way. The two blocks in ParaView — `Solid`
# colored by `Temperature`, with the `Surface` square moved aside (it
# coincides with a face of the cube). The translation is a ParaView display
# transform used only for the screenshot; the file keeps the square on the
# cube face.
#
# **What you'll learn:** how to create named blocks, organize block references
# in an assembly tree, and read both structures back.
#
# ![Collection blocks: cube and square](../assets/examples/partitioned_collection-light.png)
# ![Collection blocks: cube and square](../assets/examples/partitioned_collection-dark.png)

using VTKHDF

square = Float32[
    0 1 1 0
    0 0 1 1
    0 0 0 0
]
quad = [MeshCell(PolyData.Polys(), [1, 2, 3, 4])];

cube = Float64[
    0 1 1 0 0 1 1 0
    0 0 1 1 0 0 1 1
    0 0 0 0 1 1 1 1
]
hex = [MeshCell(VTKCellTypes.VTK_HEXAHEDRON, 1:8)];

vtkhdf_collection("composite") do col
    surface = vtkhdf_grid(col, "Surface", square, quad)
    surface["Materials", VTKCellData()] = [1]
    solid = vtkhdf_grid(col, "Solid", cube, hex)
    solid["Temperature"] = Float64.(1:8)

    add_block_ref(col, surface)
    solids = add_node(col, "solids")
    add_block_ref(solids, solid)
    both = add_node(col, "everything")
    add_block_ref(both, surface)
    add_block_ref(both, solid)
end
nothing #hide

# ## Reading it back
#
# Blocks of a composite file are readers themselves:

col = vtkhdf_open("composite")
keys(col)

#-

col["Solid"]["Temperature"]

# The assembly hierarchy comes back as a tree of node names and block
# references:

asm = VTKHDF.read_assembly(col)
[(n.name, n.blocks) for n in asm.children]

#-

close(col)
