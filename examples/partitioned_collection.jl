# # PartitionedDataSetCollection: blocks and assembly
#
# Port of the `composite.hdf` example from the VTKHDF specification: a
# collection holding a PolyData block and an UnstructuredGrid block, plus an
# Assembly hierarchy grouping them. Blocks are written with the same API as
# standalone files; assembly references are HDF5 soft links, and a block may
# appear under several nodes. `vtkhdf_multiblock` writes a
# MultiBlockDataSet the same way.

using WriteVTKHDF

square = Float32[
    0 1 1 0
    0 0 1 1
    0 0 0 0
]
quad = [MeshCell(PolyData.Polys(), [1, 2, 3, 4])]

cube = Float64[
    0 1 1 0 0 1 1 0
    0 0 1 1 0 0 1 1
    0 0 0 0 1 1 1 1
]
hex = [MeshCell(VTKCellTypes.VTK_HEXAHEDRON, 1:8)]

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
