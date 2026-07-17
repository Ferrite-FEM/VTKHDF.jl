"""
    WriteVTKHDF

Write VTK data in the [VTKHDF file format](https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/index.html)
(HDF5-based, spec version 2.x). Writing only; supports static and temporal
data for all basic VTK dataset types as well as composite datasets.

Main entry points: [`vtkhdf_grid`](@ref), [`vtkhdf_table`](@ref),
[`vtkhdf_amr`](@ref), [`vtkhdf_htg`](@ref), [`vtkhdf_collection`](@ref),
[`vtkhdf_multiblock`](@ref) and [`write_timestep`](@ref).
"""
module WriteVTKHDF

using HDF5: HDF5, h5open
using VTKBase:
    VTKBase,
    MeshCell, VTKPolyhedron, PolyData,
    VTKCellTypes, VTKCellType,
    VTKPointData, VTKCellData, VTKFieldData, AbstractFieldData,
    AbstractVTKDataset,
    VTKImageData, VTKRectilinearGrid, VTKStructuredGrid,
    VTKUnstructuredGrid, VTKPolyData

export
    vtkhdf_grid, vtkhdf_table, vtkhdf_amr, vtkhdf_htg,
    vtkhdf_collection, vtkhdf_multiblock,
    write_timestep, add_partition,
    add_level, add_box,
    add_piece,
    add_node, add_block_ref, add_empty_block,
    # re-exports from VTKBase
    MeshCell, VTKPolyhedron, PolyData, VTKCellTypes,
    VTKPointData, VTKCellData, VTKFieldData,
    VTKImageData, VTKRectilinearGrid, VTKStructuredGrid,
    VTKUnstructuredGrid, VTKPolyData

include("h5helpers.jl")
include("data_arrays.jl")
include("cells.jl")
include("file.jl")
include("temporal.jl")
include("unstructured.jl")
include("polydata.jl")
include("image.jl")
include("rectilinear.jl")
include("structured.jl")
include("table.jl")
include("amr.jl")
include("hypertreegrid.jl")
include("composite.jl")
include("api.jl")

end # module WriteVTKHDF
