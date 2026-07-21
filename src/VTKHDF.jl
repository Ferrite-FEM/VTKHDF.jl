"""
    VTKHDF

Write and read VTK data in the [VTKHDF file format](https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/index.html)
(HDF5-based, spec version 2.x). Supports static and temporal data for all
basic VTK dataset types and composite datasets.

Main entry points for writing: [`vtkhdf_grid`](@ref), [`vtkhdf_table`](@ref),
[`vtkhdf_amr`](@ref), [`vtkhdf_htg`](@ref), [`vtkhdf_collection`](@ref),
[`vtkhdf_multiblock`](@ref) and [`write_timestep`](@ref); for reading:
[`vtkhdf_open`](@ref) and [`read_timestep`](@ref).
"""
module VTKHDF

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
    # reading
    vtkhdf_open, read_timestep,
    read_points, read_cells, read_coordinates,
    VTKRowData,
    # re-exports from VTKBase
    MeshCell, VTKPolyhedron, PolyData, VTKCellTypes,
    VTKPointData, VTKCellData, VTKFieldData,
    VTKImageData, VTKRectilinearGrid, VTKStructuredGrid,
    VTKUnstructuredGrid, VTKPolyData

include("h5helpers.jl")
include("write/data_arrays.jl")
include("write/cells.jl")
include("write/core.jl")
include("write/temporal.jl")
include("write/unstructured.jl")
include("write/polydata.jl")
include("write/image.jl")
include("write/rectilinear.jl")
include("write/structured.jl")
include("write/table.jl")
include("write/amr.jl")
include("write/hypertreegrid.jl")
include("write/composite.jl")
include("write/api.jl")
include("read/core.jl")
include("read/unstructured.jl")
include("read/polydata.jl")
include("read/structured.jl")
include("read/table.jl")
include("read/amr.jl")
include("read/htg.jl")
include("read/composite.jl")

end # module VTKHDF
