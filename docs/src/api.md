# API reference

## Common

Cell types (`MeshCell`, `VTKPolyhedron`, `VTKCellTypes`, `PolyData.*`) and
data locations (`VTKPointData`, `VTKCellData`, `VTKFieldData`) are
re-exported from [VTKBase.jl](https://github.com/JuliaVTK/VTKBase.jl) and are
used for both writing and reading.

```@docs
VTKRowData
```

## Writing

### File constructors

```@docs
vtkhdf_grid
vtkhdf_table
vtkhdf_amr
vtkhdf_htg
vtkhdf_collection
vtkhdf_multiblock
```

### Data arrays

```@docs
Base.setindex!(::VTKHDF.VTKHDFFile, ::Any, ::AbstractString)
```

### Time series

```@docs
write_timestep
```

### Partitions, levels and pieces

```@docs
add_partition
add_level
add_box
add_piece
```

### Composite assembly

```@docs
add_node
add_block_ref
add_empty_block
```

## Reading

### Opening files and reading geometry

```@docs
vtkhdf_open
read_timestep
read_points
read_cells
read_coordinates
```

### Data arrays

```@docs
Base.getindex(::VTKHDF.VTKHDFReader, ::AbstractString)
Base.keys(::VTKHDF.VTKHDFReader, ::VTKHDF.VTKBase.AbstractFieldData)
Base.keys(::VTKHDF.VTKHDFCollectionReader)
Base.getindex(::VTKHDF.VTKHDFCollectionReader, ::AbstractString)
```

### Metadata and structure

Unexported (call as `VTKHDF.f`):

```@docs
VTKHDF.dataset_type
VTKHDF.file_version
VTKHDF.is_temporal
VTKHDF.nsteps
VTKHDF.time_values
VTKHDF.time_value
VTKHDF.grid_info
VTKHDF.npoints
VTKHDF.ncells
VTKHDF.nrows
VTKHDF.npartitions
VTKHDF.partition_ranges
VTKHDF.active_attributes
VTKHDF.data_attributes
VTKHDF.read_assembly
VTKHDF.nlevels
VTKHDF.amr_level
VTKHDF.level_info
VTKHDF.npieces
VTKHDF.htg_piece
```

## Module

```@docs
VTKHDF.VTKHDF
```
