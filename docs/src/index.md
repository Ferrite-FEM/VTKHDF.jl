# VTKHDF.jl

Write and read VTK data in the
[VTKHDF file format](https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/index.html),
the HDF5-based successor to the VTK XML formats. The files can be opened in
ParaView, VisIt and VTK. The writing API follows
[WriteVTK.jl](https://github.com/JuliaVTK/WriteVTK.jl) and shares its cell/data
types through [VTKBase.jl](https://github.com/JuliaVTK/VTKBase.jl).
[`vtkhdf_open`](@ref) reads those files back, and also spec-conforming files
from other writers.

## Why VTKHDF?

- Time series live in one file. The geometry of a fixed mesh is stored once
  and each step appends only its data — no more `.pvd` collections with one
  `.vtu` file per step.
- It is built on HDF5: standard tooling (`h5dump`, `h5py`, HDF5.jl) can
  inspect the files, and chunked storage and gzip compression are available.
- One format covers all the dataset types, from image data to composite
  multi-block datasets.

## Installation

```julia
import Pkg
Pkg.add(url = "https://github.com/Ferrite-FEM/VTKHDF.jl")
```

## Quick start

A grid is points plus cells connecting them. [`vtkhdf_grid`](@ref) writes
them to a file; data arrays are attached by indexing:

```@example quickstart
using VTKHDF

points = [
    0.0 1.0 0.0 0.0
    0.0 0.0 1.0 0.0
    0.0 0.0 0.0 1.0
]
cells = [MeshCell(VTKCellTypes.VTK_TETRA, [1, 2, 3, 4])]

vtkhdf_grid("output", points, cells) do vtk
    vtk["temperature"] = [0.0, 1.0, 1.0, 1.0]   # one value per point
    vtk["material", VTKCellData()] = [1]        # one value per cell
end
nothing #hide
```

This wrote `output.vtkhdf`, which can be opened in ParaView or VisIt.
[`vtkhdf_open`](@ref) reads it back with the same indexing:

```@example quickstart
vtkhdf_open("output") do r
    r["temperature"]
end
```

```@example quickstart
rm("output.vtkhdf") #hide
nothing #hide
```

The [Manual](@ref manual) covers writing and [reading](@ref reading) all
dataset types, including time series. The Examples section contains
complete programs, one per reference file of the VTKHDF specification.

## Supported dataset types

| VTKHDF type | write | read | temporal |
|---|---|---|---|
| UnstructuredGrid (incl. polyhedra, multiple partitions) | ✓ | ✓ | ✓ |
| PolyData (vertices/lines/polygons/strips) | ✓ | ✓ | ✓ |
| ImageData | ✓ | ✓ | ✓ |
| RectilinearGrid | ✓ | ✓ | ✓ |
| StructuredGrid | ✓ | ✓ | ✓ |
| Table | ✓ | ✓ | ✓ |
| OverlappingAMR | ✓ | ✓ | – |
| HyperTreeGrid | ✓ | ✓ | – |
| PartitionedDataSetCollection / MultiBlockDataSet | ✓ | ✓ | ✓ |

Temporal OverlappingAMR/HyperTreeGrid and MPI-parallel collective writing are
currently out of scope (VTK's own writer does not produce them either).

!!! note "Opening files in ParaView/VTK"
    RectilinearGrid, StructuredGrid and Table are recent additions to the
    VTKHDF specification (2.7/2.8) and require a VTK build new enough to
    read them. All other types round-trip against VTK 9.6's `vtkHDFReader`
    in this package's test suite. This package's own [reading
    support](@ref reading) handles all types listed above.

## Development

This package was mainly developed with AI assistance (Claude Code), guided by
the VTKHDF specification, with the output validated against the VTK reference
implementation.
