# VTKHDF.jl

Write and read VTK data in the
[VTKHDF file format](https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/index.html) —
the HDF5-based successor to the VTK XML formats — readable by ParaView, VisIt
and VTK. The writing API aims to feel familiar to
[WriteVTK.jl](https://github.com/JuliaVTK/WriteVTK.jl) users and shares its
cell/data types through [VTKBase.jl](https://github.com/JuliaVTK/VTKBase.jl);
[`vtkhdf_open`](@ref) reads the same files (and spec-conforming files from
other writers) back.

## Why VTKHDF?

- **One file per simulation, not per step.** VTKHDF has native time series
  support: a temporal file holds all steps, and for a fixed mesh the geometry
  is stored *once* while each step appends only its data. This replaces the
  `.pvd` + one-`.vtu`-per-step workflow.
- **HDF5 underneath**: standard tooling (`h5dump`, `h5py`, HDF5.jl) can
  inspect files; chunked storage and gzip compression come for free.
- **A single format for everything** from image data to composite
  multi-block datasets.

## Installation

```julia
pkg> add https://github.com/Ferrite-FEM/VTKHDF.jl
```

## Quick start

```julia
using VTKHDF

points = rand(3, 100)
cells = [MeshCell(VTKCellTypes.VTK_TETRA, [1, 2, 3, 4])]

vtkhdf_grid("output", points, cells) do vtk
    vtk["temperature"] = rand(100)                  # point data (auto-detected)
    vtk["material", VTKCellData()] = [1]
end

vtkhdf_open("output") do r                          # read it back
    r["temperature"], read_points(r), read_cells(r)
end
```

See the [Manual](@ref manual) for all dataset types, temporal writing and the
reading API, and the Examples section for complete programs — one per
reference file of the VTKHDF specification.

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
    support](@ref manual) handles all types listed above.

## Development

This package was mainly developed with AI assistance (Claude Code), guided by
the VTKHDF specification, with the output validated against the VTK reference
implementation.
