# WriteVTKHDF.jl

[![CI](https://github.com/Ferrite-FEM/WriteVTKHDF.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Ferrite-FEM/WriteVTKHDF.jl/actions/workflows/CI.yml)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://ferrite-fem.github.io/WriteVTKHDF.jl/dev/)

Write VTK data in the [VTKHDF file format](https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/index.html) —
the HDF5-based successor to the VTK XML formats — readable by ParaView, VisIt
and VTK. Writing only; think of it as a VTKHDF sibling of
[WriteVTK.jl](https://github.com/JuliaVTK/WriteVTK.jl), with a familiar API
(cells and data locations are shared via
[VTKBase.jl](https://github.com/JuliaVTK/VTKBase.jl)).

The headline feature over the XML formats is native **time series support in a
single file**: for a transient simulation on a fixed mesh the geometry is
stored *once*, and each time step appends only its data arrays — no more
`.pvd` collections with one `.vtu` file per step.

```julia
using WriteVTKHDF

points = rand(3, 100)                       # or vectors of SVector/tuples
cells = [MeshCell(VTKCellTypes.VTK_TETRA, [1, 2, 3, 4]), ...]

# static file
vtkhdf_grid("output", points, cells) do vtk
    vtk["temperature"] = T                   # point data (auto-detected)
    vtk["material", VTKCellData()] = mat
    vtk["velocity", VTKPointData(), attribute = :Vectors] = v  # (3, N)
end

# temporal file — geometry written once, data appended per step
vtk = vtkhdf_grid("simulation", points, cells; temporal = true)
for (t, u) in timesteps
    write_timestep(vtk, t) do frame
        frame["u"] = u
    end
end
close(vtk)
```

## Supported dataset types

| VTKHDF type | write | temporal |
|---|---|---|
| UnstructuredGrid (incl. polyhedra, multiple partitions) | ✓ | ✓ |
| PolyData (vertices/lines/polygons/strips) | ✓ | ✓ |
| ImageData | ✓ | ✓ |
| RectilinearGrid | ✓ | ✓ |
| StructuredGrid | ✓ | ✓ |
| Table | ✓ | ✓ |
| OverlappingAMR | ✓ | – |
| HyperTreeGrid | ✓ | – |
| PartitionedDataSetCollection / MultiBlockDataSet | ✓ | ✓ |

Plus gzip compression (`compress = true`), active attribute marking
(`Scalars`/`Vectors`/...), and writing datasets as blocks of composite files.
Files are validated against VTK's own `vtkHDFReader` in the test suite.

See the [documentation](https://ferrite-fem.github.io/WriteVTKHDF.jl/dev/) for
the full manual.

## Development

This package was mainly developed with AI assistance (Claude Code), guided by
the VTKHDF specification, with the file format output validated against the
VTK reference implementation.
