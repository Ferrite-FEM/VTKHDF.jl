# VTKHDF.jl

[![CI](https://github.com/Ferrite-FEM/VTKHDF.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Ferrite-FEM/VTKHDF.jl/actions/workflows/CI.yml)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://ferrite-fem.github.io/VTKHDF.jl/dev/)

Write and read VTK data in the [VTKHDF file format](https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/index.html) —
the HDF5-based successor to the VTK XML formats — readable by ParaView, VisIt
and VTK. The writing side is a VTKHDF sibling of
[WriteVTK.jl](https://github.com/JuliaVTK/WriteVTK.jl), with a familiar API
(cells and data locations are shared via
[VTKBase.jl](https://github.com/JuliaVTK/VTKBase.jl)); `vtkhdf_open` reads
the same files — and spec-conforming VTKHDF files from other writers — back
into Julia.

The headline feature over the XML formats is native **time series support in a
single file**: for a transient simulation on a fixed mesh the geometry is
stored *once*, and each time step appends only its data arrays — no more
`.pvd` collections with one `.vtu` file per step.

```julia
using VTKHDF

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

# reading files back
vtkhdf_open("output") do r
    T = r["temperature"]                     # same indexing as writing
    points, cells = read_points(r), read_cells(r)
end
vtkhdf_open("simulation") do r
    step = read_timestep(r, VTKHDF.nsteps(r))
    u_end = step["u"]
end
```

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

Plus gzip compression (`compress = true`), active attribute marking
(`Scalars`/`Vectors`/...), and writing datasets as blocks of composite files.
Files are validated against VTK's own `vtkHDFReader` in the test suite.

See the [documentation](https://ferrite-fem.github.io/VTKHDF.jl/dev/) for
the full manual.

## Development

This package was mainly developed with AI assistance (Claude Code), guided by
the VTKHDF specification, with the file format output validated against the
VTK reference implementation.
