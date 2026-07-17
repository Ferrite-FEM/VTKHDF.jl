# [Manual](@id manual)

All files are created through a small family of constructors. Every
constructor supports the do-block form (closing the file automatically), a
`temporal = true` keyword for time-dependent writing, `compress = true|0-9`
for gzip compression, and appends the `.vtkhdf` extension when the filename
has none.

## Data arrays

Data is attached with index syntax:

```julia
vtk["name"] = data                       # location inferred from the size
vtk["name", VTKPointData()] = data       # explicit location
vtk["name", VTKCellData()] = data
vtk["name", VTKFieldData()] = [1.0, 2.0] # global arrays (also strings)
vtk["v", VTKPointData(), attribute = :Vectors] = v   # mark active attribute
```

Accepted array shapes follow WriteVTK.jl's component-first convention:

- vectors of length `N` for scalar data,
- `(ncomponents, N)` matrices for vector/tensor data,
- vectors of `SVector`/`NTuple`/`Tensors.Vec`-like isbits elements,
- for ImageData/RectilinearGrid/StructuredGrid: arrays shaped like the grid,
  `(nx, ny, nz)` or `(ncomponents, nx, ny, nz)` (trailing singleton
  dimensions can be dropped).

If a size matches both the points and the cells, the location must be given
explicitly. Names must be ASCII without `/` or `.` (a VTKHDF format
restriction).

## Unstructured grids

```julia
points = rand(3, 8)        # 3×N matrix; N-vectors of point-like values also work
cells = [
    MeshCell(VTKCellTypes.VTK_HEXAHEDRON, 1:8),
]
vtkhdf_grid("mesh", points, cells) do vtk
    vtk["u"] = rand(8)
end
```

1‑ and 2‑dimensional points are zero-padded to 3 components. Polyhedra are
written by putting [`VTKPolyhedron`](https://juliavtk.github.io/WriteVTK.jl/stable/grids/unstructured/#Polyhedron-cells)
cells in the cell vector:

```julia
cube = VTKPolyhedron(
    1:8,
    (1, 4, 3, 2), (1, 5, 8, 4), (5, 6, 7, 8),
    (6, 2, 3, 7), (1, 2, 6, 5), (3, 4, 8, 7),
)
vtkhdf_grid("poly", points, [cube]) do vtk
    vtk["u"] = rand(8)
end
```

### Partitions

VTKHDF files can store multiple partitions (as produced by e.g. one MPI rank
each). Partitions are appended explicitly, with their data:

```julia
vtk = vtkhdf_grid(VTKUnstructuredGrid(), "partitioned")   # deferred geometry
add_partition(vtk, points₁, cells₁; pointdata = ("u" => u₁,))
add_partition(vtk, points₂, cells₂; pointdata = ("u" => u₂,))
close(vtk)
```

## PolyData

Cells use the `PolyData.*` cell types; pass one homogeneous vector per
category. On disk the categories are ordered Vertices, Lines, Polygons,
Strips — cell data must be supplied in that concatenated order.

```julia
polys = [MeshCell(PolyData.Polys(), [1, 2, 3, 4])]
lines = [MeshCell(PolyData.Lines(), [1, 3])]
vtkhdf_grid("surface", points, polys, lines) do vtk
    vtk["height"] = rand(4)
end
```

## Structured types

```julia
# ImageData: ranges define origin/spacing/extent
vtkhdf_grid("image", 0:0.1:1, 0:0.1:2, 0:0.5:3) do vtk
    vtk["u"] = rand(11, 21, 7)
end
# ... or explicitly
vtkhdf_grid(VTKImageData(), "image", (11, 21, 7); origin = (0, 0, 0), spacing = (0.1, 0.1, 0.5))

# RectilinearGrid: coordinate vectors
vtkhdf_grid("rect", [0.0, 1.0, 2.5], [0.0, 2.0], [0.0, 1.0])

# StructuredGrid: explicit point positions, (3, ni, nj, nk)
vtkhdf_grid("struct", xyz)
```

## Time series

Opening any constructor with `temporal = true` enables
[`write_timestep`](@ref):

```julia
vtk = vtkhdf_grid("simulation", points, cells; temporal = true)
for (t, u) in timesteps
    write_timestep(vtk, t) do frame
        frame["u"] = u
    end
end
close(vtk)
```

The geometry is written once and reused by every step (the file stores
per-step read offsets). To change the geometry, pass it to the step:

```julia
write_timestep(vtk, t; points = new_points, cells = new_cells) do frame
    frame["u"] = u
end
```

(similarly `x`/`y`/`z` for RectilinearGrid and `points` for StructuredGrid;
ImageData arrays get a time dimension automatically.)

The set of arrays and their types must be identical in every step — it is
frozen by the first step, and violations throw immediately.

## Table

```julia
vtkhdf_table("data") do tbl
    tbl["pressure"] = rand(100)
    tbl["id"] = collect(1:100)
end
```

## Overlapping AMR

```julia
vtkhdf_amr("amr"; origin = (0, 0, 0)) do amr
    lvl = add_level(amr; spacing = (1.0, 1.0, 1.0))
    add_box(lvl, (0, 4, 0, 4, 0, 4); celldata = ("ρ" => ρ,))
end
```

Boxes are inclusive cell-index extents; data sizes are validated against them.

## HyperTreeGrid

A low-level interface following the file format directly; see
[`vtkhdf_htg`](@ref) and [`add_piece`](@ref) for the field descriptions.

## Composite files

```julia
vtkhdf_collection("multi") do col          # or vtkhdf_multiblock
    mesh = vtkhdf_grid(col, "Mesh", points, cells)
    mesh["u"] = u
    surf = vtkhdf_grid(col, "Surf", spoints, polys)
    solids = add_node(col, "solids")       # assembly hierarchy
    add_block_ref(solids, mesh)
    add_block_ref(add_node(col, "surfaces"), surf)
end
```

Blocks accept the same API as standalone files, including `temporal = true`
(all temporal blocks must write the same time values).
