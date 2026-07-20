# [Manual](@id manual)

The package writes and reads all VTKHDF dataset types. A few conventions are
shared by everything: filenames without an extension get `.vtkhdf` appended,
every writer constructor and [`vtkhdf_open`](@ref) support the do-block form
(closing the file automatically), and data arrays use the same index syntax
and shapes in both directions.

Complete runnable programs, ported from the reference files shown in the
VTKHDF specification, are in the [Examples](@ref) section (and in the
`examples/` folder of the repository).

## [Data arrays](@id data-arrays)

Every data array is attached to one of a dataset's data groups, identified
by a *location* argument: `VTKPointData()` (one value per grid point),
`VTKCellData()` (one value per cell), `VTKFieldData()` (global arrays,
independent of the grid), or `VTKRowData()` (table columns). An array can
additionally be marked as an *active attribute* (`:Scalars`, `:Vectors`,
...) — the array ParaView/VTK picks by default for coloring, glyphs etc.
Writing and reading use the same indexing:

```julia
vtk["name"] = data                       # write; location inferred from the size
vtk["name", VTKCellData()] = data        # write with explicit location
vtk["v", VTKPointData(), attribute = :Vectors] = v   # mark active attribute

r["name"]                                # read; location found by name
r["name", VTKCellData()]                 # read from an explicit location
keys(r, VTKPointData())                  # array names at a location
```

Accepted array shapes when writing follow WriteVTK.jl's component-first
convention:

- vectors of length `N` for scalar data,
- `(ncomponents, N)` matrices for vector/tensor data,
- vectors of `SVector`/`NTuple`/`Tensors.Vec`-like isbits elements
  (stored as their components),
- for ImageData/RectilinearGrid/StructuredGrid: arrays shaped like the grid,
  `(nx, ny, nz)` or `(ncomponents, nx, ny, nz)` (trailing singleton
  dimensions can be dropped),
- strings and vectors of strings as field data (static files only).

Reading returns the values unchanged, as plain arrays in the same
convention regardless of which input form was used when writing: scalar
data as a `Vector`, multi-component data as a `(ncomponents, N)` matrix
(so `SVector` input comes back as its component matrix), image-like data
with its full `([ncomponents,] nx, ny, nz)` shape, and string field data as
`Vector{String}`.

If a size matches both the points and the cells (or, when reading, a name
exists at several locations), the location must be given explicitly. Names
must not contain `/` or `.` (a VTKHDF format restriction).

## Writing

All files are created through a small family of constructors. Every dataset
constructor takes `compress = true|0-9` for gzip compression, and the grid
and table constructors take `temporal = true` for time-dependent writing
(OverlappingAMR and HyperTreeGrid are static-only). The composite
constructors (`vtkhdf_collection`/`vtkhdf_multiblock`) take neither — those
options are passed per block instead.

### Unstructured grids

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

#### Partitions

VTKHDF files can store multiple partitions (as produced by e.g. one MPI rank
each). Creating the file with just the dataset tag defers the geometry;
partitions are then appended explicitly, with their data:

```julia
vtk = vtkhdf_grid(VTKUnstructuredGrid(), "partitioned")   # no geometry yet
add_partition(vtk, points₁, cells₁; pointdata = ("u" => u₁,))
add_partition(vtk, points₂, cells₂; pointdata = ("u" => u₂,))
close(vtk)
```

The same works for PolyData (`vtkhdf_grid(VTKPolyData(), ...)`), where
[`add_partition`](@ref) takes one cell vector per category:
`add_partition(vtk, points, lines, polys; ...)`.

### PolyData

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

### Structured types

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

### Time series

Opening a grid, table or composite-block constructor with `temporal = true`
enables [`write_timestep`](@ref):

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

### Table

```julia
vtkhdf_table("data") do tbl
    tbl["pressure"] = rand(100)
    tbl["id"] = collect(1:100)
end
```

### Overlapping AMR

```julia
vtkhdf_amr("amr"; origin = (0, 0, 0)) do amr
    lvl = add_level(amr; spacing = (1.0, 1.0, 1.0))
    add_box(lvl, (0, 4, 0, 4, 0, 4); celldata = ("ρ" => ρ,))
end
```

Boxes are inclusive cell-index extents; data sizes are validated against them.

### HyperTreeGrid

A low-level interface following the file format directly; see
[`vtkhdf_htg`](@ref) and [`add_piece`](@ref) for the field descriptions.

### Composite files

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

The *assembly* is an optional tree of named nodes referencing the blocks —
the grouping hierarchy ParaView shows for the file. Blocks accept the same
API as standalone files, including `temporal = true` (all temporal blocks
must write the same time values).

## Reading

[`vtkhdf_open`](@ref) opens a VTKHDF file for reading — files written by
this package or by any other spec-conforming writer such as VTK itself
(spec major versions 1 and 2 are accepted).

```julia
vtkhdf_open("mesh") do r
    VTKHDF.dataset_type(r)      # "UnstructuredGrid", ...
    points = read_points(r)     # 3×N matrix
    cells = read_cells(r)       # MeshCell/VTKPolyhedron vector, 1-based ids
    u = r["u"]                  # data arrays: see "Data arrays" above
    keys(r, VTKPointData())
end
```

Active-attribute marks (see [Data arrays](@ref data-arrays)) are queried
with `VTKHDF.active_attributes(r, loc)`. Beyond the common API, each
dataset type has a few specific accessors (written as shown — qualified
names are unexported):

| dataset type | accessors |
|---|---|
| UnstructuredGrid, PolyData, ImageData, RectilinearGrid, StructuredGrid | `VTKHDF.npoints`, `VTKHDF.ncells` |
| UnstructuredGrid, PolyData | `VTKHDF.npartitions`, `VTKHDF.partition_ranges` |
| ImageData, RectilinearGrid, StructuredGrid | `VTKHDF.grid_info` (the stored dims/origin/spacing/... attributes) |
| RectilinearGrid | [`read_coordinates`](@ref) |
| Table | `VTKHDF.nrows` |
| OverlappingAMR | `VTKHDF.grid_info`, `VTKHDF.nlevels`, `VTKHDF.amr_level` — levels support `VTKHDF.level_info`, `npoints`/`ncells`/`partition_ranges` and data indexing |
| HyperTreeGrid | `VTKHDF.grid_info`, `VTKHDF.ncells`, `VTKHDF.npieces`, `VTKHDF.htg_piece` |

### Time series

For a temporal file ([`VTKHDF.is_temporal`](@ref)), geometry and data are
read through step views:

```julia
vtkhdf_open("simulation") do r
    for i in 1:VTKHDF.nsteps(r)
        step = read_timestep(r, i)
        t = VTKHDF.time_value(step)
        u = step["u"]
        points = read_points(step)   # follows per-step geometry changes
    end
end
```

### Composite files

Opening a `PartitionedDataSetCollection`/`MultiBlockDataSet` returns a
collection reader; blocks are themselves readers and share the file handle
(closing the collection closes everything):

```julia
vtkhdf_open("multi") do col
    keys(col)                        # block names
    mesh = col["Mesh"]               # a block reader
    u = mesh["u"]
    asm = VTKHDF.read_assembly(col)  # the Assembly tree
end
```

### Partitions

Multi-partition UnstructuredGrid/PolyData files are returned as one grid:
`read_points` stacks the partitions and `read_cells` renumbers connectivity
to match, so the result is directly usable (and writable again). To recover
the partition structure, `VTKHDF.partition_ranges(r)` gives each partition's
index range into the point and cell data arrays.

One caveat is inherited from the file format, for multi-partition PolyData
only: cells are stored per category (vertices/lines/polygons/strips) but
cell *data* is stored partition by partition. The
`VTKHDF.partition_ranges(r).cells_by_category` field maps each partition's
cells per category into the cell-data arrays. With a single partition the
two orders coincide and no care is needed.
