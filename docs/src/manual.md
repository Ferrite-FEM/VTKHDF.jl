# [Manual](@id manual)

Everything builds on two ideas: geometry, and data arrays attached to it.
The manual introduces those first. [Writing](@ref writing) then covers each
dataset type; [Reading](@ref reading) does not depend on it, so skip ahead
if you only want to read existing files.

Conventions shared by everything:

- filenames without an extension get `.vtkhdf` appended,
- every writer constructor and [`vtkhdf_open`](@ref) support the do-block
  form used below, which closes the file automatically,
- data arrays use the same index syntax for writing and reading:
  `vtk["name"] = data` and `r["name"]`.

The [Examples](@ref) section contains complete runnable programs, one per
reference file of the VTKHDF specification (also in the repository's
`examples/` folder).

```@setup manual
using VTKHDF
```

## Points and cells

The simplest dataset type is the unstructured grid: point coordinates,
plus cells connecting those points. Each column of `points` is one 3D
point. The integers in a `MeshCell` are 1-based column indices:

```@example manual
points = [
    0.0 1.0 0.0 0.0
    0.0 0.0 1.0 0.0
    0.0 0.0 0.0 1.0
]
cells = [MeshCell(VTKCellTypes.VTK_TETRA, [1, 2, 3, 4])]
nothing #hide
```

Four points, one tetrahedral cell. Other dataset types describe their
geometry differently — an image grid needs no point array at all. But all
of them make the same distinction, and the rest of the API builds on it:
data is attached either to the points or to the cells.

## [Data arrays](@id data-arrays)

Every data array is attached to one of a dataset's data groups, identified
by a *location* argument:

- `VTKPointData()` — one value per grid point,
- `VTKCellData()` — one value per cell,
- `VTKFieldData()` — global arrays, independent of the grid,
- `VTKRowData()` — table columns.

Writing and reading use the same indexing. When writing, the location can
be omitted if the array size already determines it. The `attribute` keyword
marks an array as an *active attribute* (`:Scalars`, `:Vectors`, ...): the
array ParaView/VTK picks by default for coloring, glyphs etc. Using the
`points` and `cells` from above:

```@example manual
vtkhdf_grid("arrays", points, cells) do vtk
    vtk["height"] = points[3, :]          # four values: inferred as point data
    vtk["material", VTKCellData()] = [1]  # explicit location
    vtk["position", VTKPointData(), attribute = :Vectors] = points
end

vtkhdf_open("arrays") do r
    r["height"]                           # location found by name
    r["material", VTKCellData()]          # read from an explicit location
    keys(r, VTKPointData())               # array names at a location
end
nothing #hide
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

Reading returns plain arrays in this same convention, no matter which input
form was used when writing:

- scalar data as a `Vector`,
- multi-component data as a `(ncomponents, N)` matrix (`SVector` input
  comes back as its component matrix),
- image-like data with its full `([ncomponents,] nx, ny, nz)` shape,
- string field data as `Vector{String}`.

If a size matches both the points and the cells, the location must be given
explicitly. The same holds when reading a name that exists at several
locations. Names must not contain `/` or `.` (a VTKHDF format restriction).

## [Writing](@id writing)

All files are created through a small family of constructors; the dataset
type is inferred from the geometry arguments:

| dataset type | use when | created with |
|---|---|---|
| UnstructuredGrid | arbitrary cells connecting explicit points | [`vtkhdf_grid`](@ref)`(file, points, cells)` |
| PolyData | vertices/lines/polygons/strips on explicit points | `vtkhdf_grid(file, points, polys, ...)` |
| ImageData | uniformly spaced axis-aligned grid | `vtkhdf_grid(file, xrange, yrange, zrange)` |
| RectilinearGrid | axis-aligned, per-axis coordinate vectors | `vtkhdf_grid(file, xvec, yvec, zvec)` |
| StructuredGrid | curved but logically regular grid | `vtkhdf_grid(file, xyz)` |
| Table | data columns without geometry | [`vtkhdf_table`](@ref)`(file)` |
| OverlappingAMR | block-structured refinement levels | [`vtkhdf_amr`](@ref)`(file; origin)` |
| HyperTreeGrid | tree-based refinement | [`vtkhdf_htg`](@ref)`(file; dimensions)` |
| Collection / MultiBlock | several of the above in one file | [`vtkhdf_collection`](@ref) / [`vtkhdf_multiblock`](@ref) |

Every dataset constructor takes `compress = true|0-9` for gzip compression.
The grid and table constructors also take `temporal = true` for
time-dependent writing; OverlappingAMR and HyperTreeGrid are static-only.
The composite constructors (`vtkhdf_collection`/`vtkhdf_multiblock`) take
neither — both options are passed per block instead.

### Unstructured grids

Explicit points and a vector of cells, exactly as in
[Points and cells](@ref) — any mix of VTK cell types works:

```@example manual
points = [                 # unit cube corners
    0.0 1.0 1.0 0.0 0.0 1.0 1.0 0.0
    0.0 0.0 1.0 1.0 0.0 0.0 1.0 1.0
    0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0
]
cells = [
    MeshCell(VTKCellTypes.VTK_HEXAHEDRON, 1:8),
]
vtkhdf_grid("mesh", points, cells) do vtk
    vtk["height"] = points[3, :]
end
nothing #hide
```

Points can also be given as a vector of `SVector`/tuple-like values. 1‑ and
2‑dimensional points are zero-padded to 3 components.

A polyhedron is a cell described by its faces. To write them, put
[`VTKPolyhedron`](https://juliavtk.github.io/WriteVTK.jl/stable/grids/unstructured/#Polyhedron-cells)
cells in the cell vector:

```@example manual
cube = VTKPolyhedron(
    1:8,
    (1, 4, 3, 2), (1, 5, 8, 4), (5, 6, 7, 8),
    (6, 2, 3, 7), (1, 2, 6, 5), (3, 4, 8, 7),
)
vtkhdf_grid("poly", points, [cube]) do vtk   # same cube points as above
    vtk["u"] = rand(8)
end
nothing #hide
```

#### Partitions

A file can store one grid split into multiple partitions — typically one
per MPI rank of a simulation. Create the file with just the dataset tag,
without geometry. Then append each partition together with its data:

```@example manual
points1, points2 = rand(3, 4), rand(3, 4)
cells1 = [MeshCell(VTKCellTypes.VTK_TETRA, 1:4)]
cells2 = [MeshCell(VTKCellTypes.VTK_TETRA, 1:4)]

vtk = vtkhdf_grid(VTKUnstructuredGrid(), "partitioned")   # no geometry yet
add_partition(vtk, points1, cells1; pointdata = ("u" => rand(4),))
add_partition(vtk, points2, cells2; pointdata = ("u" => rand(4),))
close(vtk)
```

The same works for PolyData (`vtkhdf_grid(VTKPolyData(), ...)`), where
[`add_partition`](@ref) takes one cell vector per category:
`add_partition(vtk, points, lines, polys; ...)`.

### PolyData

PolyData describes surface-style geometry. Its cells come in four
categories: vertices, lines, polygons and triangle strips. Pass one vector
of `PolyData.*` cells per category, in any argument order:

```@example manual
points = rand(3, 4)
polys = [MeshCell(PolyData.Polys(), [1, 2, 3, 4])]
lines = [MeshCell(PolyData.Lines(), [1, 3])]
vtkhdf_grid("surface", points, polys, lines) do vtk
    vtk["height"] = rand(4)
    vtk["material", VTKCellData()] = [1, 2]   # 1 → the line, 2 → the polygon
end
nothing #hide
```

A PolyData file stores its cells in a fixed category order: Vertices,
Lines, Polygons, Strips. Cell data follows that order, not the constructor
argument order — so the first `material` value belongs to the line, even
though `polys` was passed first.

### Structured types

Choose ImageData for a uniformly spaced, axis-aligned grid; RectilinearGrid
when spacing varies independently along each axis; and StructuredGrid when
every point needs an explicit position, as in a curved grid.

```@example manual
# ImageData: ranges define a uniform origin, spacing, and extent
vtkhdf_grid("image_ranges", 0:0.1:1, 0:0.1:2, 0:0.5:3) do vtk
    vtk["u"] = rand(11, 21, 7)
end

# The equivalent geometry can be specified explicitly
vtkhdf_grid(VTKImageData(), "image_explicit", (11, 21, 7);
            origin = (0, 0, 0), spacing = (0.1, 0.1, 0.5)) do vtk
    vtk["u"] = rand(11, 21, 7)
end

# RectilinearGrid: possibly nonuniform coordinate vectors
vtkhdf_grid("rect", [0.0, 1.0, 2.5], [0.0, 2.0], [0.0, 1.0]) do vtk
    vtk["u"] = rand(3, 2, 2)
end

# StructuredGrid: arbitrary point positions as a (3, ni, nj, nk) array
xyz = rand(3, 4, 3, 2)
vtkhdf_grid("struct", xyz) do vtk
    vtk["u"] = rand(4, 3, 2)
end
nothing #hide
```

### Time series

To write a time series, create the file with `temporal = true` (grid,
table and composite-block constructors all accept it). Then add one step
at a time with [`write_timestep`](@ref):

```@example manual
points = rand(3, 4)
cells = [MeshCell(VTKCellTypes.VTK_TETRA, 1:4)]
timesteps = 0.0:0.1:1.0

vtk = vtkhdf_grid("simulation", points, cells; temporal = true)
for t in timesteps
    write_timestep(vtk, t) do frame
        frame["u"] = fill(t, 4)  # every point stores the current time
    end
end
close(vtk)
```

The geometry is written once and reused by every step (the file stores
per-step read offsets). To change the geometry, pass it to the step:

```@example manual
vtk = vtkhdf_grid("moving", points, cells; temporal = true)
for t in timesteps
    write_timestep(vtk, t; points = points .+ t, cells = cells) do frame
        frame["u"] = rand(4)
    end
end
close(vtk)
```

RectilinearGrid takes `x`/`y`/`z` the same way and StructuredGrid takes
`points`; ImageData arrays get a time dimension automatically.

Every step must write the same arrays, with the same element types, as the
first step. Adding an array, leaving one out, or changing its type in a
later step throws an error immediately.

### Table

A table has no geometry, only columns of equal length. A column is an
ordinary [data array](@ref data-arrays) with location [`VTKRowData`](@ref)
— the default for tables, so it can be left out:

```@example manual
vtkhdf_table("data") do tbl
    tbl["pressure"] = rand(100)
    tbl["id"] = collect(1:100)
end
nothing #hide
```

### Overlapping AMR

An overlapping AMR grid is a stack of refinement levels over a shared
origin. Each level has its own spacing and holds boxes of cells. A box is
given as inclusive cell-index extents:

```@example manual
ρ = rand(125)   # one value per cell of the 5×5×5 box below
vtkhdf_amr("amr"; origin = (0, 0, 0)) do amr
    lvl = add_level(amr; spacing = (1.0, 1.0, 1.0))
    add_box(lvl, (0, 4, 0, 4, 0, 4); celldata = ("ρ" => ρ,))
end
nothing #hide
```

`add_box` checks that each data array has one value per cell of the box —
here `5 × 5 × 5 = 125`, since the extents are inclusive.

### HyperTreeGrid

A HyperTreeGrid is a coarse grid whose cells can be refined recursively,
quadtree/octree style. Each coarse cell is the root of a *tree*: splitting
a cell gives `branch_factor` children per non-singleton direction, and a
child can be split again. The coordinate vectors define the coarse grid.
The interface mirrors the file format: the tree structure is passed as
flat arrays, tree by tree, level by level.

```@example manual
# 3×3×1 coordinate points => a 2×2 grid of trees:
#
#   2.0 ┌───┬───┐
#       │ 2 │ 3 │
#   1.0 ├───┼───┤      tree 0 is refined once (root + 2² children),
#       │ 0 │ 1 │      trees 1-3 stay single cells => 8 cells total
#   0.0 └───┴───┘
#      0.0  1.0  2.0
vtkhdf_htg("htg"; dimensions = (3, 3, 1), branch_factor = 2) do htg
    add_piece(
        htg;
        # one bit per cell of every non-deepest level: only tree 0's
        # root has children
        descriptors = [true],
        depth_per_tree = [2, 1, 1, 1],
        tree_ids = [0, 1, 2, 3],
        # cells per depth, tree by tree: tree 0 has 1 root + 4 children
        number_of_cells_per_tree_depth = [1, 4, 1, 1, 1],
        xcoordinates = [0.0, 1.0, 2.0],
        ycoordinates = [0.0, 1.0, 2.0],
        zcoordinates = [0.0],
        # cell data follows the same order: tree by tree, level by level
        celldata = ("depth" => Float64[0, 1, 1, 1, 1, 0, 0, 0],)
    )
end
nothing #hide
```

An optional `mask::Vector{Bool}` (same cell order) hides cells: a `true`
entry masks that cell, removing it from visualization. Refined cells
cannot be masked. See [`vtkhdf_htg`](@ref) and [`add_piece`](@ref) for all
fields.

### Composite files

A collection stores several independent datasets as named *blocks* of one
file. Blocks are created with the usual constructors, with the collection
as the first argument. Each block then accepts the same API as a standalone
file, including `temporal = true` (all temporal blocks must write the same
time values):

```@example manual
points = rand(3, 4)
cells = [MeshCell(VTKCellTypes.VTK_TETRA, 1:4)]
spoints = rand(3, 4)
polys = [MeshCell(PolyData.Polys(), [1, 2, 3, 4])]

vtkhdf_collection("multi") do col          # or vtkhdf_multiblock
    mesh = vtkhdf_grid(col, "Mesh", points, cells)
    mesh["u"] = rand(4)
    surf = vtkhdf_grid(col, "Surf", spoints, polys)
    solids = add_node(col, "solids")       # assembly hierarchy
    add_block_ref(solids, mesh)
    add_block_ref(add_node(col, "surfaces"), surf)
end
nothing #hide
```

The *assembly* built by [`add_node`](@ref)/[`add_block_ref`](@ref) is an
optional tree of named nodes referencing the blocks. ParaView shows it as
the file's grouping hierarchy.

## [Reading](@id reading)

[`vtkhdf_open`](@ref) opens a VTKHDF file for reading. The file may come
from this package or from any other spec-conforming writer, such as VTK
itself (spec major versions 1 and 2 are accepted).

```@example manual
vtkhdf_open("mesh") do r
    VTKHDF.dataset_type(r)      # "UnstructuredGrid", ...
    points = read_points(r)     # 3×N matrix
    cells = read_cells(r)       # MeshCell/VTKPolyhedron vector, 1-based ids
    height = r["height"]        # data arrays: see "Data arrays" above
    keys(r, VTKPointData())
end
nothing #hide
```

Active-attribute marks (see [Data arrays](@ref data-arrays)) are queried
with `VTKHDF.active_attributes(r, loc)`. Each dataset type also adds a few
accessors of its own. They are unexported — call them qualified, as
written:

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

```@example manual
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
collection reader. Its blocks are themselves readers and share the file
handle, so closing the collection closes everything:

```@example manual
vtkhdf_open("multi") do col
    keys(col)                        # block names
    mesh = col["Mesh"]               # a block reader
    u = mesh["u"]
    asm = VTKHDF.read_assembly(col)  # the Assembly tree
end
nothing #hide
```

### Partitions

Multi-partition UnstructuredGrid/PolyData files are returned as one grid:
`read_points` stacks the partitions, and `read_cells` renumbers the
connectivity to match. The result is usable as is, and can be written again.
`VTKHDF.partition_ranges(r)` recovers the partition structure — each
partition's index range into the point and cell data arrays.

Multi-partition PolyData has one caveat, inherited from the file format:
cells are stored per category (vertices/lines/polygons/strips), but cell
*data* is stored partition by partition. The
`VTKHDF.partition_ranges(r).cells_by_category` field maps each partition's
cells, per category, into the cell-data arrays. With a single partition the
two orders coincide and no care is needed.

```@example manual
foreach( #hide
    filename -> rm(filename; force = true), #hide
    ( #hide
        "arrays.vtkhdf", "mesh.vtkhdf", "poly.vtkhdf", "partitioned.vtkhdf", #hide
        "surface.vtkhdf", "image_ranges.vtkhdf", "image_explicit.vtkhdf", #hide
        "rect.vtkhdf", "struct.vtkhdf", #hide
        "simulation.vtkhdf", "moving.vtkhdf", "data.vtkhdf", "amr.vtkhdf", #hide
        "htg.vtkhdf", #hide
        "multi.vtkhdf", #hide
    ), #hide
) #hide
```
