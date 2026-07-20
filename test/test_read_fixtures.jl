# Hand-built HDF5 fixtures (written directly, NOT through the writer):
# spec-default fallbacks for optional datasets, legacy attribute encodings,
# and malformed files. These catch mutually-inverse writer/reader bugs that
# round trips cannot see.

function fixture(f, path; version = Int64[2, 0])
    h5open(path, "w") do file
        root = HDF5.create_group(file, "VTKHDF")
        HDF5.attrs(root)["Version"] = version
        f(root)
    end
    return path
end

# Minimal static unstructured geometry: `nparts` copies of one triangle.
function write_triangles!(root; nparts = 1)
    HDF5.attrs(root)["Type"] = "UnstructuredGrid"
    HDF5.write_dataset(root, "NumberOfPoints", fill(Int64(3), nparts))
    HDF5.write_dataset(root, "NumberOfCells", fill(Int64(1), nparts))
    HDF5.write_dataset(root, "NumberOfConnectivityIds", fill(Int64(3), nparts))
    pts = reduce(hcat, (Float64[0 1 0; 0 0 1; 0 0 0] .+ 10 * (p - 1) for p in 1:nparts))
    HDF5.write_dataset(root, "Points", pts)
    HDF5.write_dataset(root, "Connectivity", repeat(Int64[0, 1, 2], nparts))
    HDF5.write_dataset(root, "Offsets", repeat(Int64[0, 3], nparts))
    HDF5.write_dataset(root, "Types", fill(UInt8(5), nparts))
    return pts
end

@testset "fixture: temporal fallbacks (unstructured)" begin
    dir = mktempdir()
    # Two steps, one partition each; every optional Steps dataset omitted:
    # NumberOfParts/PartOffsets (constant-parts inference), PointOffsets/
    # CellOffsets (prefix-sum geometry offsets), Point/CellDataOffsets
    # (geometry fallback), FieldDataOffsets/FieldDataSizes (spec defaults).
    fn = fixture(joinpath(dir, "fallback.vtkhdf")) do root
        pts = write_triangles!(root; nparts = 2)
        pd = HDF5.create_group(root, "PointData")
        HDF5.write_dataset(pd, "u", collect(1.0:6.0))
        cd = HDF5.create_group(root, "CellData")
        HDF5.write_dataset(cd, "c", [1.0, 2.0])
        fd = HDF5.create_group(root, "FieldData")
        HDF5.write_dataset(fd, "e", [10.0, 20.0])
        HDF5.write_dataset(fd, "e2", collect(1.0:5.0))
        sg = HDF5.create_group(root, "Steps")
        HDF5.attrs(sg)["NSteps"] = Int64(2)
        HDF5.write_dataset(sg, "Values", [0.0, 1.0])
        # FieldDataSizes for e2 only: step 1 has 2 tuples, step 2 has 3
        szg = HDF5.create_group(sg, "FieldDataSizes")
        HDF5.write_dataset(szg, "e2", Int64[1 1; 2 3])
    end
    vtkhdf_open(fn) do r
        @test VTKHDF.nsteps(r) == 2
        s1, s2 = read_timestep(r, 1), read_timestep(r, 2)
        @test read_points(s1) == Float64[0 1 0; 0 0 1; 0 0 0]
        @test read_points(s2) == Float64[0 1 0; 0 0 1; 0 0 0] .+ 10.0
        @test collect(read_cells(s2)[1].connectivity) == [1, 2, 3]
        @test s1["u"] == [1.0, 2.0, 3.0]
        @test s2["u"] == [4.0, 5.0, 6.0]
        @test s1["c"] == [1.0]
        @test s2["c"] == [2.0]
        @test s1["e"] == [10.0]                  # default: 1 tuple, offset i-1
        @test s2["e"] == [20.0]
        @test s1["e2"] == [1.0, 2.0]             # sizes present, offsets cumulative
        @test s2["e2"] == [3.0, 4.0, 5.0]
    end

    # NumberOfParts missing and partitions not divisible by steps -> error
    fn = fixture(joinpath(dir, "baddiv.vtkhdf")) do root
        write_triangles!(root; nparts = 3)
        sg = HDF5.create_group(root, "Steps")
        HDF5.attrs(sg)["NSteps"] = Int64(2)
        HDF5.write_dataset(sg, "Values", [0.0, 1.0])
    end
    vtkhdf_open(fn) do r
        @test_throws ErrorException read_points(read_timestep(r, 1))
    end
end

@testset "fixture: inconsistent topology offsets are rejected" begin
    dir = mktempdir()
    # Steps/CellOffsets disagrees with the NumberOf* partition layout: the
    # reader must error rather than read the wrong cells.
    fn = fixture(joinpath(dir, "badcelloff.vtkhdf")) do root
        write_triangles!(root; nparts = 2)
        sg = HDF5.create_group(root, "Steps")
        HDF5.attrs(sg)["NSteps"] = Int64(2)
        HDF5.write_dataset(sg, "Values", [0.0, 1.0])
        HDF5.write_dataset(sg, "PartOffsets", Int64[0, 1])
        HDF5.write_dataset(sg, "NumberOfParts", Int64[1, 1])
        HDF5.write_dataset(sg, "PointOffsets", Int64[0, 3])
        HDF5.write_dataset(sg, "CellOffsets", reshape(Int64[0, 7], 1, 2))  # step 2: layout says 1
    end
    vtkhdf_open(fn) do r
        @test read_points(read_timestep(r, 1)) isa Matrix  # step 1 is consistent
        @test_throws ErrorException read_cells(read_timestep(r, 2))
    end
end

@testset "fixture: FieldDataSizes with varying components" begin
    dir = mktempdir()
    # The field dataset stores the maximum component count; per-step sizes
    # select fewer components and tuples.
    fn = fixture(joinpath(dir, "fdsizes.vtkhdf")) do root
        write_triangles!(root; nparts = 2)
        fd = HDF5.create_group(root, "FieldData")
        HDF5.write_dataset(fd, "e", Float64[1 4 7; 2 5 8; 3 6 9])  # (3, 3) tuples
        sg = HDF5.create_group(root, "Steps")
        HDF5.attrs(sg)["NSteps"] = Int64(2)
        HDF5.write_dataset(sg, "Values", [0.0, 1.0])
        szg = HDF5.create_group(sg, "FieldDataSizes")
        HDF5.write_dataset(szg, "e", Int64[3 2; 1 2])  # step 1: 3 comp × 1 tuple; step 2: 2 × 2
        og = HDF5.create_group(sg, "FieldDataOffsets")
        HDF5.write_dataset(og, "e", Int64[0, 1])
    end
    vtkhdf_open(fn) do r
        @test read_timestep(r, 1)["e"] == reshape([1.0, 2.0, 3.0], 3, 1)
        @test read_timestep(r, 2)["e"] == Float64[4 7; 5 8]  # 2 components, 2 tuples
    end
end

@testset "fixture: polydata CellDataOffsets fallback" begin
    dir = mktempdir()
    # CellOffsets present ((4, NSteps), per category) but CellDataOffsets
    # absent: the scalar cell-data offset must fall back to the column sum.
    fn = fixture(joinpath(dir, "pdfallback.vtkhdf")) do root
        HDF5.attrs(root)["Type"] = "PolyData"
        HDF5.write_dataset(root, "NumberOfPoints", Int64[3, 3])
        HDF5.write_dataset(root, "Points", Float64[0 1 0 5 6 5; 0 0 1 0 0 1; 0 0 0 0 0 0])
        pg = HDF5.create_group(root, "Polygons")
        HDF5.write_dataset(pg, "NumberOfCells", Int64[1, 1])
        HDF5.write_dataset(pg, "NumberOfConnectivityIds", Int64[3, 3])
        HDF5.write_dataset(pg, "Offsets", Int64[0, 3, 0, 3])
        HDF5.write_dataset(pg, "Connectivity", Int64[0, 1, 2, 0, 1, 2])
        cd = HDF5.create_group(root, "CellData")
        HDF5.write_dataset(cd, "c", [5.0, 6.0])
        sg = HDF5.create_group(root, "Steps")
        HDF5.attrs(sg)["NSteps"] = Int64(2)
        HDF5.write_dataset(sg, "Values", [0.0, 1.0])
        HDF5.write_dataset(sg, "PartOffsets", Int64[0, 1])
        HDF5.write_dataset(sg, "NumberOfParts", Int64[1, 1])
        HDF5.write_dataset(sg, "PointOffsets", Int64[0, 3])
        HDF5.write_dataset(sg, "CellOffsets", Int64[0 0; 0 0; 0 1; 0 0])
        HDF5.write_dataset(sg, "ConnectivityIdOffsets", Int64[0 0; 0 0; 0 3; 0 0])
    end
    vtkhdf_open(fn) do r
        s1, s2 = read_timestep(r, 1), read_timestep(r, 2)
        @test s1["c"] == [5.0]
        @test s2["c"] == [6.0]                   # offset = sum(CellOffsets[:, 2]) = 1
        @test length(read_cells(s2).polygons) == 1
        @test collect(read_cells(s2).polygons[1].connectivity) == [1, 2, 3]
    end
end

@testset "fixture: structured-kind temporal defaults" begin
    dir = mktempdir()
    # Rectilinear: no {X,Y,Z}CoordinatesOffsets -> default 0 for every step.
    fn = fixture(joinpath(dir, "rect.vtkhdf"); version = Int64[2, 7]) do root
        HDF5.attrs(root)["Type"] = "RectilinearGrid"
        HDF5.attrs(root)["WholeExtent"] = Int64[0, 1, 0, 0, 0, 0]
        HDF5.write_dataset(root, "XCoordinates", [0.0, 1.0])
        HDF5.write_dataset(root, "YCoordinates", [0.0])
        HDF5.write_dataset(root, "ZCoordinates", [0.0])
        pd = HDF5.create_group(root, "PointData")
        HDF5.write_dataset(pd, "p", reshape([1.0, 2.0, 3.0, 4.0], 2, 1, 1, 2))
        sg = HDF5.create_group(root, "Steps")
        HDF5.attrs(sg)["NSteps"] = Int64(2)
        HDF5.write_dataset(sg, "Values", [0.0, 1.0])
    end
    vtkhdf_open(fn) do r
        @test read_coordinates(read_timestep(r, 2))[1] == [0.0, 1.0]
        @test read_timestep(r, 2)["p"] == reshape([3.0, 4.0], 2, 1, 1)
    end

    # Structured: no Steps/PointOffsets -> slab 0 for every step.
    fn = fixture(joinpath(dir, "sg.vtkhdf"); version = Int64[2, 7]) do root
        HDF5.attrs(root)["Type"] = "StructuredGrid"
        HDF5.attrs(root)["WholeExtent"] = Int64[0, 1, 0, 0, 0, 0]
        HDF5.write_dataset(root, "Points", reshape(collect(1.0:6.0), 3, 2, 1, 1, 1))
        pd = HDF5.create_group(root, "PointData")
        HDF5.write_dataset(pd, "u", reshape(collect(1.0:4.0), 2, 1, 1, 2))
        sg = HDF5.create_group(root, "Steps")
        HDF5.attrs(sg)["NSteps"] = Int64(2)
        HDF5.write_dataset(sg, "Values", [0.0, 1.0])
    end
    vtkhdf_open(fn) do r
        @test read_points(read_timestep(r, 2)) == reshape(collect(1.0:6.0), 3, 2, 1, 1)
        @test read_timestep(r, 1)["u"] == reshape([1.0, 2.0], 2, 1, 1)
    end

    # a temporal array whose trailing dimension disagrees with NSteps
    fn = fixture(joinpath(dir, "badslabs.vtkhdf"); version = Int64[2, 7]) do root
        HDF5.attrs(root)["Type"] = "RectilinearGrid"
        HDF5.attrs(root)["WholeExtent"] = Int64[0, 1, 0, 0, 0, 0]
        HDF5.write_dataset(root, "XCoordinates", [0.0, 1.0])
        HDF5.write_dataset(root, "YCoordinates", [0.0])
        HDF5.write_dataset(root, "ZCoordinates", [0.0])
        pd = HDF5.create_group(root, "PointData")
        HDF5.write_dataset(pd, "p", reshape(collect(1.0:6.0), 2, 1, 1, 3))
        sg = HDF5.create_group(root, "Steps")
        HDF5.attrs(sg)["NSteps"] = Int64(2)
        HDF5.write_dataset(sg, "Values", [0.0, 1.0])
    end
    @test_throws ErrorException vtkhdf_open(fn)
end

@testset "fixture: active attribute encodings" begin
    dir = mktempdir()
    fn = fixture(joinpath(dir, "attrs.vtkhdf")) do root
        write_triangles!(root)
        pd = HDF5.create_group(root, "PointData")
        HDF5.write_dataset(pd, "u", [1.0, 2.0, 3.0])
        HDF5.write_dataset(pd, "v", rand(3, 3))
        HDF5.write_dataset(pd, "w", [4.0, 5.0, 6.0])
        HDF5.write_dataset(pd, "x", rand(3, 3))
        HDF5.attrs(pd)["Vectors"] = "v"               # group encoding (legacy + current)
        HDF5.attrs(pd["w"])["Attribute"] = "scalars"  # array-only, case-insensitive
        HDF5.attrs(pd["x"])["Attribute"] = "Vectors"  # conflicts with the group
    end
    vtkhdf_open(fn) do r
        active = @test_logs (:warn, r"conflicting active-attribute") begin
            VTKHDF.active_attributes(r, VTKPointData())
        end
        @test active == Dict(:Vectors => "v", :Scalars => "w")
        @test VTKHDF.data_attributes(r, "w", VTKPointData()) == [:Scalars]
    end

    # duplicate array-only claims: first name in sorted order wins
    fn = fixture(joinpath(dir, "attrs2.vtkhdf")) do root
        write_triangles!(root)
        pd = HDF5.create_group(root, "PointData")
        HDF5.write_dataset(pd, "a", [1.0, 2.0, 3.0])
        HDF5.write_dataset(pd, "b", [4.0, 5.0, 6.0])
        HDF5.attrs(pd["a"])["Attribute"] = "Scalars"
        HDF5.attrs(pd["b"])["Attribute"] = "Scalars"
    end
    vtkhdf_open(fn) do r
        active = @test_logs (:warn, r"conflicting active-attribute") begin
            VTKHDF.active_attributes(r, VTKPointData())
        end
        @test active == Dict(:Scalars => "a")
    end
end

@testset "fixture: PDC index order" begin
    dir = mktempdir()
    function table_block!(root, name, index)
        grp = HDF5.create_group(root, name)
        HDF5.attrs(grp)["Type"] = "Table"
        HDF5.attrs(grp)["Version"] = Int64[2, 8]
        HDF5.attrs(grp)["Index"] = Int64(index)
        HDF5.write_dataset(grp, "NumberOfRows", Int64[1])
        rd = HDF5.create_group(grp, "RowData")
        HDF5.write_dataset(rd, "x", [Float64(index)])
        return grp
    end
    fn = fixture(joinpath(dir, "pdc.vtkhdf"); version = Int64[2, 1]) do root
        HDF5.attrs(root)["Type"] = "PartitionedDataSetCollection"
        table_block!(root, "b", 1)   # created first, but Index 1
        table_block!(root, "a", 0)
    end
    vtkhdf_open(fn) do col
        @test keys(col) == ["a", "b"]  # ordered by Index, not creation/name order
        @test col["a"]["x"] == [0.0]
        @test col["b"]["x"] == [1.0]
    end

    fn = fixture(joinpath(dir, "dup.vtkhdf"); version = Int64[2, 1]) do root
        HDF5.attrs(root)["Type"] = "PartitionedDataSetCollection"
        table_block!(root, "a", 0)
        table_block!(root, "b", 0)   # duplicate Index
    end
    @test_throws ErrorException vtkhdf_open(fn)
end

@testset "fixture: malformed files" begin
    dir = mktempdir()

    # not an HDF5 file / nonexistent file
    txt = joinpath(dir, "plain.vtkhdf")
    write(txt, "not hdf5")
    @test_throws ArgumentError vtkhdf_open(txt)
    @test_throws ArgumentError vtkhdf_open(joinpath(dir, "nope.vtkhdf"))

    # HDF5 file without /VTKHDF
    fn = joinpath(dir, "novtk.vtkhdf")
    h5open(fn, "w") do file
        HDF5.write_dataset(file, "x", [1.0])
    end
    @test_throws ErrorException vtkhdf_open(fn)

    # unknown dataset type
    fn = fixture(joinpath(dir, "unknown.vtkhdf")) do root
        HDF5.attrs(root)["Type"] = "SomethingElse"
    end
    @test_throws ErrorException vtkhdf_open(fn)

    # unsupported major version
    fn = fixture(joinpath(dir, "v3.vtkhdf"); version = Int64[3, 0]) do root
        HDF5.attrs(root)["Type"] = "UnstructuredGrid"
    end
    @test_throws ErrorException vtkhdf_open(fn)

    # missing Type / missing Version
    fn = fixture(joinpath(dir, "notype.vtkhdf")) do root
    end
    @test_throws ErrorException vtkhdf_open(fn)
    fn = joinpath(dir, "nover.vtkhdf")
    h5open(fn, "w") do file
        root = HDF5.create_group(file, "VTKHDF")
        HDF5.attrs(root)["Type"] = "UnstructuredGrid"
    end
    @test_throws ErrorException vtkhdf_open(fn)

    # Steps/Values length mismatch
    fn = fixture(joinpath(dir, "badvals.vtkhdf")) do root
        write_triangles!(root)
        sg = HDF5.create_group(root, "Steps")
        HDF5.attrs(sg)["NSteps"] = Int64(2)
        HDF5.write_dataset(sg, "Values", [0.0])
    end
    @test_throws ErrorException vtkhdf_open(fn)

    # a Steps offsets dataset with the wrong length
    fn = fixture(joinpath(dir, "badoff.vtkhdf")) do root
        write_triangles!(root; nparts = 2)
        sg = HDF5.create_group(root, "Steps")
        HDF5.attrs(sg)["NSteps"] = Int64(2)
        HDF5.write_dataset(sg, "Values", [0.0, 1.0])
        HDF5.write_dataset(sg, "PointOffsets", Int64[0])
    end
    @test_throws ErrorException vtkhdf_open(fn)

    # unknown cell type id
    fn = fixture(joinpath(dir, "badcell.vtkhdf")) do root
        write_triangles!(root)
        HDF5.delete_object(root, "Types")
        HDF5.write_dataset(root, "Types", [UInt8(199)])
    end
    vtkhdf_open(fn) do r
        @test_throws ErrorException read_cells(r)
    end

    # polyhedron cell type without the polyhedron datasets
    fn = fixture(joinpath(dir, "nopoly.vtkhdf")) do root
        write_triangles!(root)
        HDF5.delete_object(root, "Types")
        HDF5.write_dataset(root, "Types", [UInt8(42)])
    end
    vtkhdf_open(fn) do r
        @test_throws ErrorException read_cells(r)
    end

    # counts inconsistent with the dataset sizes
    fn = fixture(joinpath(dir, "badcount.vtkhdf")) do root
        write_triangles!(root)
        HDF5.delete_object(root, "NumberOfPoints")
        HDF5.write_dataset(root, "NumberOfPoints", Int64[5])
    end
    @test_throws ErrorException vtkhdf_open(fn)

    # ambiguous auto-location
    fn = fixture(joinpath(dir, "ambig.vtkhdf")) do root
        write_triangles!(root)
        pd = HDF5.create_group(root, "PointData")
        HDF5.write_dataset(pd, "u", [1.0, 2.0, 3.0])
        cd = HDF5.create_group(root, "CellData")
        HDF5.write_dataset(cd, "u", [9.0])
    end
    vtkhdf_open(fn) do r
        @test_throws ErrorException r["u"]
        @test r["u", VTKPointData()] == [1.0, 2.0, 3.0]
        @test r["u", VTKCellData()] == [9.0]
    end

    # temporal AMR/HTG are rejected
    fn = fixture(joinpath(dir, "tamr.vtkhdf"); version = Int64[2, 3]) do root
        HDF5.attrs(root)["Type"] = "OverlappingAMR"
        HDF5.attrs(root)["Origin"] = Float64[0, 0, 0]
        sg = HDF5.create_group(root, "Steps")
        HDF5.attrs(sg)["NSteps"] = Int64(1)
        HDF5.write_dataset(sg, "Values", [0.0])
    end
    @test_throws ErrorException vtkhdf_open(fn)
end
