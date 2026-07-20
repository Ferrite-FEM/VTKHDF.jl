# Reading back UnstructuredGrid files (static and temporal).

@testset "read unstructured static" begin
    dir = mktempdir()
    points1 = Float64[0 1 0 0; 0 0 1 0; 0 0 0 1]
    cells1 = [MeshCell(VTKCellTypes.VTK_TETRA, [1, 2, 3, 4])]
    points2 = Float64[2 3 2 2 3 3 2 3; 0 0 1 0 1 0 1 1; 0 0 0 1 0 1 1 1]
    hexfaces = ((1, 3, 4, 2), (1, 2, 6, 5), (1, 5, 7, 3), (8, 6, 2, 4), (8, 4, 3, 7), (8, 7, 5, 6))
    cells2 = [VTKPolyhedron(1:8, hexfaces...)]
    v = rand(3, 12)

    fn = joinpath(dir, "ug")
    vtkhdf_grid(fn, points1, cells1) do vtk
        add_partition(vtk, points2, Union{MeshCell{VTKCellTypes.VTKCellType}, VTKPolyhedron}[c for c in cells2])
        vtk["u"] = collect(1.0:12.0)
        vtk["mat", VTKCellData()] = Int32[7, 9]
        vtk["v", VTKPointData(), attribute = :Vectors] = v
        vtk["note", VTKFieldData()] = "hello"
        vtk["params", VTKFieldData()] = [1.5, 2.5]
    end

    r = vtkhdf_open(fn * ".vtkhdf")
    @test r isa VTKHDF.VTKHDFReader
    @test sprint(show, r) == "VTKHDFReader{ReadUnstructured} (open)"
    @test VTKHDF.dataset_type(r) == "UnstructuredGrid"
    @test VTKHDF.file_version(r) == (2, 6)  # bumped by the attribute marking
    @test !VTKHDF.is_temporal(r)
    @test_throws ErrorException VTKHDF.nsteps(r)
    @test VTKHDF.npoints(r) == 12
    @test VTKHDF.ncells(r) == 2
    @test VTKHDF.npartitions(r) == 2
    @test read_points(r) == hcat(points1, points2)

    cells = read_cells(r)
    @test length(cells) == 2
    @test cells[1] isa MeshCell
    @test cells[1].ctype == VTKCellTypes.VTK_TETRA
    @test collect(cells[1].connectivity) == [1, 2, 3, 4]
    # the polyhedron of partition 2 is rebased to global point ids
    @test cells[2] isa VTKPolyhedron
    @test collect(cells[2].connectivity) == collect(5:12)
    @test map(f -> Tuple(f .+ 4), hexfaces) == Tuple(map(Tuple, VTKHDF.VTKBase.faces(cells[2])))

    @test r["u"] == collect(1.0:12.0)
    @test r["u"] isa Vector{Float64}
    @test r["mat"] == Int32[7, 9]
    @test r["mat"] isa Vector{Int32}
    @test r["v", VTKPointData()] == v
    @test r["note", VTKFieldData()] == ["hello"]
    @test r["params"] == [1.5, 2.5]
    @test keys(r, VTKPointData()) == ["u", "v"]
    @test keys(r, VTKCellData()) == ["mat"]
    @test sort(keys(r, VTKFieldData())) == ["note", "params"]
    @test haskey(r, "u", VTKPointData())
    @test !haskey(r, "u", VTKCellData())
    @test_throws ErrorException r["nope"]

    @test VTKHDF.active_attributes(r, VTKPointData()) == Dict(:Vectors => "v")
    @test VTKHDF.data_attributes(r, "v", VTKPointData()) == [:Vectors]
    @test VTKHDF.data_attributes(r, "u", VTKPointData()) == Symbol[]
    @test VTKHDF.active_attributes(r, VTKCellData()) == Dict{Symbol, String}()

    pr = VTKHDF.partition_ranges(r)
    @test pr.points == [1:4, 5:12]
    @test pr.cells == [1:1, 2:2]

    close(r)
    @test !isopen(r)
    @test_throws ErrorException read_points(r)
    @test_throws ErrorException r["u"]
    close(r)  # second close is a no-op

    # round trip through the writer again
    r = vtkhdf_open(fn)
    pts = read_points(r)
    cls = read_cells(r)
    u = r["u"]
    close(r)
    fn2 = joinpath(dir, "ug2")
    vtkhdf_grid(fn2, pts, cls) do vtk
        vtk["u"] = u
    end
    vtkhdf_open(fn2) do r2
        @test read_points(r2) == pts
        @test r2["u"] == u
    end
end

@testset "read unstructured static edge cases" begin
    dir = mktempdir()

    # empty file (no partitions written)
    fn = joinpath(dir, "empty")
    vtkhdf_grid(VTKUnstructuredGrid(), fn) do vtk
    end
    vtkhdf_open(fn) do r
        @test VTKHDF.npoints(r) == 0
        @test VTKHDF.ncells(r) == 0
        @test size(read_points(r)) == (3, 0)
        @test isempty(read_cells(r))
    end

    # compression and Float32 points
    fn = joinpath(dir, "comp")
    pts = rand(Float32, 3, 10)
    cells = [MeshCell(VTKCellTypes.VTK_VERTEX, [i]) for i in 1:10]
    vtkhdf_grid(fn, pts, cells; compress = true) do vtk
        vtk["u", VTKPointData()] = collect(Float32, 1:10)
    end
    vtkhdf_open(fn) do r
        @test read_points(r) == pts
        @test read_points(r) isa Matrix{Float32}
        @test r["u"] isa Vector{Float32}
        @test length(read_cells(r)) == 10
    end
end

@testset "read unstructured temporal" begin
    dir = mktempdir()
    pts_a = Float64[0 1 0; 0 0 1; 0 0 0]
    cells_a = [MeshCell(VTKCellTypes.VTK_TRIANGLE, [1, 2, 3])]
    pts_b = Float64[0 2 0 0; 0 0 2 0; 0 0 0 2]
    cells_b = [MeshCell(VTKCellTypes.VTK_TETRA, [1, 2, 3, 4])]

    fn = joinpath(dir, "tug")
    vtk = vtkhdf_grid(fn, pts_a, cells_a; temporal = true)
    write_timestep(vtk, 0.0) do f
        f["u"] = [1.0, 2.0, 3.0]
        f["c", VTKCellData()] = [10.0]
        f["e", VTKFieldData()] = [10.0]
    end
    write_timestep(vtk, 0.5) do f  # geometry reused
        f["u"] = [4.0, 5.0, 6.0]
        f["c", VTKCellData()] = [20.0]
        f["e", VTKFieldData()] = [20.0]
    end
    write_timestep(vtk, 1.0; points = pts_b, cells = cells_b) do f
        f["u"] = [7.0, 8.0, 9.0, 10.0]
        f["c", VTKCellData()] = [30.0]
        f["e", VTKFieldData()] = [30.0]
    end
    close(vtk)

    r = vtkhdf_open(fn)
    @test VTKHDF.is_temporal(r)
    @test VTKHDF.nsteps(r) == 3
    @test VTKHDF.time_values(r) == [0.0, 0.5, 1.0]
    @test occursin("temporal, 3 steps", sprint(show, r))
    # static access must point to read_timestep
    @test_throws ErrorException read_points(r)
    @test_throws ErrorException read_cells(r)
    @test_throws ErrorException r["u"]
    @test_throws ErrorException VTKHDF.npoints(r)
    @test_throws BoundsError read_timestep(r, 0)
    @test_throws BoundsError read_timestep(r, 4)

    s1, s2, s3 = (read_timestep(r, i) for i in 1:3)
    @test VTKHDF.time_value(s2) == 0.5
    @test occursin("2/3", sprint(show, s2))
    @test read_points(s1) == pts_a
    @test read_points(s2) == pts_a  # reused geometry reads identically
    @test read_points(s3) == pts_b
    @test read_cells(s1)[1].ctype == VTKCellTypes.VTK_TRIANGLE
    @test read_cells(s2)[1].ctype == VTKCellTypes.VTK_TRIANGLE
    @test read_cells(s3)[1].ctype == VTKCellTypes.VTK_TETRA
    @test collect(read_cells(s3)[1].connectivity) == [1, 2, 3, 4]
    @test s1["u"] == [1.0, 2.0, 3.0]
    @test s2["u"] == [4.0, 5.0, 6.0]
    @test s3["u"] == [7.0, 8.0, 9.0, 10.0]
    @test s1["c"] == [10.0]
    @test s3["c", VTKCellData()] == [30.0]
    @test s1["e"] == [10.0]
    @test s3["e", VTKFieldData()] == [30.0]
    @test VTKHDF.npoints(s3) == 4
    @test VTKHDF.ncells(s2) == 1
    @test VTKHDF.npartitions(s1) == 1
    @test VTKHDF.partition_ranges(s3).points == [1:4]
    close(r)

    # zero-step temporal file
    fn = joinpath(dir, "zero")
    @test_logs (:warn, r"without any time steps") close(
        vtkhdf_grid(fn, pts_a, cells_a; temporal = true)
    )
    vtkhdf_open(fn) do r
        @test VTKHDF.nsteps(r) == 0
        @test VTKHDF.time_values(r) == Float64[]
        @test_throws BoundsError read_timestep(r, 1)
    end
end

@testset "read unstructured temporal polyhedra" begin
    dir = mktempdir()
    hexfaces = ((1, 3, 4, 2), (1, 2, 6, 5), (1, 5, 7, 3), (8, 6, 2, 4), (8, 4, 3, 7), (8, 7, 5, 6))
    pts1 = Float64[0 1 0 0 0 1 0 1; 0 0 1 0 1 0 1 1; 0 0 0 1 1 1 1 1]
    pts2 = pts1 .+ 5
    poly = VTKPolyhedron(1:8, hexfaces...)

    fn = joinpath(dir, "tpoly")
    vtk = vtkhdf_grid(fn, pts1, [poly]; temporal = true)
    write_timestep(f -> f["u"] = collect(1.0:8.0), vtk, 0.0)
    write_timestep(f -> f["u"] = collect(9.0:16.0), vtk, 1.0; points = pts2, cells = [poly])
    write_timestep(f -> f["u"] = collect(17.0:24.0), vtk, 2.0)  # reuse
    close(vtk)

    vtkhdf_open(fn) do r
        for (i, expected_pts) in ((1, pts1), (2, pts2), (3, pts2))
            s = read_timestep(r, i)
            @test read_points(s) == expected_pts
            cells = read_cells(s)
            @test length(cells) == 1
            @test cells[1] isa VTKPolyhedron
            @test collect(cells[1].connectivity) == collect(1:8)
            @test Tuple(map(Tuple, VTKHDF.VTKBase.faces(cells[1]))) == hexfaces
        end
        @test read_timestep(r, 3)["u"] == collect(17.0:24.0)
    end
end

@testset "read unstructured temporal multi-partition" begin
    dir = mktempdir()
    p1 = Float64[0 1 0; 0 0 1; 0 0 0]
    c1 = [MeshCell(VTKCellTypes.VTK_TRIANGLE, [1, 2, 3])]
    p2 = Float64[5 6 5; 0 0 1; 0 0 0]
    c2 = [MeshCell(VTKCellTypes.VTK_TRIANGLE, [1, 2, 3])]

    fn = joinpath(dir, "mp")
    vtk = vtkhdf_grid(VTKUnstructuredGrid(), fn; temporal = true)
    write_timestep(vtk, 0.0) do f
        add_partition(f, p1, c1)
        add_partition(f, p2, c2)
        f["u"] = collect(1.0:6.0)
    end
    write_timestep(vtk, 1.0) do f
        f["u"] = collect(7.0:12.0)
    end
    close(vtk)

    vtkhdf_open(fn) do r
        s1 = read_timestep(r, 1)
        s2 = read_timestep(r, 2)
        @test VTKHDF.npartitions(s1) == 2
        @test read_points(s1) == hcat(p1, p2)
        @test read_points(s2) == hcat(p1, p2)
        cells = read_cells(s1)
        @test length(cells) == 2
        @test collect(cells[1].connectivity) == [1, 2, 3]
        @test collect(cells[2].connectivity) == [4, 5, 6]  # rebased into step points
        @test s1["u"] == collect(1.0:6.0)
        @test s2["u"] == collect(7.0:12.0)
        @test VTKHDF.partition_ranges(s2) == (points = [1:3, 4:6], cells = [1:1, 2:2])
    end
end
