# Reading back PolyData files, including the partition-major cell-data
# ordering exposed by partition_ranges.

@testset "read polydata static multi-partition" begin
    dir = mktempdir()
    # Adversarial layout: two partitions, several nonempty categories per
    # partition, distinct cell-data markers to pin down the on-disk order.
    p1 = Float64[0 1 2; 0 1 2; 0 0 0]
    verts1 = [MeshCell(PolyData.Verts(), [1])]
    lines1 = [MeshCell(PolyData.Lines(), [1, 2]), MeshCell(PolyData.Lines(), [2, 3])]
    p2 = Float64[5 6 7 8; 0 0 1 1; 0 0 0 0]
    lines2 = [MeshCell(PolyData.Lines(), [1, 2])]
    polys2 = [MeshCell(PolyData.Polys(), [1, 2, 3, 4])]

    fn = joinpath(dir, "pd")
    vtkhdf_grid(fn, p1, verts1, lines1) do vtk
        add_partition(vtk, p2, lines2, polys2)
        # partition-major: [vert1, line1, line1, line2, poly2]
        vtk["marker", VTKCellData()] = Float64[1, 2, 3, 4, 5]
        vtk["u"] = collect(1.0:7.0)
    end

    r = vtkhdf_open(fn)
    @test VTKHDF.dataset_type(r) == "PolyData"
    @test VTKHDF.npoints(r) == 7
    @test VTKHDF.ncells(r) == 5
    @test VTKHDF.npartitions(r) == 2
    @test read_points(r) == hcat(p1, p2)

    cells = read_cells(r)
    @test cells isa NamedTuple{(:vertices, :lines, :polygons, :strips)}
    @test length(cells.vertices) == 1
    @test length(cells.lines) == 3
    @test length(cells.polygons) == 1
    @test isempty(cells.strips)
    @test collect(cells.vertices[1].connectivity) == [1]
    @test collect(cells.lines[1].connectivity) == [1, 2]
    @test collect(cells.lines[2].connectivity) == [2, 3]
    @test collect(cells.lines[3].connectivity) == [4, 5]        # partition 2, rebased
    @test collect(cells.polygons[1].connectivity) == [4, 5, 6, 7]

    # partition-major cell data with per-partition category ranges
    marker = r["marker", VTKCellData()]
    pr = VTKHDF.partition_ranges(r)
    @test pr.points == [1:3, 4:7]
    @test pr.cells == [1:3, 4:5]
    @test pr.cells_by_category[1].vertices == 1:1
    @test pr.cells_by_category[1].lines == 2:3
    @test isempty(pr.cells_by_category[1].polygons)
    @test pr.cells_by_category[2].lines == 4:4
    @test pr.cells_by_category[2].polygons == 5:5
    @test marker[pr.cells_by_category[1].vertices] == [1.0]
    @test marker[pr.cells_by_category[1].lines] == [2.0, 3.0]
    @test marker[pr.cells_by_category[2].lines] == [4.0]
    @test marker[pr.cells_by_category[2].polygons] == [5.0]
    @test r["u"] == collect(1.0:7.0)
    close(r)

    # the cells read back can be written again
    fn2 = joinpath(dir, "pd2")
    vtkhdf_open(fn) do r
        c = read_cells(r)
        vtkhdf_grid(fn2, read_points(r), c.vertices, c.lines, c.polygons) do vtk
        end
    end
    vtkhdf_open(fn2) do r2
        @test VTKHDF.ncells(r2) == 5
        @test length(read_cells(r2).lines) == 3
    end

    # empty polydata
    fn3 = joinpath(dir, "empty")
    vtkhdf_grid(VTKPolyData(), fn3) do vtk
    end
    vtkhdf_open(fn3) do r
        @test VTKHDF.npoints(r) == 0
        @test all(isempty, read_cells(r))
    end
end

@testset "read polydata temporal" begin
    dir = mktempdir()
    p1 = Float64[0 1 2; 0 1 2; 0 0 0]
    lines1 = [MeshCell(PolyData.Lines(), [1, 2]), MeshCell(PolyData.Lines(), [2, 3])]
    p2 = Float64[5 6 7 8; 0 0 1 1; 0 0 0 0]
    lines2 = [MeshCell(PolyData.Lines(), [1, 2])]
    polys2 = [MeshCell(PolyData.Polys(), [1, 2, 3, 4])]

    fn = joinpath(dir, "tpd")
    vtk = vtkhdf_grid(fn, p1, lines1; temporal = true)
    write_timestep(vtk, 1.0) do f
        f["u"] = [1.0, 2.0, 3.0]
        f["c", VTKCellData()] = [1.0, 2.0]
    end
    write_timestep(vtk, 2.0; points = p2, cells = (lines2, polys2)) do f
        f["u"] = [4.0, 5.0, 6.0, 7.0]
        f["c", VTKCellData()] = [3.0, 4.0]
    end
    write_timestep(vtk, 3.0) do f  # geometry reused
        f["u"] = [8.0, 9.0, 10.0, 11.0]
        f["c", VTKCellData()] = [5.0, 6.0]
    end
    close(vtk)

    vtkhdf_open(fn) do r
        @test VTKHDF.nsteps(r) == 3
        @test VTKHDF.time_values(r) == [1.0, 2.0, 3.0]
        s1, s2, s3 = (read_timestep(r, i) for i in 1:3)
        @test read_points(s1) == p1
        @test read_points(s2) == p2
        @test read_points(s3) == p2
        @test length(read_cells(s1).lines) == 2
        @test isempty(read_cells(s1).polygons)
        @test length(read_cells(s2).lines) == 1
        @test length(read_cells(s2).polygons) == 1
        @test length(read_cells(s3).polygons) == 1
        @test s1["u"] == [1.0, 2.0, 3.0]
        @test s2["u"] == [4.0, 5.0, 6.0, 7.0]
        @test s3["u"] == [8.0, 9.0, 10.0, 11.0]
        @test s1["c"] == [1.0, 2.0]
        @test s2["c"] == [3.0, 4.0]
        @test s3["c"] == [5.0, 6.0]
        @test VTKHDF.npoints(s2) == 4
        @test VTKHDF.ncells(s3) == 2
        @test VTKHDF.partition_ranges(s2).cells_by_category[1].polygons == 2:2
    end
end
