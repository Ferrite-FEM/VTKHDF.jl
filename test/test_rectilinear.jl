# RectilinearGrid was added in spec 2.7; VTK's reader may not support it yet
# (VTK 9.6 predates it), so validation is structural.
@testset "RectilinearGrid" begin
    dir = mktempdir()
    x, y, z = [0.0, 1.0, 2.5, 4.0], [0.0, 2.0, 3.0], [0.0, 10.0]

    @testset "static" begin
        fn = joinpath(dir, "rect.vtkhdf")
        u = rand(4, 3, 2)
        vtkhdf_grid(fn, x, y, z) do vtk
            vtk["u"] = u
            vtk["c", VTKCellData()] = rand(3, 2, 1)
        end
        h5open(fn) do f
            g = f["VTKHDF"]
            @test attrs(g)["Type"] == "RectilinearGrid"
            @test attrs(g)["Version"] == [2, 7]
            @test attrs(g)["WholeExtent"] == [0, 3, 0, 2, 0, 1]
            @test read(g["XCoordinates"]) == x
            @test read(g["YCoordinates"]) == y
            @test read(g["ZCoordinates"]) == z
            @test read(g["PointData/u"]) == u
        end
        HAS_H5DUMP && @test h5dump_shape(fn, "/VTKHDF/PointData/u") == (2, 3, 4)
    end

    @testset "mixed ranges and vectors dispatch to rectilinear" begin
        fn = joinpath(dir, "rect2.vtkhdf")
        vtkhdf_grid(identity, fn, 0:1.0:3, [0.0, 2.0, 3.0])
        h5open(fn) do f
            @test attrs(f["VTKHDF"])["Type"] == "RectilinearGrid"
            @test read(f["VTKHDF/ZCoordinates"]) == [0.0]
        end
    end

    @testset "temporal with changing coordinates" begin
        fn = joinpath(dir, "rect_t.vtkhdf")
        vtk = vtkhdf_grid(VTKRectilinearGrid(), fn, x, y, z; temporal = true)
        write_timestep(vtk, 0.0) do frame
            frame["u"] = rand(4, 3, 2)
        end
        write_timestep(vtk, 1.0; x = x .+ 1) do frame
            frame["u"] = rand(4, 3, 2)
        end
        write_timestep(vtk, 2.0) do frame
            frame["u"] = rand(4, 3, 2)
        end
        close(vtk)
        h5open(fn) do f
            g = f["VTKHDF"]
            @test read(g["XCoordinates"]) == [x; x .+ 1]
            @test read(g["YCoordinates"]) == y
            @test read(g["Steps/XCoordinatesOffsets"]) == [0, 4, 4]
            @test read(g["Steps/YCoordinatesOffsets"]) == [0, 0, 0]
            @test read(g["Steps/ZCoordinatesOffsets"]) == [0, 0, 0]
        end
        HAS_H5DUMP && @test h5dump_shape(fn, "/VTKHDF/PointData/u") == (3, 2, 3, 4)
    end
end
