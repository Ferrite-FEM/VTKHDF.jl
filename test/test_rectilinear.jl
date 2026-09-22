# RectilinearGrid was added in spec 2.7 and is read by VTK >= 9.7; with an
# older VTK the validation is structural only.
@testset "RectilinearGrid" begin
    dir = mktempdir()
    x, y, z = [0.0, 1.0, 2.5, 4.0], [0.0, 2.0, 3.0], [0.0, 10.0]

    @testset "static" begin
        fn = joinpath(dir, "rect.vtkhdf")
        u = rand(4, 3, 2)
        c = rand(3, 2, 1)
        vtkhdf_grid(fn, x, y, z) do vtk
            vtk["u"] = u
            vtk["c", VTKCellData()] = c
        end
        h5open(fn) do f
            g = f["VTKHDF"]
            @test attrs(g)["Type"] == "RectilinearGrid"
            @test attrs(g)["Version"] == [2, 7]
            @test attrs(g)["Dimensions"] == [4, 3, 2]
            @test attrs(g)["WholeExtent"] == [0, 3, 0, 2, 0, 1]
            @test read(g["XCoordinates"]) == x
            @test read(g["YCoordinates"]) == y
            @test read(g["ZCoordinates"]) == z
            @test read(g["PointData/u"]) == u
        end
        if HAS_H5DUMP
            @test h5dump_shape(fn, "/VTKHDF/PointData/u") == (2, 3, 4)
            @test h5dump_shape(fn, "/VTKHDF/CellData/c") == (1, 2, 3)
        end
        if HAS_VTK_97
            s = only(vtkdump(fn)["steps"])
            @test s["class"] == "vtkRectilinearGrid"
            @test s["dimensions"] == [4, 3, 2]
            @test s["x"] == x && s["y"] == y && s["z"] == z
            @test s["point_data"]["u"] ≈ vec(u)
            @test s["cell_data"]["c"] ≈ vec(c)
        end
    end

    @testset "2-D" begin
        fn = joinpath(dir, "rect2d.vtkhdf")
        u = rand(4, 3)
        c = rand(3, 2)
        vtkhdf_grid(fn, x, y) do vtk
            vtk["u"] = u
            vtk["c", VTKCellData()] = c
        end
        h5open(fn) do f
            @test attrs(f["VTKHDF"])["Dimensions"] == [4, 3, 1]
        end
        if HAS_H5DUMP
            @test h5dump_shape(fn, "/VTKHDF/PointData/u") == (3, 4)
            @test h5dump_shape(fn, "/VTKHDF/CellData/c") == (2, 3)
        end
        if HAS_VTK_97
            s = only(vtkdump(fn)["steps"])
            @test s["dimensions"] == [4, 3, 1]
            @test s["point_data"]["u"] ≈ vec(u)
            @test s["cell_data"]["c"] ≈ vec(c)
        end
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
        us = [rand(4, 3, 2) for _ in 1:3]
        write_timestep(vtk, 0.0) do frame
            frame["u"] = us[1]
        end
        write_timestep(vtk, 1.0; x = x .+ 1) do frame
            frame["u"] = us[2]
        end
        write_timestep(vtk, 2.0) do frame
            frame["u"] = us[3]
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
        if HAS_VTK_97
            d = vtkdump(fn)
            @test d["time_steps"] == [0.0, 1.0, 2.0]
            @test d["steps"][1]["x"] == x
            @test d["steps"][2]["x"] == x .+ 1
            @test d["steps"][3]["x"] == x .+ 1
            @test d["steps"][3]["y"] == y
            for s in 1:3
                @test d["steps"][s]["point_data"]["u"] ≈ vec(us[s])
            end
        end
    end
end
