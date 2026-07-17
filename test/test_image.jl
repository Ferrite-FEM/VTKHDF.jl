@testset "ImageData" begin
    dir = mktempdir()

    @testset "static" begin
        fn = joinpath(dir, "img.vtkhdf")
        u = rand(5, 4, 2)
        v = rand(3, 5, 4, 2)
        c = rand(4, 3, 1)
        vtkhdf_grid(fn, 0:0.5:2, 1:1.0:4, 0:2.0:2) do vtk
            vtk["u"] = u
            vtk["v"] = v
            vtk["c"] = c
        end
        h5open(fn) do f
            g = f["VTKHDF"]
            @test attrs(g)["Type"] == "ImageData"
            @test attrs(g)["WholeExtent"] == [0, 4, 0, 3, 0, 1]
            @test attrs(g)["Origin"] == [0.0, 1.0, 0.0]
            @test attrs(g)["Spacing"] == [0.5, 1.0, 2.0]
            @test attrs(g)["Direction"] == [1, 0, 0, 0, 1, 0, 0, 0, 1]
            @test read(g["PointData/u"]) == u
            @test read(g["PointData/v"]) == v
            @test read(g["CellData/c"]) == c
        end
        if HAS_H5DUMP
            # on-disk (z, y, x) / (z, y, x, ncomp)
            @test h5dump_shape(fn, "/VTKHDF/PointData/u") == (2, 4, 5)
            @test h5dump_shape(fn, "/VTKHDF/PointData/v") == (2, 4, 5, 3)
            @test h5dump_shape(fn, "/VTKHDF/CellData/c") == (1, 3, 4)
        end
        if HAS_VTK
            d = vtkdump(fn)
            s = only(d["steps"])
            @test s["class"] == "vtkImageData"
            @test s["dimensions"] == [5, 4, 2]
            @test s["origin"] == [0, 1, 0]
            @test s["spacing"] == [0.5, 1, 2]
            # VTK's flat point order is x fastest, then y, then z — identical
            # to the Julia memory order of u
            @test s["point_data"]["u"] ≈ vec(u)
            @test s["point_data"]["v"] ≈ nested(reshape(v, 3, :))
            @test s["cell_data"]["c"] ≈ vec(c)
        end
    end

    @testset "direction and extent" begin
        fn = joinpath(dir, "img_dir.vtkhdf")
        dir9 = (0.0, 1, 0, -1, 0, 0, 0, 0, 1)  # 90° rotation
        vtkhdf_grid(
            identity, VTKImageData(), fn, (3, 3, 1);
            origin = (1.0, 2.0, 3.0), spacing = (0.1, 0.2, 0.3), direction = dir9,
            whole_extent = (10, 12, 20, 22, 5, 5)
        )
        h5open(fn) do f
            @test attrs(f["VTKHDF"])["WholeExtent"] == [10, 12, 20, 22, 5, 5]
        end
        if HAS_VTK
            s = only(vtkdump(fn)["steps"])
            @test s["direction"] == collect(dir9)
            @test s["origin"] == [1, 2, 3]
        end
    end

    @testset "2-D and 1-D" begin
        fn = joinpath(dir, "img2d.vtkhdf")
        vtkhdf_grid(fn, 0:1.0:3, 0:1.0:2) do vtk
            vtk["u"] = rand(4, 3)          # nz = 1 implied
            vtk["c", VTKCellData()] = rand(3, 2)
        end
        HAS_VTK && @test only(vtkdump(fn)["steps"])["dimensions"] == [4, 3, 1]
        fn1 = joinpath(dir, "img1d.vtkhdf")
        vtkhdf_grid(fn1, 0:0.25:1) do vtk
            vtk["u"] = rand(5)
            vtk["c", VTKCellData()] = rand(4)
        end
        HAS_VTK && @test only(vtkdump(fn1)["steps"])["dimensions"] == [5, 1, 1]
    end

    @testset "temporal" begin
        fn = joinpath(dir, "img_t.vtkhdf")
        us = [rand(4, 3, 2) for _ in 1:3]
        vs = [rand(3, 4, 3, 2) for _ in 1:3]
        vtk = vtkhdf_grid(fn, 0:1.0:3, 0:1.0:2, 0:1.0:1; temporal = true)
        for s in 1:3
            write_timestep(vtk, s / 4) do frame
                frame["u"] = us[s]
                frame["v"] = vs[s]
            end
        end
        close(vtk)
        if HAS_H5DUMP
            # time-prepended on disk: (t, z, y, x[, ncomp])
            @test h5dump_shape(fn, "/VTKHDF/PointData/u") == (3, 2, 3, 4)
            @test h5dump_shape(fn, "/VTKHDF/PointData/v") == (3, 2, 3, 4, 3)
        end
        if HAS_VTK
            d = vtkdump(fn)
            @test d["time_steps"] == [0.25, 0.5, 0.75]
            for s in 1:3
                @test d["steps"][s]["point_data"]["u"] ≈ vec(us[s])
                @test d["steps"][s]["point_data"]["v"] ≈ nested(reshape(vs[s], 3, :))
            end
        end
    end
end
