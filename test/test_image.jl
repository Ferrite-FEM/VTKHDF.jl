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
        # VTK expects (and writes) the arrays without the directions that
        # have a single point: rank 2 for a 2-D image, rank 1 for a 1-D one
        fn = joinpath(dir, "img2d.vtkhdf")
        u = rand(4, 3)
        v = rand(2, 4, 3)
        c = rand(3, 2)
        vtkhdf_grid(fn, 0:1.0:3, 0:1.0:2) do vtk
            vtk["u"] = u                   # nz = 1 implied
            vtk["v"] = v
            vtk["c", VTKCellData()] = c
        end
        h5open(fn) do f
            @test read(f["VTKHDF/PointData/u"]) == u
            @test read(f["VTKHDF/PointData/v"]) == v
            @test read(f["VTKHDF/CellData/c"]) == c
        end
        if HAS_H5DUMP
            @test h5dump_shape(fn, "/VTKHDF/PointData/u") == (3, 4)
            @test h5dump_shape(fn, "/VTKHDF/PointData/v") == (3, 4, 2)
            @test h5dump_shape(fn, "/VTKHDF/CellData/c") == (2, 3)
        end
        if HAS_VTK
            s = only(vtkdump(fn)["steps"])
            @test s["dimensions"] == [4, 3, 1]
            @test s["point_data"]["u"] ≈ vec(u)
            @test s["point_data"]["v"] ≈ nested(reshape(v, 2, :))
            @test s["cell_data"]["c"] ≈ vec(c)
        end
        # the same image from a (nx, ny, 1) array
        fn3 = joinpath(dir, "img2d_3.vtkhdf")
        vtkhdf_grid(fn3, 0:1.0:3, 0:1.0:2) do vtk
            vtk["u"] = reshape(u, 4, 3, 1)
        end
        HAS_H5DUMP && @test h5dump_shape(fn3, "/VTKHDF/PointData/u") == (3, 4)

        fn1 = joinpath(dir, "img1d.vtkhdf")
        u1 = rand(5)
        c1 = rand(4)
        vtkhdf_grid(fn1, 0:0.25:1) do vtk
            vtk["u"] = u1
            vtk["c", VTKCellData()] = c1
        end
        if HAS_H5DUMP
            @test h5dump_shape(fn1, "/VTKHDF/PointData/u") == (5,)
            @test h5dump_shape(fn1, "/VTKHDF/CellData/c") == (4,)
        end
        if HAS_VTK
            s = only(vtkdump(fn1)["steps"])
            @test s["dimensions"] == [5, 1, 1]
            @test s["point_data"]["u"] ≈ u1
            @test s["cell_data"]["c"] ≈ c1
        end

        # temporal 2-D: (t, y, x) on disk
        fnt = joinpath(dir, "img2d_t.vtkhdf")
        vtk = vtkhdf_grid(fnt, 0:1.0:3, 0:1.0:2; temporal = true)
        us = [rand(4, 3) for _ in 1:2]
        for s in 1:2
            write_timestep(frame -> frame["u"] = us[s], vtk, Float64(s))
        end
        close(vtk)
        HAS_H5DUMP && @test h5dump_shape(fnt, "/VTKHDF/PointData/u") == (2, 3, 4)
        if HAS_VTK
            d = vtkdump(fnt)
            @test d["time_steps"] == [1.0, 2.0]
            for s in 1:2
                @test d["steps"][s]["point_data"]["u"] ≈ vec(us[s])
            end
        end
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
