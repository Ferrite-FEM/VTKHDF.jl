# StructuredGrid was added in spec 2.7 and is read by VTK >= 9.7; with an
# older VTK the validation is structural only.
@testset "StructuredGrid" begin
    dir = mktempdir()
    nx, ny, nz = 4, 3, 2
    xyz = Array{Float64}(undef, 3, nx, ny, nz)
    for k in 1:nz, j in 1:ny, i in 1:nx
        xyz[:, i, j, k] .= (i - 1.0, 2j - 2.0, 10k - 10.0)
    end

    @testset "static" begin
        fn = joinpath(dir, "struct.vtkhdf")
        u = rand(nx, ny, nz)
        vtkhdf_grid(fn, xyz) do vtk
            vtk["u"] = u
        end
        h5open(fn) do f
            g = f["VTKHDF"]
            @test attrs(g)["Type"] == "StructuredGrid"
            @test attrs(g)["Dimensions"] == [4, 3, 2]
            @test attrs(g)["WholeExtent"] == [0, 3, 0, 2, 0, 1]
            @test read(g["Points"]) == xyz
            @test read(g["PointData/u"]) == u
        end
        HAS_H5DUMP && @test h5dump_shape(fn, "/VTKHDF/Points") == (2, 3, 4, 3)
        if HAS_VTK_97
            s = only(vtkdump(fn)["steps"])
            @test s["class"] == "vtkStructuredGrid"
            @test s["npoints"] == nx * ny * nz
            @test s["ncells"] == (nx - 1) * (ny - 1) * (nz - 1)
            @test s["points"] ≈ nested(reshape(xyz, 3, :))
            @test s["point_data"]["u"] ≈ vec(u)
        end
    end

    @testset "2-D" begin
        # Points keep all three directions, data arrays drop the degenerate one
        fn = joinpath(dir, "struct2d.vtkhdf")
        u = rand(nx, ny)
        vtkhdf_grid(fn, xyz[:, :, :, 1:1]) do vtk
            vtk["u"] = u
        end
        h5open(fn) do f
            @test attrs(f["VTKHDF"])["Dimensions"] == [4, 3, 1]
        end
        if HAS_H5DUMP
            @test h5dump_shape(fn, "/VTKHDF/Points") == (1, 3, 4, 3)
            @test h5dump_shape(fn, "/VTKHDF/PointData/u") == (3, 4)
        end
        if HAS_VTK_97
            s = only(vtkdump(fn)["steps"])
            @test s["npoints"] == nx * ny
            @test s["ncells"] == (nx - 1) * (ny - 1)
            @test s["point_data"]["u"] ≈ vec(u)
        end
    end

    @testset "three coordinate arrays" begin
        fn = joinpath(dir, "struct2.vtkhdf")
        vtkhdf_grid(
            identity, fn,
            xyz[1, :, :, :], xyz[2, :, :, :], xyz[3, :, :, :]
        )
        h5open(fn) do f
            @test read(f["VTKHDF/Points"]) == xyz
        end
    end

    @testset "temporal with moving points" begin
        fn = joinpath(dir, "struct_t.vtkhdf")
        vtk = vtkhdf_grid(fn, xyz; temporal = true)
        us = [rand(nx, ny, nz) for _ in 1:2]
        write_timestep(vtk, 0.0) do frame
            frame["u"] = us[1]
        end
        write_timestep(vtk, 1.0; points = xyz .+ 1) do frame
            frame["u"] = us[2]
        end
        close(vtk)
        h5open(fn) do f
            g = f["VTKHDF"]
            @test read(g["Steps/PointOffsets"]) == [0, 1]
        end
        # on-disk 5-D points: (t, nz, ny, nx, 3)
        HAS_H5DUMP && @test h5dump_shape(fn, "/VTKHDF/Points") == (2, 2, 3, 4, 3)
        if HAS_VTK_97
            d = vtkdump(fn)
            @test d["time_steps"] == [0.0, 1.0]
            @test d["steps"][1]["points"] ≈ nested(reshape(xyz, 3, :))
            @test d["steps"][2]["points"] ≈ nested(reshape(xyz .+ 1, 3, :))
            @test d["steps"][2]["point_data"]["u"] ≈ vec(us[2])
        end
    end
end
