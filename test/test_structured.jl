# StructuredGrid was added in spec 2.7; VTK's reader may not support it yet
# (VTK 9.6 predates it), so validation is structural.
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
            @test attrs(g)["WholeExtent"] == [0, 3, 0, 2, 0, 1]
            @test read(g["Points"]) == xyz
            @test read(g["PointData/u"]) == u
        end
        HAS_H5DUMP && @test h5dump_shape(fn, "/VTKHDF/Points") == (2, 3, 4, 3)
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
        write_timestep(vtk, 0.0) do frame
            frame["u"] = rand(nx, ny, nz)
        end
        write_timestep(vtk, 1.0; points = xyz .+ 1) do frame
            frame["u"] = rand(nx, ny, nz)
        end
        close(vtk)
        h5open(fn) do f
            g = f["VTKHDF"]
            @test read(g["Steps/PointOffsets"]) == [0, 1]
        end
        # on-disk 5-D points: (t, nz, ny, nx, 3)
        HAS_H5DUMP && @test h5dump_shape(fn, "/VTKHDF/Points") == (2, 2, 3, 4, 3)
    end
end
