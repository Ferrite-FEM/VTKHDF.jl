# Reading back ImageData / RectilinearGrid / StructuredGrid files.

@testset "read image data" begin
    dir = mktempdir()
    fn = joinpath(dir, "img")
    T = reshape(collect(1.0:24.0), 4, 3, 2)
    vec = rand(3, 4, 3, 2)
    cd = reshape(collect(1.0:6.0), 3, 2, 1)
    vtkhdf_grid(fn, 0.0:0.5:1.5, 0.0:1.0:2.0, 0.0:2.0:2.0) do vtk
        vtk["T"] = T
        vtk["vec", VTKPointData()] = vec
        vtk["cd", VTKCellData()] = cd
        vtk["e", VTKFieldData()] = [4.0]
    end
    vtkhdf_open(fn) do r
        @test VTKHDF.dataset_type(r) == "ImageData"
        gi = VTKHDF.grid_info(r)
        @test gi.dims == (4, 3, 2)
        @test gi.origin == (0.0, 0.0, 0.0)
        @test gi.spacing == (0.5, 1.0, 2.0)
        @test gi.direction == (1.0, 0, 0, 0, 1, 0, 0, 0, 1)
        @test gi.whole_extent == (0, 3, 0, 2, 0, 1)
        @test VTKHDF.npoints(r) == 24
        @test VTKHDF.ncells(r) == 6
        @test r["T"] == T
        @test r["vec"] == vec
        @test r["cd", VTKCellData()] == cd
        @test r["e"] == [4.0]
        @test keys(r, VTKPointData()) == ["T", "vec"]
    end

    # nonzero whole_extent, degenerate dimension, Int data
    fn = joinpath(dir, "img2")
    vtkhdf_grid(VTKImageData(), fn, (3, 2); whole_extent = (2, 4, 5, 6, 0, 0)) do vtk
        vtk["a"] = reshape(Int32.(1:6), 3, 2, 1)
    end
    vtkhdf_open(fn) do r
        @test VTKHDF.grid_info(r).whole_extent == (2, 4, 5, 6, 0, 0)
        @test VTKHDF.grid_info(r).dims == (3, 2, 1)
        @test r["a"] == reshape(Int32.(1:6), 3, 2, 1)
        @test r["a"] isa Array{Int32, 3}
    end

    # temporal
    fn = joinpath(dir, "timg")
    vtk = vtkhdf_grid(VTKImageData(), fn, (2, 2); temporal = true, compress = true)
    A1 = [1.0 2.0; 3.0 4.0]
    A2 = [5.0 6.0; 7.0 8.0]
    write_timestep(vtk, 0.0) do f
        f["A"] = A1
        f["e", VTKFieldData()] = [1.0, 2.0]
    end
    write_timestep(vtk, 1.0) do f
        f["A"] = A2
        f["e", VTKFieldData()] = [3.0, 4.0]
    end
    close(vtk)
    vtkhdf_open(fn) do r
        @test VTKHDF.nsteps(r) == 2
        @test_throws ErrorException r["A"]
        @test read_timestep(r, 1)["A"] == reshape(A1, 2, 2, 1)
        @test read_timestep(r, 2)["A"] == reshape(A2, 2, 2, 1)
        @test read_timestep(r, 1)["e"] == [1.0, 2.0]
        @test read_timestep(r, 2)["e"] == [3.0, 4.0]
    end
end

@testset "read rectilinear" begin
    dir = mktempdir()
    x = [0.0, 0.5, 2.0]
    y = [0.0, 1.0]
    fn = joinpath(dir, "rect")
    p = rand(3, 2, 1)
    vtkhdf_grid(fn, x, y) do vtk
        vtk["p"] = p
    end
    vtkhdf_open(fn) do r
        @test VTKHDF.dataset_type(r) == "RectilinearGrid"
        @test VTKHDF.grid_info(r) == (dims = (3, 2, 1), whole_extent = (0, 2, 0, 1, 0, 0))
        cx, cy, cz = read_coordinates(r)
        @test cx == x && cy == y && cz == [0.0]
        @test r["p"] == p
        @test VTKHDF.npoints(r) == 6
        @test VTKHDF.ncells(r) == 2
    end

    # temporal with a coordinate change and reuse
    fn = joinpath(dir, "trect")
    vtk = vtkhdf_grid(VTKRectilinearGrid(), fn, x, y; temporal = true)
    x2 = [0.0, 1.0, 4.0]
    write_timestep(f -> f["p"] = 1 .* ones(3, 2, 1), vtk, 0.0)
    write_timestep(f -> f["p"] = 2 .* ones(3, 2, 1), vtk, 1.0; x = x2)
    write_timestep(f -> f["p"] = 3 .* ones(3, 2, 1), vtk, 2.0)
    close(vtk)
    vtkhdf_open(fn) do r
        @test_throws ErrorException read_coordinates(r)
        @test read_coordinates(read_timestep(r, 1))[1] == x
        @test read_coordinates(read_timestep(r, 2))[1] == x2
        @test read_coordinates(read_timestep(r, 3))[1] == x2
        @test read_coordinates(read_timestep(r, 3))[2] == y
        @test read_timestep(r, 3)["p"] == 3 .* ones(3, 2, 1)
    end
end

@testset "read structured" begin
    dir = mktempdir()
    xyz = rand(3, 3, 2, 2)
    fn = joinpath(dir, "sg")
    u = rand(3, 2, 2)
    vtkhdf_grid(fn, xyz) do vtk
        vtk["u"] = u
    end
    vtkhdf_open(fn) do r
        @test VTKHDF.dataset_type(r) == "StructuredGrid"
        @test VTKHDF.grid_info(r).dims == (3, 2, 2)
        @test read_points(r) == xyz
        @test r["u"] == u
        @test VTKHDF.npoints(r) == 12
        @test VTKHDF.ncells(r) == 2
    end

    # temporal: changed then reused points
    fn = joinpath(dir, "tsg")
    xyz2 = xyz .+ 1
    vtk = vtkhdf_grid(fn, xyz; temporal = true)
    write_timestep(f -> f["u"] = 1 .* ones(3, 2, 2), vtk, 0.0)
    write_timestep(f -> f["u"] = 2 .* ones(3, 2, 2), vtk, 1.0; points = xyz2)
    write_timestep(f -> f["u"] = 3 .* ones(3, 2, 2), vtk, 2.0)
    close(vtk)
    vtkhdf_open(fn) do r
        @test_throws ErrorException read_points(r)
        @test read_points(read_timestep(r, 1)) == xyz
        @test read_points(read_timestep(r, 2)) == xyz2
        @test read_points(read_timestep(r, 3)) == xyz2
        @test read_timestep(r, 2)["u"] == 2 .* ones(3, 2, 2)
    end
end
