@testset "errors and edge cases" begin
    dir = mktempdir()
    cube = Float64[0 1 1 0 0 1 1 0; 0 0 1 1 0 0 1 1; 0 0 0 0 1 1 1 1]
    hex = [MeshCell(VTKCellTypes.VTK_HEXAHEDRON, 1:8)]

    # invalid names
    vtk = vtkhdf_grid(joinpath(dir, "e1"), cube, hex)
    @test_throws ArgumentError vtk["bad/name"] = rand(8)
    @test_throws ArgumentError vtk["bad.name"] = rand(8)
    @test_throws ArgumentError vtk[""] = rand(8)
    vtk["ok"] = rand(8)
    close(vtk)

    # ambiguous / mismatched lengths (4 points and 4 vertex cells)
    vertcells = [MeshCell(VTKCellTypes.VTK_VERTEX, [i]) for i in 1:4]
    vtk = vtkhdf_grid(joinpath(dir, "e2"), cube[:, 1:4], vertcells)
    @test_throws ErrorException vtk["amb"] = rand(4)  # npoints == ncells
    vtk["amb", VTKPointData()] = rand(4)
    vtk["ambc", VTKCellData()] = rand(4)
    close(vtk)
    vtk = vtkhdf_grid(joinpath(dir, "e3"), cube, hex)
    @test_throws ErrorException vtk["wrong"] = rand(7)
    vtk["good"] = rand(8)
    close(vtk)

    # static file rejects write_timestep; temporal rejects direct writes
    vtk = vtkhdf_grid(joinpath(dir, "e4"), cube, hex)
    @test_throws ErrorException write_timestep(identity, vtk, 0.0)
    vtk["u"] = rand(8)
    close(vtk)
    vtk = vtkhdf_grid(joinpath(dir, "e5"), cube, hex; temporal = true)
    @test_throws ErrorException vtk["u"] = rand(8)
    write_timestep(f -> (f["u"] = rand(8)), vtk, 0.0)
    close(vtk)

    # temporal schema is frozen by the first step
    vtk = vtkhdf_grid(joinpath(dir, "e6"), cube, hex; temporal = true)
    write_timestep(f -> (f["u"] = rand(8)), vtk, 0.0)
    @test_throws ErrorException write_timestep(f -> (f["v"] = rand(8)), vtk, 1.0)
    @test_throws ErrorException write_timestep(identity, vtk, 2.0)  # missing array

    # incomplete static data
    vtk = vtkhdf_grid(joinpath(dir, "e7"), cube, hex)
    vtk["u", VTKPointData()] = rand(8)
    vtk["u2", VTKPointData()] = rand(4)  # explicit location skips length check...
    @test_throws ErrorException close(vtk) # ...but close validates totals

    # non-ascii and unknown attribute kinds
    vtk = vtkhdf_grid(joinpath(dir, "e8"), cube, hex)
    @test_throws ArgumentError vtk["ünicode"] = rand(8)
    @test_throws ArgumentError vtk["u", VTKPointData(), attribute = :NotAThing] = rand(8)
    close(vtk)

    # do-block closes on exception
    local leaked
    @test_throws ErrorException vtkhdf_grid(joinpath(dir, "e9"), cube, hex) do vtk
        leaked = vtk
        error("boom")
    end
    @test !leaked.isopen

    # closing twice is fine; writing after close is not
    vtk = vtkhdf_grid(joinpath(dir, "e10"), cube, hex)
    close(vtk)
    close(vtk)
    @test_throws ErrorException vtk["u"] = rand(8)

    # filenames get the .vtkhdf extension only when no extension is present
    vtk = vtkhdf_grid(joinpath(dir, "name"), cube, hex)
    close(vtk)
    @test isfile(joinpath(dir, "name.vtkhdf"))
    vtk = vtkhdf_grid(joinpath(dir, "name.hdf"), cube, hex)
    close(vtk)
    @test isfile(joinpath(dir, "name.hdf"))

    # vectors of static-vector-like eltypes
    fn = joinpath(dir, "svec.vtkhdf")
    vtkhdf_grid(fn, cube, hex) do vtk
        vtk["nt"] = [(1.0, 2.0, 3.0) for _ in 1:8]
    end
    h5open(fn) do f
        @test size(read(f["VTKHDF/PointData/nt"])) == (3, 8)
    end
end
