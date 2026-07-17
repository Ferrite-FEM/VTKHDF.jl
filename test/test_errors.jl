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

    # NUL bytes and unknown attribute kinds are rejected (UTF-8 itself is fine)
    vtk = vtkhdf_grid(joinpath(dir, "e8"), cube, hex)
    @test_throws ArgumentError vtk["a\0b"] = rand(8)
    @test_throws ArgumentError vtk["u", VTKPointData(), attribute = :NotAThing] = rand(8)
    vtk["ünicode"] = rand(8)
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

    # a failed timestep marks the file as failed instead of wedging it
    vtk = vtkhdf_grid(joinpath(dir, "efail"), cube, hex; temporal = true)
    write_timestep(f -> (f["u"] = rand(8)), vtk, 0.0)
    @test_throws ErrorException write_timestep(vtk, 1.0) do f
        f["u"] = rand(8)
        error("user error")
    end
    @test_throws ErrorException write_timestep(f -> (f["u"] = rand(8)), vtk, 2.0)
    @test_logs (:warn, r"failed write") close(vtk)   # close still works

    # appends cannot silently change the element type
    vtk = vtkhdf_grid(VTKUnstructuredGrid(), joinpath(dir, "eeltype"))
    add_partition(vtk, cube, hex; pointdata = ("u" => rand(Float64, 8),))
    @test_throws ErrorException add_partition(vtk, Float32.(cube), hex)
    @test_throws ErrorException vtk["u", VTKPointData()] = rand(Float32, 8)
    add_partition(vtk, cube .+ 2, hex; pointdata = ("u" => rand(Float64, 8),))
    close(vtk)

    # out-of-range connectivity is rejected before anything is written
    @test_throws ArgumentError vtkhdf_grid(
        joinpath(dir, "econn"), cube,
        [MeshCell(VTKCellTypes.VTK_HEXAHEDRON, 2:9)]
    )
    @test_throws ArgumentError vtkhdf_grid(
        joinpath(dir, "econn0"), cube,
        [MeshCell(VTKCellTypes.VTK_HEXAHEDRON, 0:7)]
    )

    # UTF-8 names are allowed (only '/' and '.' are forbidden by the format)
    fnu = joinpath(dir, "utf8.vtkhdf")
    vtkhdf_grid(fnu, cube, hex) do vtk
        vtk["ρ"] = rand(8)
    end
    h5open(fnu) do f
        @test haskey(f, "VTKHDF/PointData/ρ")
    end
    if HAS_VTK
        @test "ρ" in keys(only(vtkdump(fnu)["steps"])["point_data"])
    end

    # closing a deferred-geometry file without partitions writes an empty layout
    fne = joinpath(dir, "empty.vtkhdf")
    close(vtkhdf_grid(VTKUnstructuredGrid(), fne))
    h5open(fne) do f
        @test read(f["VTKHDF/NumberOfPoints"]) == [0]
        @test read(f["VTKHDF/Offsets"]) == [0]
    end
    if HAS_VTK
        s = only(vtkdump(fne)["steps"])
        @test s["npoints"] == 0
        @test s["ncells"] <= 0  # VTK reports -1 cells for a fully empty grid
    end

    # temporal file closed without steps still materializes Steps
    fnt = joinpath(dir, "nosteps.vtkhdf")
    vtk = vtkhdf_grid(fnt, cube, hex; temporal = true)
    @test_logs (:warn, r"without any time steps") close(vtk)
    h5open(fnt) do f
        @test attrs(f["VTKHDF/Steps"])["NSteps"] == 0
        @test length(f["VTKHDF/Steps/Values"]) == 0
    end

    # mismatched whole_extent
    @test_throws ArgumentError vtkhdf_grid(
        VTKImageData(), joinpath(dir, "eext"),
        (4, 3, 2); whole_extent = (0, 4, 0, 2, 0, 1)
    )

    # a 1×N matrix is written as a scalar (1-D) array
    fn1 = joinpath(dir, "mat1.vtkhdf")
    vtkhdf_grid(fn1, cube, hex) do vtk
        vtk["s"] = reshape(Float64.(1:8), 1, 8)
    end
    h5open(fn1) do f
        @test read(f["VTKHDF/PointData/s"]) == 1:8
        @test ndims(f["VTKHDF/PointData/s"]) == 1
    end

    # vectors of static-vector-like eltypes
    fn = joinpath(dir, "svec.vtkhdf")
    vtkhdf_grid(fn, cube, hex) do vtk
        vtk["nt"] = [(1.0, 2.0, 3.0) for _ in 1:8]
    end
    h5open(fn) do f
        @test size(read(f["VTKHDF/PointData/nt"])) == (3, 8)
    end
end
