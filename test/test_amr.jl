@testset "OverlappingAMR" begin
    dir = mktempdir()
    fn = joinpath(dir, "amr.vtkhdf")
    g0 = rand(125)
    p0 = rand(216)
    g1a, g1b = rand(240), rand(240)
    vtkhdf_amr(fn; origin = (-2.0, -2.0, 0.0)) do amr
        l0 = add_level(amr; spacing = (0.5, 0.5, 0.5))
        add_box(l0, (0, 4, 0, 4, 0, 4); celldata = ("g" => g0,), pointdata = ("p" => p0,))
        l1 = add_level(amr; spacing = (0.25, 0.25, 0.25))
        add_box(l1, (0, 3, 0, 5, 0, 9); celldata = ("g" => g1a,))
        add_box(l1, (6, 9, 4, 9, 0, 9); celldata = ("g" => g1b,))
    end
    h5open(fn) do f
        g = f["VTKHDF"]
        @test attrs(g)["Type"] == "OverlappingAMR"
        @test attrs(g)["Origin"] == [-2.0, -2.0, 0.0]
        @test attrs(g)["GridDescription"] == "XYZ"
        @test attrs(g["Level0"])["Spacing"] == [0.5, 0.5, 0.5]
        @test read(g["Level0/AMRBox"]) == reshape([0, 4, 0, 4, 0, 4], 6, 1)
        @test size(read(g["Level1/AMRBox"])) == (6, 2)
        @test read(g["Level1/CellData/g"]) == [g1a; g1b]
    end
    HAS_H5DUMP && @test h5dump_shape(fn, "/VTKHDF/Level1/AMRBox") == (2, 6)
    if HAS_VTK
        d = vtkdump(fn)
        s = only(d["steps"])
        @test s["class"] == "vtkOverlappingAMR"
        @test s["bounds"][[1, 3, 5]] == [-2, -2, 0]  # min corner == origin here
        @test length(s["levels"]) == 2
        @test s["levels"][1]["spacing"] == [0.5, 0.5, 0.5]
        @test length(s["levels"][2]["boxes"]) == 2
        @test s["levels"][1]["boxes"][1]["cell_data"]["g"] ≈ g0
        @test s["levels"][2]["boxes"][2]["cell_data"]["g"] ≈ g1b
    end

    # size validation
    amr = vtkhdf_amr(joinpath(dir, "amr_bad.vtkhdf"); origin = (0, 0, 0))
    lvl = add_level(amr; spacing = (1.0, 1.0, 1.0))
    @test_throws Exception add_box(lvl, (0, 1, 0, 1, 0, 1); celldata = ("g" => rand(7),))
end
