# Reading back OverlappingAMR and HyperTreeGrid files (static only).

@testset "read amr" begin
    dir = mktempdir()
    fn = joinpath(dir, "amr")
    rho0 = collect(1.0:64.0)
    vel1 = rand(3, 16)
    vtkhdf_amr(fn; origin = (1.0, 2.0, 3.0)) do amr
        l0 = add_level(amr; spacing = (1.0, 1.0, 1.0))
        add_box(l0, (0, 3, 0, 3, 0, 3); celldata = ("rho" => rho0,))
        l1 = add_level(amr; spacing = (0.5, 0.5, 0.5))
        add_box(l1, (0, 1, 0, 1, 0, 1); celldata = ("rho" => collect(1.0:8.0),))
        add_box(l1, (2, 3, 2, 3, 2, 3); celldata = ("rho" => collect(9.0:16.0),))
        l1["vel", VTKCellData()] = vel1
    end

    vtkhdf_open(fn) do r
        @test VTKHDF.dataset_type(r) == "OverlappingAMR"
        @test VTKHDF.nlevels(r) == 2
        @test VTKHDF.grid_info(r) == (origin = (1.0, 2.0, 3.0), grid_description = "XYZ")
        @test_throws ErrorException r["rho"]  # data is per level
        @test_throws ErrorException VTKHDF.amr_level(r, 2)
        @test_throws ErrorException VTKHDF.amr_level(r, -1)

        l0 = VTKHDF.amr_level(r, 0)
        li = VTKHDF.level_info(l0)
        @test li.spacing == (1.0, 1.0, 1.0)
        @test li.boxes == [(0, 3, 0, 3, 0, 3)]
        @test l0["rho"] == rho0
        @test VTKHDF.ncells(l0) == 64
        @test VTKHDF.npoints(l0) == 125

        l1 = VTKHDF.amr_level(r, 1)
        @test VTKHDF.level_info(l1).boxes == [(0, 1, 0, 1, 0, 1), (2, 3, 2, 3, 2, 3)]
        @test l1["rho", VTKCellData()] == collect(1.0:16.0)
        @test l1["vel"] == vel1
        @test keys(l1, VTKCellData()) == ["rho", "vel"]
        @test haskey(l1, "rho", VTKCellData())
        @test !haskey(l1, "rho", VTKPointData())
        @test_throws ErrorException l1["nope"]
        pr = VTKHDF.partition_ranges(l1)
        @test pr.cells == [1:8, 9:16]
        @test pr.points == [1:27, 28:54]
    end
end

@testset "read htg" begin
    dir = mktempdir()
    fn = joinpath(dir, "htg")
    vtkhdf_htg(fn; dimensions = (3, 2, 1), branch_factor = 2) do htg
        add_piece(
            htg;
            descriptors = [true],
            depth_per_tree = [2, 1],
            tree_ids = [0, 1],
            number_of_cells_per_tree_depth = [1, 4, 1],
            xcoordinates = [0.0, 1.0, 2.0], ycoordinates = [0.0, 1.0], zcoordinates = [0.0],
            mask = [false, false, true, false, false, false],
            celldata = ("q" => collect(1.0:6.0),),
        )
        add_piece(
            htg;
            descriptors = Bool[],
            depth_per_tree = [1],
            tree_ids = [1],
            number_of_cells_per_tree_depth = [1],
            xcoordinates = [5.0, 6.0, 7.0], ycoordinates = [0.0, 1.0], zcoordinates = [0.0],
            mask = [true],
            celldata = ("q" => [7.0],),
        )
    end

    vtkhdf_open(fn) do r
        @test VTKHDF.dataset_type(r) == "HyperTreeGrid"
        @test VTKHDF.npieces(r) == 2
        @test VTKHDF.ncells(r) == 7
        gi = VTKHDF.grid_info(r)
        @test gi.dimensions == (3, 2, 1)
        @test gi.branch_factor == 2
        @test gi.transposed_root_indexing == false

        p1 = VTKHDF.htg_piece(r, 1)
        @test p1.descriptors == [true]
        @test p1.tree_ids == [0, 1]
        @test p1.depth_per_tree == [2, 1]
        @test p1.number_of_cells_per_tree_depth == [1, 4, 1]
        @test p1.xcoordinates == [0.0, 1.0, 2.0]
        @test p1.ycoordinates == [0.0, 1.0]
        @test p1.zcoordinates == [0.0]
        @test p1.mask == [false, false, true, false, false, false]
        @test p1.celldata["q"] == collect(1.0:6.0)

        p2 = VTKHDF.htg_piece(r, 2)
        @test p2.descriptors == Bool[]
        @test p2.tree_ids == [1]
        @test p2.xcoordinates == [5.0, 6.0, 7.0]
        @test p2.mask == [true]
        @test p2.celldata["q"] == [7.0]

        @test r["q"] == collect(1.0:7.0)   # whole-file cell data
        @test_throws ErrorException VTKHDF.htg_piece(r, 3)
    end

    # no mask, no cell data
    fn = joinpath(dir, "htg2")
    vtkhdf_htg(fn; dimensions = (2, 1, 1)) do htg
        add_piece(
            htg;
            descriptors = Bool[],
            depth_per_tree = [1],
            tree_ids = [0],
            number_of_cells_per_tree_depth = [1],
            xcoordinates = [0.0, 1.0], ycoordinates = [0.0], zcoordinates = [0.0],
        )
    end
    vtkhdf_open(fn) do r
        p = VTKHDF.htg_piece(r, 1)
        @test p.mask === nothing
        @test isempty(p.celldata)
    end
end
