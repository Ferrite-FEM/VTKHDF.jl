@testset "HyperTreeGrid" begin
    dir = mktempdir()
    fn = joinpath(dir, "htg.vtkhdf")
    # 2×2 grid of trees (dimensions (3,3,1)), branch factor 2; tree 0 is
    # refined once (5 cells), the others are single cells => 8 cells total.
    lvl = Float64[0, 1, 1, 1, 1, 0, 0, 0]
    vtkhdf_htg(fn; dimensions = (3, 3, 1), branch_factor = 2) do htg
        add_piece(
            htg;
            descriptors = [true, false, false, false, false],
            depth_per_tree = [2, 1, 1, 1],
            tree_ids = [0, 1, 2, 3],
            number_of_cells_per_tree_depth = [1, 4, 1, 1, 1],
            xcoordinates = [0.0, 1.0, 2.0],
            ycoordinates = [0.0, 1.0, 2.0],
            zcoordinates = [0.0],
            celldata = ("lvl" => lvl,)
        )
    end
    h5open(fn) do f
        g = f["VTKHDF"]
        @test attrs(g)["Type"] == "HyperTreeGrid"
        @test attrs(g)["Version"] == [2, 4]
        @test attrs(g)["BranchFactor"] == 2
        @test attrs(g)["Dimensions"] == [3, 3, 1]
        @test read(g["Descriptors"]) == UInt8[0x80]  # MSB-first bit packing
        @test read(g["DescriptorsSize"]) == [5]
        @test read(g["NumberOfCells"]) == [8]
        @test read(g["NumberOfDepths"]) == [5]
        @test read(g["NumberOfTrees"]) == [4]
    end
    if HAS_VTK
        d = vtkdump(fn)
        s = only(d["steps"])
        @test s["class"] == "vtkHyperTreeGrid"
        @test s["ncells"] == 8
        @test s["cell_data"]["lvl"] == lvl
        @test s["x"] == [0, 1, 2]
    end

    # a mask; masked cell must disappear from the reader's perspective
    fn2 = joinpath(dir, "htg_mask.vtkhdf")
    vtkhdf_htg(fn2; dimensions = (2, 2, 1), branch_factor = 2) do htg
        add_piece(
            htg;
            descriptors = [false],
            depth_per_tree = [1],
            tree_ids = [0],
            number_of_cells_per_tree_depth = [1],
            xcoordinates = [0.0, 1.0], ycoordinates = [0.0, 1.0], zcoordinates = [0.0],
            mask = [true],
            celldata = ("v" => [1.0],)
        )
    end
    h5open(fn2) do f
        @test read(f["VTKHDF/Mask"]) == UInt8[0x80]
    end

    # point data is rejected
    htg = vtkhdf_htg(joinpath(dir, "htg_bad.vtkhdf"); dimensions = (2, 2, 1))
    @test_throws ArgumentError htg["p", VTKPointData()] = rand(3)
    close(htg)
end
