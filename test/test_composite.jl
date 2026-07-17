@testset "Composite" begin
    dir = mktempdir()
    cube = Float64[0 1 1 0 0 1 1 0; 0 0 1 1 0 0 1 1; 0 0 0 0 1 1 1 1]
    hex = [MeshCell(VTKCellTypes.VTK_HEXAHEDRON, 1:8)]
    sq = Float64[0 1 1 0; 0 0 1 1; 0 0 0 0]
    polys = [MeshCell(PolyData.Polys(), [1, 2, 3, 4])]

    @testset "PartitionedDataSetCollection" begin
        fn = joinpath(dir, "pdc.vtkhdf")
        u = rand(8)
        vtkhdf_collection(fn) do col
            b0 = vtkhdf_grid(col, "Mesh", cube, hex)
            b0["u"] = u
            b1 = vtkhdf_grid(col, "Surf", sq, polys)
            b1["w"] = rand(4)
            img = vtkhdf_grid(col, "Img", VTKImageData(), (2, 2, 2))
            img["s"] = rand(2, 2, 2)
            empty = add_empty_block(col, "Empty")
            add_block_ref(col, empty)
            solids = add_node(col, "solids")
            add_block_ref(solids, b0)
            inner = add_node(solids, "inner")
            add_block_ref(inner, b1)
            add_block_ref(col, img)
        end
        h5open(fn) do f
            g = f["VTKHDF"]
            @test attrs(g)["Type"] == "PartitionedDataSetCollection"
            @test attrs(g["Mesh"])["Index"] == 0
            @test attrs(g["Surf"])["Index"] == 1
            @test attrs(g["Img"])["Index"] == 2
            @test attrs(g["Mesh"])["Type"] == "UnstructuredGrid"
            @test !haskey(attrs(g["Empty"]), "Type")
            @test attrs(g["Empty"])["Index"] == 3  # PDC blocks all carry an Index
            @test haskey(g, "Assembly/Empty")
            @test haskey(g, "Assembly/solids/Mesh")
            @test haskey(g, "Assembly/solids/inner/Surf")
        end
        if HAS_VTK
            d = vtkdump(fn)
            s = only(d["steps"])
            @test s["class"] == "vtkPartitionedDataSetCollection"
            classes = [b["class"] for b in s["blocks"]]
            @test "vtkUnstructuredGrid" in classes && "vtkPolyData" in classes && "vtkImageData" in classes
            @test s["blocks"][1]["point_data"]["u"] ≈ u
            names = [c["name"] for c in s["assembly"]["children"]]
            @test "solids" in names
        end
    end

    @testset "MultiBlockDataSet" begin
        fn = joinpath(dir, "mb.vtkhdf")
        vtkhdf_multiblock(fn) do col
            b0 = vtkhdf_grid(col, "A", cube, hex)
            b0["u"] = rand(8)
            b1 = vtkhdf_grid(col, "B", sq, polys)
            add_block_ref(add_node(col, "g1"), b0)
            add_block_ref(add_node(col, "g2"), b1)
        end
        h5open(fn) do f
            g = f["VTKHDF"]
            @test attrs(g)["Type"] == "MultiBlockDataSet"
            @test !haskey(attrs(g["A"]), "Index")  # MB blocks have no Index
        end
        if HAS_VTK
            s = only(vtkdump(fn)["steps"])
            @test s["class"] == "vtkMultiBlockDataSet"
        end
    end

    @testset "temporal blocks" begin
        fn = joinpath(dir, "tcomp.vtkhdf")
        vtkhdf_collection(fn) do col
            bA = vtkhdf_grid(col, "A", cube, hex; temporal = true)
            bB = vtkhdf_grid(col, "B", sq, polys; temporal = true)
            for s in 1:3
                write_timestep(bA, Float64(s)) do f
                    f["u"] = fill(Float64(s), 8)
                end
                write_timestep(bB, Float64(s)) do f
                    f["w"] = fill(Float64(2s), 4)
                end
            end
            add_block_ref(add_node(col, "all"), bA)
            add_block_ref(add_node(col, "all"), bB)
        end
        h5open(fn) do f
            @test attrs(f["VTKHDF/A/Steps"])["NSteps"] == 3
            @test attrs(f["VTKHDF/B/Steps"])["NSteps"] == 3
        end
        if HAS_VTK
            d = vtkdump(fn)
            @test d["time_steps"] == 1:3
            @test d["steps"][3]["blocks"][1]["point_data"]["u"] == fill(3, 8)
            @test d["steps"][2]["blocks"][2]["point_data"]["w"] == fill(4, 4)
        end
    end

    @testset "mismatched block time values throw" begin
        fn = joinpath(dir, "tbad.vtkhdf")
        col = vtkhdf_collection(fn)
        bA = vtkhdf_grid(col, "A", cube, hex; temporal = true)
        bB = vtkhdf_grid(col, "B", sq, polys; temporal = true)
        write_timestep(f -> (f["u"] = rand(8)), bA, 1.0)
        write_timestep(f -> (f["w"] = rand(4)), bB, 2.0)
        @test_throws Exception close(col)
    end
end
