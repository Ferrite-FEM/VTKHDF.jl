@testset "PolyData" begin
    dir = mktempdir()
    # unit square + diagonal line + two vertices
    pts = Float64[
        0 1 1 0
        0 0 1 1
        0 0 0 0
    ]
    polys = [MeshCell(PolyData.Polys(), [1, 2, 3, 4])]
    lines = [MeshCell(PolyData.Lines(), [1, 3])]
    verts = [MeshCell(PolyData.Verts(), [1]), MeshCell(PolyData.Verts(), [2])]

    @testset "static" begin
        fn = joinpath(dir, "pd.vtkhdf")
        vtkhdf_grid(fn, pts, polys, lines, verts) do vtk
            # npoints == ncells here, so the location must be explicit
            @test_throws ErrorException vtk["h"] = Float64.(1:4)
            vtk["h", VTKPointData()] = Float64.(1:4)
            # cell data in Vertices, Lines, Polygons, Strips order
            vtk["c", VTKCellData()] = Float64[1, 2, 3, 4]
        end
        h5open(fn) do f
            g = f["VTKHDF"]
            @test attrs(g)["Type"] == "PolyData"
            for cat in ("Vertices", "Lines", "Polygons", "Strips")
                @test haskey(g, cat)
                @test length(read(g["$cat/NumberOfCells"])) == 1
            end
            @test read(g["Vertices/NumberOfCells"]) == [2]
            @test read(g["Lines/NumberOfCells"]) == [1]
            @test read(g["Polygons/NumberOfCells"]) == [1]
            @test read(g["Strips/NumberOfCells"]) == [0]
            @test read(g["Polygons/Connectivity"]) == [0, 1, 2, 3]
            @test read(g["Lines/Connectivity"]) == [0, 2]
            @test read(g["Strips/Offsets"]) == [0]
        end
        if HAS_VTK
            d = vtkdump(fn)
            s = only(d["steps"])
            @test s["class"] == "vtkPolyData"
            @test s["npoints"] == 4 && s["ncells"] == 4
            @test s["points"] == nested(pts)
            @test s["point_data"]["h"] == 1:4
            @test s["cell_data"]["c"] == 1:4
        end
    end

    @testset "mixed categories in one vector are rejected" begin
        fn = joinpath(dir, "pd_mixed.vtkhdf")
        @test_throws ArgumentError vtkhdf_grid(
            fn, pts,
            [MeshCell(PolyData.Polys(), [1, 2, 3, 4]), MeshCell(PolyData.Lines(), [1, 2])]
        )
    end

    @testset "temporal" begin
        fn = joinpath(dir, "pd_t.vtkhdf")
        vtk = vtkhdf_grid(fn, pts, polys, lines; temporal = true)
        for s in 1:3
            write_timestep(vtk, Float64(s)) do frame
                frame["w"] = fill(Float64(s), 4)
            end
        end
        close(vtk)
        h5open(fn) do f
            g = f["VTKHDF"]
            @test size(read(g["Steps/CellOffsets"])) == (4, 3)
            @test read(g["Steps/PointDataOffsets/w"]) == [0, 4, 8]
            @test length(g["Points"]) == 3 * 4  # HDF5.length = prod(dims); geometry once
        end
        if HAS_VTK
            d = vtkdump(fn)
            @test d["time_steps"] == 1:3
            @test [st["point_data"]["w"] for st in d["steps"]] == [fill(s, 4) for s in 1:3]
        end
    end
end
