# Table was added in spec 2.8; VTK's reader may not support it yet
# (VTK 9.6 predates it), so validation is structural.
@testset "Table" begin
    dir = mktempdir()

    @testset "static" begin
        fn = joinpath(dir, "tbl.vtkhdf")
        a = rand(5)
        vtkhdf_table(fn) do tbl
            tbl["a"] = a
            tbl["b"] = collect(1:5)
            tbl["m"] = rand(3, 5)   # multi-component column
        end
        h5open(fn) do f
            g = f["VTKHDF"]
            @test attrs(g)["Type"] == "Table"
            @test attrs(g)["Version"] == [2, 8]
            @test read(g["NumberOfRows"]) == [5]
            @test read(g["RowData/a"]) == a
            @test read(g["RowData/b"]) == 1:5
        end
        HAS_H5DUMP && @test h5dump_shape(fn, "/VTKHDF/RowData/m") == (5, 3)
    end

    @testset "inconsistent columns throw" begin
        fn = joinpath(dir, "tbl_bad.vtkhdf")
        tbl = vtkhdf_table(fn)
        tbl["a"] = rand(5)
        tbl["b"] = rand(4)
        @test_throws Exception close(tbl)
    end

    @testset "temporal" begin
        fn = joinpath(dir, "tbl_t.vtkhdf")
        vtk = vtkhdf_table(fn; temporal = true)
        for s in 1:3
            write_timestep(vtk, Float64(s)) do frame
                frame["a"] = fill(Float64(s), 4 + s)
            end
        end
        close(vtk)
        h5open(fn) do f
            g = f["VTKHDF"]
            @test read(g["NumberOfRows"]) == [5, 6, 7]
            @test read(g["Steps/RowDataOffsets/a"]) == [0, 5, 11]
            @test read(g["Steps/Values"]) == 1:3
        end
    end
end
