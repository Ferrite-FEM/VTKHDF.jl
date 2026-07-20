# Reading back Table files.

@testset "read table static" begin
    dir = mktempdir()
    fn = joinpath(dir, "tbl")
    b = Int32[4 5 6; 7 8 9]
    vtkhdf_table(fn) do tbl
        tbl["a"] = [1.0, 2.0, 3.0]
        tbl["b"] = b
        tbl["meta", VTKFieldData()] = "info"
    end
    vtkhdf_open(fn) do r
        @test VTKHDF.dataset_type(r) == "Table"
        @test VTKHDF.nrows(r) == 3
        @test r["a"] == [1.0, 2.0, 3.0]
        @test r["a", VTKRowData()] == [1.0, 2.0, 3.0]
        @test r["b"] == b
        @test r["b"] isa Matrix{Int32}
        @test r["meta"] == ["info"]
        @test keys(r, VTKRowData()) == ["a", "b"]
        @test_throws MethodError VTKHDF.npoints(r)  # tables have no points
    end

    # empty table
    fn = joinpath(dir, "empty")
    vtkhdf_table(fn) do tbl
    end
    vtkhdf_open(fn) do r
        @test VTKHDF.nrows(r) == 0
        @test isempty(keys(r, VTKRowData()))
    end
end

@testset "read table temporal" begin
    dir = mktempdir()
    fn = joinpath(dir, "ttbl")
    vtk = vtkhdf_table(fn; temporal = true)
    write_timestep(vtk, 0.0) do f
        f["a"] = [1.0, 2.0]
        f["b"] = Int64[10, 20]
    end
    write_timestep(vtk, 1.0) do f
        f["a"] = [3.0]
        f["b"] = Int64[30]
    end
    write_timestep(vtk, 2.0) do f
        f["a"] = [4.0, 5.0, 6.0]
        f["b"] = Int64[40, 50, 60]
    end
    close(vtk)

    vtkhdf_open(fn) do r
        @test VTKHDF.nsteps(r) == 3
        @test_throws ErrorException VTKHDF.nrows(r)
        @test VTKHDF.nrows(read_timestep(r, 1)) == 2
        @test VTKHDF.nrows(read_timestep(r, 2)) == 1
        @test VTKHDF.nrows(read_timestep(r, 3)) == 3
        @test read_timestep(r, 1)["a"] == [1.0, 2.0]
        @test read_timestep(r, 2)["a"] == [3.0]
        @test read_timestep(r, 3)["a"] == [4.0, 5.0, 6.0]
        @test read_timestep(r, 3)["b"] == Int64[40, 50, 60]
        @test_throws ErrorException r["a"]
    end
end
