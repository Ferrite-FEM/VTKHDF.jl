using Test

# JET does not support Julia nightlies (it loads as an empty stub there), and
# on Julia 1.10 only pre-0.10 JET versions resolve; run the analysis on
# supported stable/pre-release versions only.
const RUN_JET = v"1.11" <= VERSION && !occursin("DEV", string(VERSION))

if !RUN_JET
    @info "Skipping JET analysis on Julia $VERSION"
else
    using JET
    using VTKHDF

    @testset "JET.jl" begin
        JET.test_package(VTKHDF; target_modules = (VTKHDF,), toplevel_logger = nothing)
    end

    @testset "JET entry points" begin
        points = rand(3, 8)
        cells = [MeshCell(VTKCellTypes.VTK_HEXAHEDRON, 1:8)]
        JET.test_call(
            vtkhdf_grid, (String, Matrix{Float64}, Vector{eltype(cells)});
            target_modules = (VTKHDF,)
        )
        JET.test_call(vtkhdf_table, (String,); target_modules = (VTKHDF,))
        JET.test_call(vtkhdf_collection, (String,); target_modules = (VTKHDF,))
    end

    @testset "JET reader entry points" begin
        UGReader = VTKHDF.VTKHDFReader{VTKHDF.ReadUnstructured}
        JET.test_call(vtkhdf_open, (String,); target_modules = (VTKHDF,))
        JET.test_call(getindex, (UGReader, String); target_modules = (VTKHDF,))
        JET.test_call(read_timestep, (UGReader, Int); target_modules = (VTKHDF,))
        JET.test_call(read_points, (UGReader,); target_modules = (VTKHDF,))
        JET.test_call(read_cells, (UGReader,); target_modules = (VTKHDF,))
        JET.test_call(
            read_coordinates, (VTKHDF.VTKHDFReader{VTKHDF.ReadRectilinear},);
            target_modules = (VTKHDF,)
        )
        JET.test_call(
            getindex, (VTKHDF.VTKHDFCollectionReader, String);
            target_modules = (VTKHDF,)
        )
        JET.test_call(close, (UGReader,); target_modules = (VTKHDF,))
    end
end
