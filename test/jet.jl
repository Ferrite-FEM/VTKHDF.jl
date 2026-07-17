using JET
using Test
using WriteVTKHDF

@testset "JET.jl" begin
    JET.test_package(WriteVTKHDF; target_modules = (WriteVTKHDF,), toplevel_logger = nothing)
end

@testset "JET entry points" begin
    points = rand(3, 8)
    cells = [MeshCell(VTKCellTypes.VTK_HEXAHEDRON, 1:8)]
    JET.test_call(
        vtkhdf_grid, (String, Matrix{Float64}, Vector{eltype(cells)});
        target_modules = (WriteVTKHDF,)
    )
    JET.test_call(vtkhdf_table, (String,); target_modules = (WriteVTKHDF,))
    JET.test_call(vtkhdf_collection, (String,); target_modules = (WriteVTKHDF,))
end
