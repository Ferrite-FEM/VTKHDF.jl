using Aqua
using Test
using WriteVTKHDF

@testset "Aqua.jl" begin
    Aqua.test_all(WriteVTKHDF)
end
