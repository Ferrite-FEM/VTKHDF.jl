using Aqua
using Test
using VTKHDF

@testset "Aqua.jl" begin
    Aqua.test_all(VTKHDF)
end
