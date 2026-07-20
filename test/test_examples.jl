# Run every script in examples/ and read the results back through VTK.
@testset "examples" begin
    exdir = joinpath(dirname(@__DIR__), "examples")
    scripts = sort(filter(endswith(".jl"), readdir(exdir)))
    mktempdir() do dir
        cd(dir) do
            for f in scripts
                @testset "$f" begin
                    Base.include(Module(gensym(f)), joinpath(exdir, f))
                end
            end
        end
        files = filter(endswith(".vtkhdf"), readdir(dir; join = true))
        @test length(files) == length(scripts)
        if HAS_VTK
            for fn in files
                @test !isempty(vtkdump(fn)["steps"])
            end
        end
    end
end
