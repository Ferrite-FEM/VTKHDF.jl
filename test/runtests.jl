using WriteVTKHDF
using ParallelTestRunner

const testsuite = find_tests(@__DIR__)
delete!(testsuite, "setup")  # shared helpers, not a test file

runtests(
    WriteVTKHDF, ARGS;
    testsuite,
    init_code = :(include($(joinpath(@__DIR__, "setup.jl"))))
)
