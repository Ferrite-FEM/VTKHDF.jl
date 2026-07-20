using Documenter
using Literate
using WriteVTKHDF

# One example per dataset type of the VTKHDF specification's reference files.
const EXAMPLES = [
    "image_data.jl",
    "unstructured_grid.jl",
    "poly_data.jl",
    "overlapping_amr.jl",
    "partitioned_collection.jl",
    "temporal_poly_data.jl",
]

for ex in EXAMPLES
    Literate.markdown(
        joinpath(@__DIR__, "..", "examples", ex), joinpath(@__DIR__, "src", "examples");
        credit = false, codefence = "```julia" => "```"
    )
end

makedocs(
    sitename = "WriteVTKHDF.jl",
    modules = [WriteVTKHDF],
    pages = [
        "Home" => "index.md",
        "Manual" => "manual.md",
        "Examples" => ["examples/$(splitext(ex)[1]).md" for ex in EXAMPLES],
        "API" => "api.md",
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://kristofferc.github.io/WriteVTKHDF.jl",
    ),
    checkdocs = :exports,
)

deploydocs(
    repo = "github.com/KristofferC/WriteVTKHDF.jl.git",
    push_preview = true,
)
