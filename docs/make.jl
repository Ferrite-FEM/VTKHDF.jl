using Documenter
using WriteVTKHDF

makedocs(
    sitename = "WriteVTKHDF.jl",
    modules = [WriteVTKHDF],
    pages = [
        "Home" => "index.md",
        "Manual" => "manual.md",
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
