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

prettyurls = get(ENV, "CI", nothing) == "true"

# The example's title is its first Literate markdown line: `# # <title>`.
example_title(ex) = strip(chopprefix(readline(joinpath(@__DIR__, "..", "examples", ex)), "# #"))

# Gallery of clickable cards. Documenter does not rewrite paths in raw HTML,
# so they depend on whether the page ends up as overview.html or overview/index.html.
open(joinpath(@__DIR__, "src", "examples", "overview.md"), "w") do io
    println(io, "# Examples\n")
    println(io, "One example per dataset type of the VTKHDF specification, ported from its reference files. Click an example to see the code that writes the file.\n")
    println(io, "```@raw html")
    println(io, "<div class=\"example-gallery\">")
    for ex in EXAMPLES
        name = splitext(ex)[1]
        title = example_title(ex)
        href = prettyurls ? "../$name/" : "$name.html"
        assets = prettyurls ? "../../assets" : "../assets"
        println(io, "<a class=\"example-card\" href=\"$href\">")
        println(io, "  <img src=\"$assets/examples/$name-light.png\" alt=\"$title\">")
        println(io, "  <img src=\"$assets/examples/$name-dark.png\" alt=\"$title\">")
        println(io, "  <div class=\"example-card-title\">$title</div>")
        println(io, "</a>")
    end
    println(io, "</div>")
    println(io, "```")
end

makedocs(
    sitename = "WriteVTKHDF.jl",
    modules = [WriteVTKHDF],
    pages = [
        "Home" => "index.md",
        "Manual" => "manual.md",
        "Examples" => [
            "Overview" => "examples/overview.md",
            ["examples/$(splitext(ex)[1]).md" for ex in EXAMPLES]...,
        ],
        "API" => "api.md",
    ],
    format = Documenter.HTML(
        prettyurls = prettyurls,
        canonical = "https://kristofferc.github.io/WriteVTKHDF.jl",
        assets = ["assets/custom.css"],
    ),
    checkdocs = :exports,
)

deploydocs(
    repo = "github.com/KristofferC/WriteVTKHDF.jl.git",
    push_preview = true,
)
