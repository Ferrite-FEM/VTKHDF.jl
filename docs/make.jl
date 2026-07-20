using Documenter
using Literate
using VTKHDF

# One example per reference file shown in the VTKHDF specification.
const EXAMPLES = [
    "image_data.jl",
    "unstructured_grid.jl",
    "poly_data.jl",
    "overlapping_amr.jl",
    "partitioned_collection.jl",
    "temporal_poly_data.jl",
]

for ex in EXAMPLES
    # `@example` blocks: Documenter runs the example during the docs build,
    # so the reading sections show their actual output.
    Literate.markdown(
        joinpath(@__DIR__, "..", "examples", ex), joinpath(@__DIR__, "src", "examples");
        credit = false, codefence = "```@example $(splitext(ex)[1])" => "```"
    )
end

prettyurls = get(ENV, "CI", nothing) == "true"

# The example's title is its first Literate markdown line: `# # <title>`.
example_title(ex) = strip(chopprefix(readline(joinpath(@__DIR__, "..", "examples", ex)), "# #"))

# Gallery of clickable cards. Documenter does not rewrite paths in raw HTML,
# so they depend on whether the page ends up as overview.html or overview/index.html.
open(joinpath(@__DIR__, "src", "examples", "overview.md"), "w") do io
    println(io, "# Examples\n")
    println(io, "One example per reference file shown in the VTKHDF specification, ported to this package. Click an example to see the code that writes (and reads back) the file.\n")
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
    sitename = "VTKHDF.jl",
    modules = [VTKHDF],
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
        canonical = "https://ferrite-fem.github.io/VTKHDF.jl",
        assets = ["assets/custom.css"],
    ),
    checkdocs = :exports,
)

deploydocs(
    repo = "github.com/Ferrite-FEM/VTKHDF.jl.git",
    push_preview = true,
)
