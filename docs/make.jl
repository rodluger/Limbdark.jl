push!(LOAD_PATH, "../src/")
using Documenter, Limbdark
include("../src/transit_structure.jl")

makedocs(
    sitename="Limbdark",
    pages = [
        "Home" => "index.md",
        "Installation" => "install.md",
        "Quick start" => "quickstart.md",
        "API" => "api.md"
    ],
    format = Documenter.HTML(prettyurls = false)
)

