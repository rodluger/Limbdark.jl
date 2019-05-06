push!(LOAD_PATH, "../src/")
using Documenter, Limbdark
include("../src/transit_structure.jl")

makedocs(sitename="Limbdark")