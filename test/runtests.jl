using Test, Limbdark, LinearAlgebra, Statistics, ForwardDiff, DiffResults, Random

include("loglinspace.jl")

include("test_cel_vec.jl")
include("test_s2_gradient.jl")
include("test_transit_poly.jl")
include("test_Mn.jl")

# Skip plots for CI tests
skip_plots = true
include("test_transit_poly_gradient2.jl")
include("test_integrate_transit_simpson_vec.jl")
