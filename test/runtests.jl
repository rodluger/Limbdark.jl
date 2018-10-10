using Test
include("loglinspace.jl")
include("test_cel_vec.jl")
include("test_s2_gradient.jl")
include("test_transit_poly.jl")

#include("test_IJv_derivative.jl")
#include("test_sn_jacobian2.jl")

skip_plots = false
include("test_transit_poly_gradient2.jl")
