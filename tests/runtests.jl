using Base.Test

include("test_IJv_derivative.jl")
include("test_sn_jacobian2.jl")
include("test_s2_gradient.jl")
include("test_transit_poly.jl")
#skip_plots = true
skip_plots = false
include("test_transit_poly_gradient2.jl")
