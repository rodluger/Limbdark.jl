if VERSION >= v"0.7"
  using Test
  using LinearAlgebra
else
  using Base.Test
  include("random.jl")
end
include("../src/define_constants.jl")
include("loglinspace.jl")
include("test_cel_vec.jl")
include("test_s2_gradient.jl")
include("test_transit_poly.jl")

#include("test_IJv_derivative.jl")
#include("test_Mn.jl")
#include("test_sn_jacobian2.jl")

skip_plots = false
#skip_plots = true
include("test_transit_poly_gradient2.jl")
include("test_integrate_transit_simpson_vec.jl")


