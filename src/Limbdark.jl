"""
    Limbdark

Fast, analytic, and numerically stable transit light curves for stars with
arbitrary order limb darkening.

"""
module Limbdark

using SpecialFunctions, LinearAlgebra

export transit_init, transit_poly, transit_poly!, integrate_lightcurve!

include("transit_structure.jl")
include("Mn_coeff.jl")
include("Nn_coeff.jl")
include("compute_g_n_struct.jl")
include("define_constants.jl")
include("cel_bulirsch.jl")
include("s2.jl")
include("Mn_compute.jl")
include("Nn_compute.jl")
include("transit_poly_struct.jl")
include("simpson_vec.jl")
include("integrate_transit_simpson_vec.jl")
include("integrate_lightcurve.jl")

end