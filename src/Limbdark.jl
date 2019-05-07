"""
    Limbdark

Fast, analytic, and numerically stable transit light curves for stars with
arbitrary order limb darkening.

"""
module Limbdark

include("transit_structure.jl")
include("transit_poly_struct.jl")
include("integrate_lightcurve.jl")

export transit_init, transit_poly, integrate_lightcurve

end