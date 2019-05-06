10/11/2018

Julia implementation of the STARRY/limbdark equations, with
derivatives.

Quick start:  to compute the flux of a star with a transiting
planet of radius ratio `r = 0.1` and impact parameter `b=0.5`,
with limb-darkening coefficients `u_n = [0.1,0.2,0.3]`, without
derivatives, call:

```
julia> include("transit_poly_struct.jl")
julia> f = transit_poly(0.1,0.5,[0.1,0.2,0.3])
0.9891353686760945

```
Although this is the quickest entry into the code, it is
*not* the fastest when computing an entire light curve.
