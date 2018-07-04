4/26/2018

Julia implementation of the STARRY/limbdark equations, with
derivatives.

So far the s_n(r,b) vector has been computed and
auto-diffed to give the jacobian of s_n with respect
to r and b up to l_max.  To run this from the Julia prompt,
go to ../tests/ directory and run:

julia> include("test_sn_jacobian.jl")

Routines in this directory:

related to limbdark paper:
- transit_poly.jl:  Computes the polynomial limb-darkening transits
  with derivatives in terms of r, b, and u_n/c_n.
- cel_bulirsch.jl:  Computes the generalized complete elliptic integral
  from Bulirsch (1969).
- IJv_derivative.jl:  Computes the I_v and J_v functions and derivatives
  for use with transit_poly.jl.
- s2.jl: Computes the all important linear limb-darkening function.
- area_triangle.jl: Accurate algorithm for computing the area of a triangle.

related to STARRY paper:
- sn.jl: Computes the s_n functions from the starry paper, which
  are a basis set for transit/occultation of spherical harmonics.
- sn_jacobian.jl:  Computes Jacobian of s_n with respect to r, b.
