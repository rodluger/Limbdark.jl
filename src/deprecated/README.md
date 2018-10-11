7/4/2018

This is a graveyard of code which was used in development
of STARRY/limbdark trying out:

- different ways of computing the H_uv function
- different ways of computing M_pq
- s2 (linear limb-darkening) in terms of standard complete
  elliptic integrals
- derivatives of standard elliptic integrals
- estimating light travel time delays

sn.jl:  Computed the s_n functions needed for starry (a different
  set of Green's functions is used in limbdark).
sn_jacobian.jl:  Computed the derivatives of s_n with auto-diff.
area_triangle.jl:  Computes are of triangle using Kahan method.
