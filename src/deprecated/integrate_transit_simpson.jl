# This is code for computing a transit model and derivatives integrated over
# a time step, giving the fluence in units of time (since flux is normalized to unity).

using Cubature
include("../transit_poly_struct.jl")
include("simpson.jl")

# Now the version with derivatives:
function integrate_timestep_gradient!(param::Array{T,1},trans::Transit_Struct{T},t1::T,t2::T,tol::T,maxdepth::Int64) where {T <: Real}
  # Keep track of maximum number of depths for recursion:
  depthmax = 0
  # Define third to save time in division:
  third = inv(3)
  # Keep track of how many evaluations of the transit function are carried out:
  neval = 0
  
  # This routine gives a simple model for the impact parameter versus time:
  function solver(time::T) where {T <: Real}
  # This predicts the impact parameter as a function of time in
  # terms of the transit_struct.
  # For now, transit is approximated as a straight line, where params are given as:
  # param[1] = t0 (mid-transit time; units of days/JD)
  # param[2] = v/R_* (sky velocity in terms of radius of the star - units are day^{-1}
  # param[3] = b (impact parameter at time t0)
  return sqrt(param[3]^2+(param[2]*(time-param[1]))^2)
  end

  # Function to define the vector integration in integrator:
#  function transit_flux(tmid::T,fmid::T) where {T <: Real}
  function transit_flux(tmid::T) where {T <: Real}
  # Add the impact parameter to the transit structure:
  trans.b = solver(tmid)
  # Carry out the transit computation:
  fmid = transit_poly_d!(trans); neval += 1
  return fmid
  end

  # Function to define the vector integration in cubature:
  function transit_flux(tmid::T,fmid::T) where {T <: Real}
  trans.b = solver(tmid)
  fmid = transit_poly_d!(trans)
  return fmid
  end

#  fint,ferr = hquadrature(trans.n+6,transit_flux,t1,t2,abstol=tol)
  # fint contains time-integrated flux (fluence) and derivatives:
# simpson(a::T, b::T, f::Function, I_of_f::Array{T,1}, i::T, eps::T, N::Int64, nf::Int64) where {T <: Real}
  tend = zero(T)
  fint = zero(T)
  fint = simpson(t1,t2,transit_flux,fint,tend,tol,maxdepth)
  fint2,ferr = hquadrature(transit_flux,t1,t2,abstol=tol)
  if abs(fint-fint2) > 1e-8
    println("itg, r: ",trans.r," b: ",trans.b," f/df: ",fint," tol: ",tol," t2-t1: ",t2-t1)
    println("cub, r: ",trans.r," b: ",trans.b," f/df: ",fint2)
    println("f[t1]: ",transit_flux(t1)," f[t2]: ",transit_flux(t2))
    read(STDIN,Char)
  end
  # Return the number of evalutions and maximum depth for record-keeping:
return fint::T,neval::Int64,depthmax::Int64
end
