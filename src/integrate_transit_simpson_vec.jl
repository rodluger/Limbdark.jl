
# This is code for computing a transit model and derivatives integrated over
# a time step, giving the fluence in units of time (since flux is normalized to unity).

#using Cubature
include("transit_poly_struct.jl")
include("simpson_vec.jl")

# Now the version with derivatives:
function integrate_timestep_gradient!(param::Array{T,1},trans::Transit_Struct{T},t1::T,t2::T,tol::T,maxdepth::Int64,fint::Array{T,1}) where {T <: Real}
  # Keep track of maximum number of depths for recursion:
  depthmax = 0
  # Define third to save time in division:
  third = inv(3)
  # Keep track of how many evaluations of the transit function are carried out:
  neval = 0
  
  # When transit computation is carried out, we need to fill in the vector for the flux and derivatives
  # at the midpoint of the current refinement time range:
  function fill_flux!(tmid::T,fmid0::T,trans_mid::Transit_Struct{T},fmid::Array{T,1}) where {T <: Real}
#  fmid = Array{T}(undef,6+trans_mid.n)
  fmid[1] = fmid0   # flux
  fmid[2] = trans_mid.dfdrb[1]  # r derivative
  binv = inv(trans_mid.b)
  fmid[3] = trans_mid.dfdrb[2]*binv*param[2]^2*(param[1]-tmid)  # t0 derivative
  fmid[4] = trans_mid.dfdrb[2]*param[2]*binv*(tmid-param[1])^2  # v derivative
  fmid[5] = trans_mid.dfdrb[2]*param[3]*binv  # b0 derivative
  @inbounds for i=1:trans_mid.n+1
    fmid[5+i] = trans_mid.dfdg[i]  # d_i derivatives
  end
  return fmid
  end

  function fill_flux!(tmid::T,fmid0::T,trans_mid::Transit_Struct{T}) where {T <: Real}
  fmid = Array{T}(6+trans_mid.n)
  binv = inv(trans_mid.b)
  fmid[1] = fmid0   # flux
  fmid[2] = trans_mid.dfdrb[1]  # r derivative
  fmid[3] = trans_mid.dfdrb[2]*binv*param[2]^2*(param[1]-tmid)  # t0 derivative
  fmid[4] = trans_mid.dfdrb[2]*param[2]*binv*(tmid-param[1])^2  # v derivative
  fmid[5] = trans_mid.dfdrb[2]*param[3]*binv  # b0 derivative
  @inbounds for i=1:trans_mid.n+1
    fmid[5+i] = trans_mid.dfdg[i]  # d_i derivatives
  end
  return fmid
  end

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
#  function transit_flux_derivative!(tmid::T,fmid::Array{T,1}) where {T <: Real}
  function transit_flux_derivative(tmid::T) where {T <: Real}
#  if VERSION >=  v"0.7"
#    fmid = Array{T}(undef,6+trans.n)
#  else
#    fmid = Array{T}(6+trans.n)  # This line is taking a lot of time
#  end
  # Add the impact parameter to the transit structure:
  trans.b = solver(tmid)
  # Carry out the transit computation:
  fmid0 = transit_poly_g!(trans); neval += 1
  # Then, fillin the flux and derivatives in the fmid vector:
#  fill_flux!(tmid,fmid0,trans,fmid)  # So is this line
  fill_flux!(tmid,fmid0-1,trans,fmid)  # So is this line
#  println("tfd, r: ",trans.r," b: ",trans.b," f/df: ",fmid)
  return
  end

  # Function to define the vector integration in cubature:
  function transit_flux_derivative!(tmid::T,fmid::Array{T,1}) where {T <: Real}
  trans.b = solver(tmid)
  fmid0 = transit_poly_g!(trans); neval += 1
#  fmid .= fill_flux!(tmid,fmid0,trans)
  fill_flux!(tmid,fmid0-1,trans,fmid)
#  println("b: ",trans.b," f/df: ",fmid)
  return
  end

#  fint,ferr = hquadrature(trans.n+6,transit_flux_derivative,t1,t2,abstol=tol)
  # fint contains time-integrated flux (fluence) and derivatives:
# simpson(a::T, b::T, f::Function, I_of_f::Array{T,1}, i::T, epsilon::T, N::Int64, nf::Int64) where {T <: Real}
  tend = zero(T)
  fint .= simpson_vec(t1,t2,transit_flux_derivative!,fint,tend,tol,maxdepth,trans.n+6)
#  println("fint: ",fint)
#  fint2,ferr = hquadrature(trans.n+6,transit_flux_derivative,t1,t2,abstol=tol)
#  if maximum(abs,fint-fint2) > 1e-8
#    println("itg, r: ",trans.r," b: ",trans.b," f/df: ",fint," tol: ",tol," t2-t1: ",t2-t1)
#    println("cub, r: ",trans.r," b: ",trans.b," f/df: ",fint2)
#    println("f[t1]: ",transit_flux_derivative(t1)," f[t2]: ",transit_flux_derivative(t2))
#    read(STDIN,Char)
#  end
  # Return the number of evalutions and maximum depth for record-keeping:
return neval::Int64,depthmax::Int64
end
