# This is code for computing a transit model and derivatives integrated over
# a time step, giving the fluence in units of time (since 
# flux is normalized to unity).

using Cubature

include("transit_poly_struct.jl")

# Now the version with derivatives:
function integrate_timestep_gradient!(param::Array{T,1},trans::Transit_Struct{T},t1::T,t2::T,tol::T,maxdepth::Int64) where {T <: Real}

  function fill_flux!(tmid::T,fmid0::T,trans_mid::Transit_Struct{T}) where {T <: Real}
  fmid = Array{T}(6+trans_mid.n)
  fmid[1] = fmid0   # flux
  fmid[2] = trans_mid.dfdrb[1]  # r derivative
  fmid[3] = trans_mid.dfdrb[2]/trans_mid.b*param[2]^2*(param[1]-tmid)  # t0 derivative
  fmid[4] = trans_mid.dfdrb[2]*param[2]/trans_mid.b*(tmid-param[1])^2  # v derivative
  fmid[5] = trans_mid.dfdrb[2]*param[3]/trans_mid.b  # b0 derivative
  @inbounds for i=1:trans_mid.n+1
    fmid[5+i] = trans_mid.dfdd[i]  # d_i derivatives
  end
  return fmid
  end

  function solver(time::T) where {T <: Real}
  # This predicts the impact parameter as a function of time in
  # terms of the transit_struct.
  # For now, transit is approximated as a straight line, where params are given as:
  # param[1] = t0 (mid-transit time; units of days/JD)
  # param[2] = v/R_* (sky velocity in terms of radius of the star - units are day^{-1}
  # param[3] = b (impact parameter at time t0)
  return sqrt(param[3]^2+(param[2]*(time-param[1]))^2)
  end

  # Function to define the vector integration in cubature:
  function transit_flux_derivative(tmid::T,fmid::Array{T,1}) where {T <: Real}
  trans.b = solver(tmid)
  fmid0 = transit_poly_d!(trans)
  fmid .= fill_flux!(tmid,fmid0,trans)
#  println("b: ",trans.b," f/df: ",fmid)
  return
  end
  
  fint,ferr = hquadrature(trans.n+6,transit_flux_derivative,t1,t2,abstol=tol)
return fint
end
