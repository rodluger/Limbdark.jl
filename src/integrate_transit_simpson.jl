# This is code for computing a transit model and derivatives integrated over
# a time step, giving the fluence in units of time (since flux is normalized to unity).

include("transit_poly_struct.jl")

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
    fmid[5+i] = trans_mid.dfdd[i]  # d_i derivatives
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
#  function transit_flux_derivative(tmid::T,fmid::Array{T,1}) where {T <: Real}
  function transit_flux_derivative(tmid::T) where {T <: Real}
  if VERSION >=  v"0.7"
    fmid = Array{T}(undef,6+trans.n)
  else
    fmid = Array{T}(6+trans.n)
  end
  # Add the impact parameter to the transit structure:
  trans.b = solver(tmid)
  # Carry out the transit computation:
  fmid0 = transit_poly_d!(trans); neval += 1
  # Then, fillin the flux and derivatives in the fmid vector:
  fill_flux!(tmid,fmid0,trans,fmid)
#  println("b: ",trans.b," f/df: ",fmid)
  return fmid
  end

  # Inner integrator which uses adaptive Simpson's method thanks to DFM:
  function inner(func::Function,x0::T,dx::T,ym::Array{T,1},y0::Array{T,1},yp::Array{T,1},tol::T,max_depth::Int64,min_depth::Int64,depth::Int64) where {T <: Real}
  # For each depth we need a left and a right copy of: x_m, val_m; x_p, val_p, int_mp, pred
  x_m = x0 - 0.5*dx; val_m = func(x_m)
  x_p = x0 + 0.5*dx; val_p = func(x_p)
  int_mp = 0.5*dx * (2*y0+ ym + yp + 4*(val_m+val_p))*third
  pred   =     dx * (4*y0+ ym + yp)*third
  if (depth < min_depth || (depth < max_depth && maximum(abs,pred-int_mp)> tol))
    return inner(func, x_m, 0.5*dx, ym, val_m, y0, tol, max_depth, min_depth, depth+1) +
           inner(func, x_p, 0.5*dx, y0, val_p, yp, tol, max_depth, min_depth, depth+1)
  end
  return int_mp
  end

#  # Inner integrator (this version avoids broadcast operations; however, we don't need this for now).
#  function inner(func::Function,x0::T,dx::T,ym::Array{T,1},y0::Array{T,1},yp::Array{T,1},tol::T,max_depth::Int64,min_depth::Int64,depth::Int64) where {T <: Real}
#  if VERSION >= v"0.7"
#    int_p = Array{T}(undef,6+trans.n); int_m = Array{T}(undef,6+trans.n); pred = Array{T}(undef,6+trans.n)
#  else
#    int_p = Array{T}(6+trans.n); int_m = Array{T}(6+trans.n); pred = Array{T}(6+trans.n)
#  end
#  x_m = x0 - 0.5*dx;
#  val_m = func(x_m);
#  x_p = x0 + 0.5*dx;
#  val_p = func(x_p);
#  maxdiff = zero(T)
#  @inbounds for i=1:6+trans.n
#    int_m[i] = 0.5*dx * (4*val_m[i] + ym[i] + y0[i])*third 
#    int_p[i] = 0.5*dx * (4*val_p[i] + yp[i] + y0[i])*third
#    pred[i]  =     dx * (4*y0[i]    + ym[i] + yp[i])*third
#    diffabs = abs(pred[i]-int_m[i]-int_p[i])
#    if  diffabs > maxdiff
#      maxdiff = diffabs
#    end
#  end
#  if (depth < min_depth || (depth < max_depth && maxdiff > tol))
#    int_m = inner(func, x_m, 0.5*dx, ym, val_m, y0, tol, max_depth, min_depth, depth+1);
#    int_p = inner(func, x_p, 0.5*dx, y0, val_p, yp, tol, max_depth, min_depth, depth+1);
#  end
#  @inbounds for i=1:6+trans.n
#     int_m[i] += int_p[i]
#  end
#  return int_m::Array{T,1}
#  end

  # Outer integrators - this dispatch computes function at upper and lower limits:
  function outer!(func::Function,lower::T,upper::T,tol::T;max_depth=50,min_depth=0) where {T <: Real}
  ym = func(lower)
  yp = func(upper)
  return outer!(func,lower,upper,ym,yp,tol,max_depth,min_depth)
  end
 
  # This version first compares trapzoid and Simpson's rule - if they are not precise
  # enough in comparison, then it begins refinement in inner recursion:
  function outer!(func::Function,lower::T,upper::T,ym::Array{T,1},yp::Array{T,1},tol::T,max_depth::Int64,min_depth::Int64) where {T <: Real}
  x0 = 0.5 * (upper + lower);
  dx = 0.5 * (upper - lower);
  y0 = func(x0);
  if VERSION >= v"0.7"
    int_trap = Array{T}(undef,6+trans.n)
    int_simp = Array{T}(undef,6+trans.n)
  else
    int_trap = Array{T}(6+trans.n)
    int_simp = Array{T}(6+trans.n)
  end
  maxdiff = zero(T)
  @inbounds for i=1:6+trans.n
    int_trap[i] = 0.5 * dx * (ym[i] + 2 * y0[i] + yp[i]);
    int_simp[i] = dx * (ym[i] + 4 * y0[i] + yp[i]) / 3;
    diffabs = int_trap[i]-int_simp[i]
    if diffabs > maxdiff
      maxdiff = diffabs
    end
  end 
  if (min_depth == 0 && maxdiff < tol)
    return int_simp;
  end
  return inner(func, x0, dx, ym, y0, yp, tol, max_depth, min_depth, 0);
  end

#  fint,ferr = hquadrature(trans.n+6,transit_flux_derivative,t1,t2,abstol=tol)
  # fint contains time-integrated flux (fluence) and derivatives:
  fint .= outer!(transit_flux_derivative,t1,t2,tol)
  # Return the number of evalutions and maximum depth for record-keeping:
return neval::Int64,depthmax::Int64
end
