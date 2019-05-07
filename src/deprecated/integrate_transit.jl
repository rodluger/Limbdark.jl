# This is code for computing a transit model integrated over
# a time step, giving the fluence in units of time (since 
# flux is normalized to unity).

include("transit_poly_struct.jl")

# First the version without derivatives, which tracks the number of depths:
function integrate_timestep_track(param::Array{T,1},trans::Transit_Struct{T},time::T,dt::T,tol::T,maxdepth::Int64) where {T <: Real}

  function solver(param::Array{T,1},time::T) where {T <: Real}
  # This predicts the impact parameter as a function of time in
  # terms of the transit_struct.
  # For now, transit is approximated as a straight line, where params are given as:
  # param[1] = t0 (mid-transit time; units of days/JD)
  # param[2] = v/R_* (sky velocity in terms of radius of the star - units are day^{-1}
  # param[3] = b (impact parameter at time t0)
  return sqrt(param[3]^2+(param[2]*(time-param[1]))^2)
  end

  # Now, adopt Dan Foreman-Mackey's recursive integration approach:
  function integrate(param::Array{T,1}, f1::T, f2::T, t1::T, t2::T, depth::Int64) where {T <: Real}
  tmid = 0.5*(t1+t2);
  bmid = solver(param,tmid)
  trans.b = bmid
  fmid = transit_poly_d!(trans)
  fapprox = 0.5*(f1+f2)
  d = abs(fmid-fapprox)
  if d > tol && depth < maxdepth
     a,dmaxa,da = integrate(param,f1,fmid,t1,tmid,depth+1)
     b,dmaxb,db = integrate(param,fmid,f2,tmid,t2,depth+1)
     return a+b,maximum([dmaxa,dmaxb]),maximum([da,db])
  end
  return 0.5*(fapprox+fmid)*(t2-t1),depth,d
  end
  
  t1 = time - 0.5*dt
  b1 = solver(param,t1)
  trans.b = b1
  f1 = transit_poly_d!(trans)
  t2 = time + 0.5*dt
  b2 = solver(param,t2)
  trans.b = b2
  f2 = transit_poly_d!(trans)
  fint,dmax,diff = integrate(param,f1,f2,t1,t2,0)
return fint,dmax,diff
end

# First the version without derivatives, which doesn't track depth:
function integrate_timestep(param::Array{T,1},trans::Transit_Struct{T},time::T,dt::T,tol::T,maxdepth::Int64) where {T <: Real}

  function solver(param::Array{T,1},time::T) where {T <: Real}
  # This predicts the impact parameter as a function of time in
  # terms of the transit_struct.
  # For now, transit is approximated as a straight line, where params are given as:
  # param[1] = t0 (mid-transit time; units of days/JD)
  # param[2] = v/R_* (sky velocity in terms of radius of the star - units are day^{-1}
  # param[3] = b (impact parameter at time t0)
  return sqrt(param[3]^2+(param[2]*(time-param[1]))^2)
  end

  # Now, adopt Dan Foreman-Mackey's recursive integration approach:
  function integrate(param::Array{T,1}, f1::T, f2::T, t1::T, t2::T, depth::Int64) where {T <: Real}
  tmid = 0.5*(t1+t2);
  bmid = solver(param,tmid)
  trans.b = bmid
  fmid = transit_poly_d!(trans)
  fapprox = 0.5*(f1+f2)
  d = abs(fmid-fapprox)
  if d > tol && depth < maxdepth
     a = integrate(param,f1,fmid,t1,tmid,depth+1)
     b = integrate(param,fmid,f2,tmid,t2,depth+1)
     return a+b
  end
  return 0.5*(fapprox+fmid)*(t2-t1)
  end
  
  t1 = time - 0.5*dt
  b1 = solver(param,t1)
  trans.b = b1
  f1 = transit_poly_d!(trans)
  t2 = time + 0.5*dt
  b2 = solver(param,t2)
  trans.b = b2
  f2 = transit_poly_d!(trans)
  fint = integrate(param,f1,f2,t1,t2,0)
return fint
end

# Now the version with derivatives:
function integrate_timestep_gradient(param::Array{T,1},trans::Transit_Struct{T},time::T,dt::T,tol::T,maxdepth::Int64) where {T <: Real}

  function solver(param::Array{T,1},time::T) where {T <: Real}
  # This predicts the impact parameter as a function of time in
  # terms of the transit_struct.
  # For now, transit is approximated as a straight line, where params are given as:
  # param[1] = t0 (mid-transit time; units of days/JD)
  # param[2] = v/R_* (sky velocity in terms of radius of the star - units are day^{-1}
  # param[3] = b (impact parameter at time t0)
  return sqrt(param[3]^2+(param[2]*(time-param[1]))^2)
  end

  # Now, adopt Dan Foreman-Mackey's recursive integration approach:
  function integrate(param::Array{T,1}, f1::Array{T,1}, f2::Array{T,1}, t1::T, t2::T, depth::Int64) where {T <: Real}
  tmid = 0.5*(t1+t2);
  bmid = solver(param,tmid)
  trans.b = bmid
  fmid0 = transit_poly_d!(trans)
  # Compute: flux, df/dr, df/dt0, df/dv, df/db0, df/d d_n:
  fmid = [fmid0;trans.dfdrb[1];trans.dfdrb[2]/trans.b*param[2]*(param[1]-tmid);
     trans.dfdrb[2]*param[2]/trans.b*(tmid-param[1])^2;trans.dfdrb[2]*param[3]/trans.b;trans.dfdg]
  # Extrapolate the flux at the boundaries to the midpoint:
  fapprox = 0.5*(f1+f2)
  # Tolerance is set by the flux term only:
  d = abs(fmid[1]-fapprox[1])
  if d > tol && depth < maxdepth
     a = integrate(param,f1,fmid,t1,tmid,depth+1)
     b = integrate(param,fmid,f2,tmid,t2,depth+1)
     return a+b
  end
  return 0.5*(fapprox+fmid)*(t2-t1)
  end
  
  t1 = time - 0.5*dt
  b1 = solver(param,t1)
  trans.b = b1
  f10 = transit_poly_d!(trans)
  f1 = [f10;trans.dfdrb;trans.dfdg]
  f1 = [f10;trans.dfdrb[1];trans.dfdrb[2]/trans.b*param[2]*(param[1]-t1);
     trans.dfdrb[2]*param[2]/trans.b*(t1-param[1])^2;trans.dfdrb[2]*param[3]/trans.b;trans.dfdg]
  t2 = time + 0.5*dt
  b2 = solver(param,t2)
  trans.b = b2
  f20 = transit_poly_d!(trans)
  f2 = [f20;trans.dfdrb[1];trans.dfdrb[2]/trans.b*param[2]^2*(param[1]-t2);
     trans.dfdrb[2]*param[2]/trans.b*(t2-param[1])^2;trans.dfdrb[2]*param[3]/trans.b;trans.dfdg]
  fint = integrate(param,f1,f2,t1,t2,0)
return fint
end
