# This is code for computing a transit model integrated over
# a time step.

# First the version without derivatives:
function integrate_timestep(param::Array{T,1},trans::Transit_Struct{T},time::T,dt::T,tol::T,maxdepth::Int64) where {T <: Real}

  function solver(param::Array{T,1},time::T) where {T <: Real}
  # This predicts the impact parameter as a function of time in
  # terms of the transit_struct.
  # For now, transit is approximated as a straight line, where params are given as:
  # param[1] = t0 (mid-transit time; units of days/JD)
  # param[2] = v/R_* (sky velocity in terms of radius of the star - units are day^{-1}
  # param[3] = r (radius-ratio)
  # param[4] = b (impact parameter at time t0)
  btime = sqrt(param[4]^2+(param[2]*(time-param[1]))^2)
  return 
  end

  # Now, adopt Dan Foreman-Mackey's recursive integration approach:
  function integrate(param::Array{T,1}, f1::T, f2::T, t1::T, t2::T, trans::Transit_Struct{T}) where {T <: Real}
  tmid = 0.5*(t1+t2);
  bmid = solver(param,tmid)
  trans.b = bmid
  fmid = transit_poly!(trans)
  fapprox = 0.5*(f1+f2)
  d = abs(fmid-fapprox)
  if d > tol && depth < maxdepth
     a = integrate(param,f1,fmid,t1,tmid,depth+1)
     b = integrate(param,fmid,f2,tmid,t2,depth+1)
     return a+b
  end
  return fapprox*(t2-t1)
  end
  
  t1 = time - 0.5*dt
  b1 = solver(param,t1)
  trans.b = b1
  f1 = transit_poly!(trans)
  t2 = time + 0.5*dt
  b2 = solver(param,t2)
  trans.b = b2
  f2 = transit_poly!(trans)
  fint = integrate(param,f1,f2,t1,t2,trans)
return fint
end

