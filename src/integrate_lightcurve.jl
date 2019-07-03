
# This computes a straight-line transit lightcurve, integrated over a timestep.

include("integrate_transit_simpson_vec.jl")

function integrate_lightcurve!(trans::Transit_Struct{T},param::Array{T,1},t::Array{T,1},dt::T,favg1::Array{T,2},nt::Int64,tol::T,maxdepth::Int64,neval_t::Array{Int64,1},depthmax::Array{Int64,1}) where {T <: Real}
# This routine integrates a lightcurve with a constant velocity, v,  
# and impact parameter, b0, which are specified in the "param" vector.
# The trans structure needs to be initialized before calling this routine,
# which will contain workspace for the transit computation, the radius ratio,
# and the limb-darkening parameters.

# Input quantities:
# trans:  structure which holds the parameters describing the transit
# param:  parameters specifying the orbital path of the transit (straight line for now)
#     t:  vector of transit times
#    dt:  exposure time
# favg1:  vector which contains the time-averaged flux and derivatives
#    nt:  number of timesteps
#   tol:  tolerance for computing the integration
# maxdepth:  maximum depth allowed for the adaptive Simpson integration
# neval_t:  number of evaluations for each time step
# depthmax:  maximum depth reached for each timestep
#
# Upon return, favg1, neval_t, and depthmax will be updated.
#
# First, find the points of contact:
b0 = abs(param[3]); r = trans.r; t0 = param[1]; v = param[2]
# Check to see whether the two bodies overlap:
if b0 > (1.0+r)
  println("b0: ",b0," r: ",r)
  println("No transit")
  # No transit occurs, so returns ones for lightcurves, and zeros for derivatives:
#  favg1[1,:] = one(T)
  favg1[1,:] = zero(T)
  favg1[2:5+trans.n,:] = zero(T)
  return
# If a grazing transit, then there are only two points of contact:
elseif (b0+r) > 1.0
  # Times of two points of contact:
  tc = [t0 - sqrt((1.0+r)^2-b0^2)/v,t0 + sqrt((1.0+r)^2-b0^2)/v]
  nc = 2
# If a full transit, then there are four points of contact:
else
  # Times of four points of contact:
  tc = [t0 - sqrt((1.0+r)^2-b0^2)/v,t0 - sqrt((1.0-r)^2-b0^2)/v,t0 + sqrt((1.0-r)^2-b0^2)/v,t0 + sqrt((1.0+r)^2-b0^2)/v]
  nc = 4
end
# Invers of time step:
dtinv = inv(dt)
# Set up array to hold the integrated time step:
fint  = zeros(T,6+trans.n)
# Loop over the number of time steps:
@inbounds for i=1:nt
  # Keep track of the number of evaluations and the maximum depth:
  neval_t[i] = 0
  depthmax[i] = 0
  # Specify the starting and final timestep endpoints:
  t1 = t[i]-0.5*dt ; t2 = t[i]+0.5*dt
  # Check to see whether these midpoints contain the contact times:
  if t2 < tc[1] || t1 > tc[nc]
    # No points lie within the transit:
    ftmp=zeros(T,6+trans.n)
#    ftmp[1]=one(T)
    ftmp[1]=zero(T)
  else
    # Loop over points of contact:
    tlim = [t1]
    for j=1:nc
      if t1 < tc[j] && tc[j] < t2
        push!(tlim,tc[j])
      end
    end
    push!(tlim,t2)
    # ftmp is an array which contains the derivatives with respect to g_n:
    ftmp=zeros(T,6+trans.n)
    # Carry out the integration between the start & end time of the exposure,
    # as well as up to each contact point:
    for j=1:length(tlim)-1
      nevali,depthmaxi =  integrate_timestep_gradient!(param,trans,tlim[j],tlim[j+1],tol,maxdepth,fint)
#      println("r: ",trans.r," b: ",trans.b," fint: ",fint," dt: ",tlim[j+1]-tlim[j])
      neval_t[i] += nevali
      if depthmaxi > depthmax[i]
        depthmax[i] = depthmaxi
      end
      ftmp += fint
    end
    # Now divide by the exposure time to get the average flux over the integration:
    ftmp *= dtinv
#    println("ftmp: ",ftmp)
  end
  # Pass the flux and derivatives with respect to the lightcurve variables to
  # the average flux vector:
  favg1[1:5,i]=ftmp[1:5]
  # Convert from g_n to u_n derivatives for the remaining parts of the vector:
  favg1[6:5+trans.n,i]=BLAS.gemv('T',1.0,trans.dgdu,ftmp[6:6+trans.n])
#  println("i: ",i," t: ",t[i]," result: ",ftmp)
end
return
end
