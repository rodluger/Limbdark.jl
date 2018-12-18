# This is test of the code for computing a transit model integrated over a time step.
# Computes derivatives over the timestep.
using PyPlot

#include("../src/integrate_transit_cubature.jl")
include("../src/integrate_transit_simpson_vec.jl")

# Test it out:

include("../test/loglinspace.jl")

t1 = -1.5; t2 = 1.5; nt = 1000; dt = 0.3
t = zeros(nt)
t .= linearspace(t1,t2,nt)
#t = [t[1235]] ; nt =1
# The following compares different tolerances and maxdepths:
r = 0.1; b0 = 0.5; u_n = [0.3,0.3]; nu = length(u_n)
favg0 = zeros(nt,5+nu)
favg1 = zeros(nt,5+nu)
neval1= zeros(Int64,nt)
depthmax1= zeros(Int64,nt)
favg2 = zeros(nt,5+nu)
neval2= zeros(Int64,nt)
depthmax2= zeros(Int64,nt)
favg3 = zeros(nt,5+nu)
neval3= zeros(Int64,nt)
depthmax3= zeros(Int64,nt)
trans = transit_init(r,b0,u_n,true)
param = [0.0,1.0,b0]   # [t_0,v,b_0]

# Now, carry out speed test:

function compute_lightcurve!(trans::Transit_Struct{T},param::Array{T,1},t::Array{T,1},favg0::Array{T,2},nt::Int64) where {T <: Real}
@inbounds for i=1:nt
  trans.b = sqrt(param[3]^2+(param[2]*(t[i]-param[1]))^2)
  favg0[i,1] =  transit_poly!(trans)
  favg0[i,2] =  trans.dfdrb[1]
  favg0[i,3] =  trans.dfdrb[2]/trans.b*param[2]^2*(param[1]-t[i])
  favg0[i,4] =  trans.dfdrb[2]*param[2]/trans.b*(t[i]-param[1])^2
  favg0[i,5] =  trans.dfdrb[2]*param[3]/trans.b
  favg0[i,6:5+trans.n] = trans.dfdu
end
return
end

compute_lightcurve!(trans,param,t,favg0,nt)
@time compute_lightcurve!(trans,param,t,favg0,nt)

function integrate_lightcurve!(trans::Transit_Struct{T},param::Array{T,1},t::Array{T,1},dt::T,favg1::Array{T,2},nt::Int64,tol::T,maxdepth::Int64,neval_t::Array{Int64,1},depthmax::Array{Int64,1}) where {T <: Real}
# First, find the points of contact:
b0 = abs(param[3]); r = trans.r; t0 = param[1]; v = param[2]
if b0 > (1.0+r)
  println("No transit")
  favg1[:,1] = one(T)
  favg1[:,2:5+trans.n] = zero(T)
  return
elseif (b0+r) > 1.0
  # Two points of contact:
  tc = [t0 - sqrt((1.0+r)^2-b0^2)/v,t0 + sqrt((1.0+r)^2-b0^2)/v]
  nc = 2
else
  # Four points of contact:
  tc = [t0 - sqrt((1.0+r)^2-b0^2)/v,t0 - sqrt((1.0-r)^2-b0^2)/v,t0 + sqrt((1.0-r)^2-b0^2)/v,t0 + sqrt((1.0+r)^2-b0^2)/v]
  nc = 4
end
dtinv = inv(dt)
fint  = zeros(T,6+trans.n)
@inbounds for i=1:nt
  neval_t[i] = 0
  depthmax[i] = 0
  t1 = t[i]-0.5*dt ; t2 = t[i]+0.5*dt
  if t2 < tc[1] || t1 > tc[nc]
    # No points lie within the transit:
    ftmp=zeros(T,6+trans.n)
    ftmp[1]=one(T)
  else
    # Loop over points of conjunction:
    tlim = [t1]
    for j=1:nc
      if t1 < tc[j] && tc[j] < t2
        push!(tlim,tc[j])
      end
    end
    push!(tlim,t2)
    ftmp=zeros(T,6+trans.n)
    for j=1:length(tlim)-1
      nevali,depthmaxi =  integrate_timestep_gradient!(param,trans,tlim[j],tlim[j+1],tol,maxdepth,fint)
      neval_t[i] += nevali
      if depthmaxi > depthmax[i]
        depthmax[i] = depthmaxi
      end
      ftmp += fint
    end
    ftmp *= dtinv
#    println("ftmp: ",ftmp)
  end
  favg1[i,1:5]=ftmp[1:5]
  # Convert from d_n to u_n derivatives:
  favg1[i,6:5+trans.n]=BLAS.gemv('T',1.0,trans.dddu,ftmp[6:6+trans.n])
#  println("i: ",i," t: ",t[i]," result: ",ftmp)
end
return
end

#integrate_lightcurve!(trans,param,t,dt,favg1,nt,1e-4,3,neval1)
#@time integrate_lightcurve!(trans,param,t,dt,favg1,nt,1e-4,8,neval1)

ntol = 12
tol = logarithmspace(-13.0,-2.0,ntol)
K = 50
neval_mean = zeros(K,ntol)
precision = zeros(K,ntol)
timing = zeros(K,ntol)
for k=1:K
  texp = 0.01*10.0^(rand())
  q1,q2 = rand(2)
  u_n=[2sqrt(q1)*q2,q1*(1-2*q2)]
  r=0.01*10.0^rand()
  b=rand()*(1+r)
  trans = transit_init(r,b,u_n,true)
  tau = 1.0 # Duration of transit
  t .= linearspace(-0.5*(tau+texp),0.5*(tau+texp),nt)
  v = 2/tau*sqrt((1+r)^2-b^2)
#  t  .= x/v
  param = [0.0,v,b]   # [t_0,v,b_0]
  t_nobin1 = time_ns()
  compute_lightcurve!(trans,param,t,favg0,nt)
  t_nobin = time_ns()-t_nobin1
  integrate_lightcurve!(trans,param,t,texp,favg2,nt,1e-14,50,neval2,depthmax2)
  for j=1:ntol
    t_bin1 = time_ns()
    integrate_lightcurve!(trans,param,t,texp,favg1,nt,tol[j],50,neval1,depthmax1)
    t_bin = time_ns()-t_bin1
    precision[k,j] += maximum(abs.(favg1-favg2))
    neval_mean[k,j] += mean(neval1) 
    println("k: ",k," j: ",j," tol: ",tol[j]," neval: ",mean(neval1)," relative time: ",t_bin/t_nobin/mean(neval1)," max depth: ",maximum(depthmax1))
  end
#  clf()
#  plot(t,favg2[:,1])
#  read(STDIN,Char)
end
precision_tol = zeros(ntol)
neval_tol = zeros(ntol)
for j=1:ntol
  precision_tol[j]=median(precision[:,j])
  neval_tol[j] = median(neval_mean[:,j])
end
clf()
loglog(precision_tol,neval_tol)
