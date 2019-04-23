# This is test of the code for computing a transit model integrated over a time step.
# Computes derivatives over the timestep.
using PyPlot

if VERSION >= v"0.7"
  using Statistics
end

#include("../src/integrate_transit_cubature.jl")
#include("../src/integrate_transit_simpson_vec.jl")
include("../src/integrate_lightcurve.jl")
#include("integrate_transit_simpson_vec.jl")

# Test it out:

#include("../../test/loglinspace.jl")
include("loglinspace.jl")

t1 = -1.5; t2 = 1.5; nt = 10000; dt = 0.3
t = zeros(nt)
t .= linearspace(t1,t2,nt)
#t = [t[1235]] ; nt =1
# The following compares different tolerances and maxdepths:
r = 0.1; b0 = 0.5; u_n = [0.3,0.3]; nu = length(u_n)
favg0 = zeros(5+nu,nt)
neval0= ones(Int64,nt)
favg1 = zeros(5+nu,nt)
neval1= zeros(Int64,nt)
depthmax1= zeros(Int64,nt)
favg2 = zeros(5+nu,nt)
neval2= zeros(Int64,nt)
depthmax2= zeros(Int64,nt)
favg3 = zeros(5+nu,nt)
neval3= zeros(Int64,nt)
depthmax3= zeros(Int64,nt)
trans = transit_init(r,b0,u_n,true)
param = [0.0,1.0,b0]   # [t_0,v,b_0]

# Now, carry out speed test:

function compute_lightcurve!(trans::Transit_Struct{T},param::Array{T,1},t::Array{T,1},favg0::Array{T,2},nt::Int64,neval::Array{Int64,1}) where {T <: Real}
@inbounds for i=1:nt
  if neval[i] == 0
    neval[i] = 1
  end
  for j=1:neval[i]
    trans.b = sqrt(param[3]^2+(param[2]*(t[i]-param[1]))^2)
    favg0[1,i] =  transit_poly!(trans)
    binv = inv(trans.b)
    favg0[2,i] =  trans.dfdrb[1]
    favg0[3,i] =  trans.dfdrb[2]*binv*param[2]^2*(param[1]-t[i])
    favg0[4,i] =  trans.dfdrb[2]*param[2]*binv*(t[i]-param[1])^2
    favg0[5,i] =  trans.dfdrb[2]*param[3]*binv
    favg0[6:5+trans.n,i] = trans.dfdu
  end
end
return
end

compute_lightcurve!(trans,param,t,favg0,nt,neval0)
@time compute_lightcurve!(trans,param,t,favg0,nt,neval0)


#integrate_lightcurve!(trans,param,t,dt,favg1,nt,1e-4,3,neval1)
#@time integrate_lightcurve!(trans,param,t,dt,favg1,nt,1e-4,8,neval1)

ntol = 9
tol = logarithmspace(-10.0,-2.0,ntol)
K = 1
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
#  t_nobin1 = time_ns()
#  compute_lightcurve!(trans,param,t,favg0,nt,neval0)
#  t_nobin = time_ns()-t_nobin1
  integrate_lightcurve!(trans,param,t,texp,favg2,nt,1e-12,30,neval2,depthmax2)
  println("r: ",r," b: ",b," texp: ",texp)
  for j=1:ntol
    t_bin1 = time_ns()
    integrate_lightcurve!(trans,param,t,texp,favg1,nt,tol[j],30,neval1,depthmax1)
    t_bin = time_ns()-t_bin1
    t_nobin1 = time_ns()
    compute_lightcurve!(trans,param,t,favg0,nt,neval1)
    t_nobin = time_ns()-t_nobin1
    # Now time an "unbinned" version:
    precision[k,j] += maximum(abs.(favg1-favg2))
    neval_mean[k,j] += mean(neval1) 
    println("k: ",k," j: ",j," tol: ",tol[j]," neval: ",mean(neval1)," relative time: ",t_bin/t_nobin," max depth: ",maximum(depthmax1))
  end
  clf()
  plot(t,favg2[1,:])
  plot(reverse(t),favg2[1,:])
  if VERSION >= v"0.7"
    read(stdin,Char)
  else
    read(STDIN,Char)
  end
end
precision_tol = zeros(ntol)
neval_tol = zeros(ntol)
for j=1:ntol
  precision_tol[j]=median(precision[:,j])
  neval_tol[j] = median(neval_mean[:,j])
end
clf()
loglog(precision_tol,neval_tol)
