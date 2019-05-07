 # This is test of the code for computing a transit model integrated over a time step2
# Computes derivatives over the timestep.
using PyPlot

#include("../src/integrate_transit_cubature.jl")
include("../src/integrate_transit_simpson.jl")

# Test it out:

include("../test/loglinspace.jl")

t1 = -1.5; t2 = 1.5; nt = 10000; dt = 0.3
t = zeros(nt)
t .= linearspace(t1,t2,nt)
#t = [t[1235]] ; nt =1
# The following compares different tolerances and maxdepths:
r = 0.1; b0 = 0.5; u_n = [0.3,0.3]; nu = length(u_n)
favg0 = zeros(nt,5+nu)
favg1 = zeros(nt,5+nu)
neval1= zeros(Int64,nt)
favg2 = zeros(nt,5+nu)
neval2= zeros(Int64,nt)
favg3 = zeros(nt,5+nu)
neval3= zeros(Int64,nt)
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

function integrate_lightcurve!(trans::Transit_Struct{T},param::Array{T,1},t::Array{T,1},dt::T,favg1::Array{T,2},nt::Int64,tol::T,maxdepth::Int64,neval_t::Array{Int64,1}) where {T <: Real}
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
#      neval_t[i] += integrate_timestep_gradient!(param,trans,tlim[j],tlim[j+1],tol*trans.r^2,maxdepth,fint)
       nevali,depthmaxi =  integrate_timestep_gradient!(param,trans,tlim[j],tlim[j+1],tol*trans.r^2,maxdepth,fint)
       neval_t[i] += nevali
      ftmp += fint
    end
    ftmp *= dtinv
  end
  favg1[i,1:5]=ftmp[1:5]
  # Convert from d_n to u_n derivatives:
  favg1[i,6:5+trans.n]=BLAS.gemv('T',1.0,trans.dgdu,ftmp[6:6+trans.n])
#  println("i: ",i," t: ",t[i]," result: ",ftmp)
end
return
end

integrate_lightcurve!(trans,param,t,dt,favg1,nt,1e-4,3,neval1)
@time integrate_lightcurve!(trans,param,t,dt,favg1,nt,1e-4,8,neval1)
@time integrate_lightcurve!(trans,param,t,dt,favg2,nt,1e-6,128,neval2)
@time integrate_lightcurve!(trans,param,t,dt,favg3,nt,1e-8,128,neval3)

#@time for i=1:nt
#        ftmp = integrate_timestep_gradient(param,trans,t[i],dt,1e-6*r^2,32)*dtinv
#        favg1[i,1:5]=ftmp[1:5]
#        # Convert from d_n to u_n derivatives:
#        favg1[i,6:5+nu]=BLAS.gemv('T',1.0,trans.dgdu,ftmp[6:6+nu])
#      end

# Now plot the results:
fig,axes = subplots(2,3,figsize=(12,8))
ax = axes[1]
ax[:plot](t,favg0[:,1],linestyle="--",linewidth=2,label=L"$\mathcal{F}$")
ax[:plot](t,favg1[:,1],linewidth=2,label=L"$\overline{\mathcal{F}}$")
ax[:plot]([-dt/2,-dt/2,-dt/2,dt/2,dt/2,dt/2],[0.9925,0.9935,0.993,0.993,0.9925,0.9935],color="b")
ax[:annotate](L"$\Delta t$",xy=[-0.08,0.991])
ax[:legend](loc="upper center")
ax[:axis]([-1.5,1.5,0.988,1.002])
ax[:set_xlabel]("Time")
ax[:set_ylabel]("Flux")
ax = axes[2]
ax[:plot](t,favg0[:,2]*100.,linestyle="--",linewidth=2,label=L"$\frac{\partial\mathcal{F}}{\partial r}$")
ax[:plot](t,favg1[:,2]*100.,linewidth=2,label=L"$\frac{\partial\overline{\mathcal{F}}}{\partial r}$")
ax[:legend](loc="upper center")
ax[:axis]([-1.5,1.5,-30,2])
ax[:set_xlabel]("Time")
ax[:set_ylabel]("Flux derivative [pct]")
ax = axes[3]
ax[:plot](t,favg0[:,3]*100.,linestyle="--",linewidth=2,label=L"$\frac{\partial\mathcal{F}}{\partial t_0}$")
ax[:plot](t,favg1[:,3]*100.,linewidth=2,label=L"$\frac{\partial\overline{\mathcal{F}}}{\partial t_0}$")
ax[:set_xlabel]("Time")
ax[:set_ylabel]("Flux derivative [pct]")
#ax[:axis]([-1.5,1.5,-1,6])
ax[:legend](loc="lower left")
ax = axes[4]
ax[:plot](t,favg0[:,4]*100.,linestyle="--",linewidth=2,label=L"$\frac{\partial\mathcal{F}}{\partial v}$")
ax[:plot](t,favg1[:,4]*100.,linewidth=2,label=L"$\frac{\partial\overline{\mathcal{F}}}{\partial v}$")
ax[:set_xlabel]("Time")
ax[:set_ylabel]("Flux derivative [pct]")
ax[:axis]([-1.5,1.5,-1,6])
ax[:legend](loc="upper center")
ax = axes[5]
ax[:plot](t,favg0[:,5]*100.,linestyle="--",linewidth=2,label=L"$\frac{\partial\mathcal{F}}{\partial b_0}$")
ax[:plot](t,favg1[:,5]*100.,linewidth=2,label=L"$\frac{\partial\overline{\mathcal{F}}}{\partial b_0}$")
ax[:set_xlabel]("Time")
ax[:set_ylabel]("Flux derivative [pct]")
ax[:legend](loc="upper center")
ax = axes[6]
ax[:plot](t,favg0[:,6]*1000.,linestyle="--",linewidth=2,label=L"$\frac{\partial\mathcal{F}}{\partial u_1}$")
ax[:plot](t,favg1[:,6]*1000.,linewidth=2,label=L"$\frac{\partial\overline{\mathcal{F}}}{\partial u_1}$")
ax[:plot](t,favg0[:,7]*1000.,linestyle="--",linewidth=2,label=L"$\frac{\partial\mathcal{F}}{\partial u_2}$")
ax[:plot](t,favg1[:,7]*1000.,linewidth=2,label=L"$\frac{\partial\overline{\mathcal{F}}}{\partial u_2}$")
ax[:set_xlabel]("Time")
ax[:set_ylabel]("Flux derivative [ppt]")
ax[:legend](loc="upper center")
#b = sqrt.(param[3]^2+(param[2]*(t-param[1])).^2)
#plot(dt^2*(d2fdb2.*(param[2]*(t-param[1])./b).^2+dfdb.*param[2]^2./b)/24.0)
savefig("integrate_transit_gradient.pdf", bbox_inches="tight")

read(STDIN,Char)

for i=1:7
#  clf(); plot(t,favg1[:,i]); plot(t,favg2[:,i]); plot(t,favg3[:,i]); read(STDIN,Char); clf(); plot(t,abs.(favg1[:,i]-favg3[:,i]));semilogy(t,abs.(favg2[:,i]-favg3[:,i])); read(STDIN,Char)
  clf(); plot(t,favg1[:,i]); plot(t,favg2[:,i]); plot(t,favg3[:,i]); read(STDIN,Char); clf(); plot(t,favg1[:,i]-favg3[:,i]);plot(t,favg2[:,i]-favg3[:,i]); read(STDIN,Char)
end


