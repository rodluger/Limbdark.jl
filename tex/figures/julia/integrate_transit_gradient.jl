
# Computes derivatives over the timestep.
using PyPlot

include("../../../src/integrate_lightcurve.jl")

# Test it out:

include("../../../test/loglinspace.jl")

t1 = -1.5; t2 = 1.5; nt = 10000; dt = 0.3
t = zeros(nt)
t .= linearspace(t1,t2,nt)
# The following compares different tolerances and maxdepths:
r = 0.1; b0 = 0.5; u_n = [0.3,0.3]; nu = length(u_n)
favg0 = zeros(5+nu,nt)
favg1 = zeros(5+nu,nt)
neval1= zeros(Int64,nt)
depthmax1= zeros(Int64,nt)
trans = transit_init(r,b0,u_n,true)
param = [0.0,1.0,b0]   # [t_0,v,b_0]

# Now, carry out speed test:

function compute_lightcurve!(trans::Transit_Struct{T},param::Array{T,1},t::Array{T,1},favg0::Array{T,2},nt::Int64) where {T <: Real}
@inbounds for i=1:nt
  trans.b = sqrt(param[3]^2+(param[2]*(t[i]-param[1]))^2)
  favg0[1,i] =  transit_poly!(trans)-1
  favg0[2,i] =  trans.dfdrb[1]
  favg0[3,i] =  trans.dfdrb[2]/trans.b*param[2]^2*(param[1]-t[i])
  favg0[4,i] =  trans.dfdrb[2]*param[2]/trans.b*(t[i]-param[1])^2
  favg0[5,i] =  trans.dfdrb[2]*param[3]/trans.b
  favg0[6:5+trans.n,i] = trans.dfdu
end
return
end

compute_lightcurve!(trans,param,t,favg0,nt)
@time compute_lightcurve!(trans,param,t,favg0,nt)

integrate_lightcurve!(trans,param,t,dt,favg1,nt,1e-10,30,neval1,depthmax1)

# Now plot the results:
fig,axes = subplots(2,3,figsize=(12,8))
if VERSION >= v"0.7"
  ax = axes[1]
  ax.plot(t,favg0[1,:].+1,linestyle="--",linewidth=2,label=L"$\mathcal{F}$")
  ax.plot(t,favg1[1,:].+1,linewidth=2,label=L"$\overline{\mathcal{F}}$")
  ax.plot([-dt/2,-dt/2,-dt/2,dt/2,dt/2,dt/2],[0.9925,0.9935,0.993,0.993,0.9925,0.9935],color="b")
  ax.annotate(L"$\Delta t$",xy=[-0.08,0.991])
  ax.legend(loc="upper center")
  ax.axis([-1.5,1.5,0.988,1.002])
  ax.set_xlabel("Time")
  ax.set_ylabel("Flux")
  ax = axes[2]
  ax.plot(t,favg0[2,:]*100.,linestyle="--",linewidth=2,label=L"$\frac{\partial\mathcal{F}}{\partial r}$")
  ax.plot(t,favg1[2,:]*100.,linewidth=2,label=L"$\frac{\partial\overline{\mathcal{F}}}{\partial r}$")
  ax.legend(loc="upper center")
  ax.axis([-1.5,1.5,-30,2])
  ax.set_xlabel("Time")
  ax.set_ylabel("Flux derivative [pct]")
  ax = axes[3]
  ax.plot(t,favg0[3,:]*100.,linestyle="--",linewidth=2,label=L"$\frac{\partial\mathcal{F}}{\partial t_0}$")
  ax.plot(t,favg1[3,:]*100.,linewidth=2,label=L"$\frac{\partial\overline{\mathcal{F}}}{\partial t_0}$")
  ax.set_xlabel("Time")
  ax.set_ylabel("Flux derivative [pct]")
  ax.legend(loc="lower left")
  ax = axes[4]
  ax.plot(t,favg0[4,:]*100.,linestyle="--",linewidth=2,label=L"$\frac{\partial\mathcal{F}}{\partial v}$")
  ax.plot(t,favg1[4,:]*100.,linewidth=2,label=L"$\frac{\partial\overline{\mathcal{F}}}{\partial v}$")
  ax.set_xlabel("Time")
  ax.set_ylabel("Flux derivative [pct]")
  ax.axis([-1.5,1.5,-1,6])
  ax.legend(loc="upper center")
  ax = axes[5]
  ax.plot(t,favg0[5,:]*100.,linestyle="--",linewidth=2,label=L"$\frac{\partial\mathcal{F}}{\partial b_0}$")
  ax.plot(t,favg1[5,:]*100.,linewidth=2,label=L"$\frac{\partial\overline{\mathcal{F}}}{\partial b_0}$")
  ax.set_xlabel("Time")
  ax.set_ylabel("Flux derivative [pct]")
  ax.legend(loc="upper center")
  ax = axes[6]
  ax.plot(t,favg0[6,:]*1000.,linestyle="--",linewidth=2,label=L"$\frac{\partial\mathcal{F}}{\partial u_1}$")
  ax.plot(t,favg1[6,:]*1000.,linewidth=2,label=L"$\frac{\partial\overline{\mathcal{F}}}{\partial u_1}$")
  ax.plot(t,favg0[7,:]*1000.,linestyle="--",linewidth=2,label=L"$\frac{\partial\mathcal{F}}{\partial u_2}$")
  ax.plot(t,favg1[7,:]*1000.,linewidth=2,label=L"$\frac{\partial\overline{\mathcal{F}}}{\partial u_2}$")
  ax.set_xlabel("Time")
  ax.set_ylabel("Flux derivative [ppt]")
  ax.legend(loc="upper center")
  #b = sqrt.(param[3]^2+(param[2]*(t-param[1])).^2)
  #plot(dt^2*(d2fdb2.*(param[2]*(t-param[1])./b).^2+dfdb.*param[2]^2./b)/24.0)
else
  # Now plot the results:
  fig,axes = subplots(2,3,figsize=(12,8))
  ax = axes[1]
  ax[:plot](t,favg0[1,:].+1,linestyle="--",linewidth=2,label=L"$\mathcal{F}$")
  ax[:plot](t,favg1[1,:].+1,linewidth=2,label=L"$\overline{\mathcal{F}}$")
  ax[:plot]([-dt/2,-dt/2,-dt/2,dt/2,dt/2,dt/2],[0.9925,0.9935,0.993,0.993,0.9925,0.9935],color="b")
  ax[:annotate](L"$\Delta t$",xy=[-0.08,0.991])
  ax[:legend](loc="upper center")
  ax[:axis]([-1.5,1.5,0.988,1.002])
  ax[:set_xlabel]("Time")
  ax[:set_ylabel]("Flux")
  ax = axes[2]
  ax[:plot](t,favg0[2,:]*100.,linestyle="--",linewidth=2,label=L"$\frac{\partial\mathcal{F}}{\partial r}$")
  ax[:plot](t,favg1[2,:]*100.,linewidth=2,label=L"$\frac{\partial\overline{\mathcal{F}}}{\partial r}$")
  ax[:legend](loc="upper center")
  ax[:axis]([-1.5,1.5,-30,2])
  ax[:set_xlabel]("Time")
  ax[:set_ylabel]("Flux derivative [pct]")
  ax = axes[3]
  ax[:plot](t,favg0[3,:]*100.,linestyle="--",linewidth=2,label=L"$\frac{\partial\mathcal{F}}{\partial t_0}$")
  ax[:plot](t,favg1[3,:]*100.,linewidth=2,label=L"$\frac{\partial\overline{\mathcal{F}}}{\partial t_0}$")
  ax[:set_xlabel]("Time")
  ax[:set_ylabel]("Flux derivative [pct]")
  ax[:legend](loc="lower left")
  ax = axes[4]
  ax[:plot](t,favg0[4,:]*100.,linestyle="--",linewidth=2,label=L"$\frac{\partial\mathcal{F}}{\partial v}$")
  ax[:plot](t,favg1[4,:]*100.,linewidth=2,label=L"$\frac{\partial\overline{\mathcal{F}}}{\partial v}$")
  ax[:set_xlabel]("Time")
  ax[:set_ylabel]("Flux derivative [pct]")
  ax[:axis]([-1.5,1.5,-1,6])
  ax[:legend](loc="upper center")
  ax = axes[5]
  ax[:plot](t,favg0[5,:]*100.,linestyle="--",linewidth=2,label=L"$\frac{\partial\mathcal{F}}{\partial b_0}$")
  ax[:plot](t,favg1[5,:]*100.,linewidth=2,label=L"$\frac{\partial\overline{\mathcal{F}}}{\partial b_0}$")
  ax[:set_xlabel]("Time")
  ax[:set_ylabel]("Flux derivative [pct]")
  ax[:legend](loc="upper center")
  ax = axes[6]
  ax[:plot](t,favg0[6,:]*1000.,linestyle="--",linewidth=2,label=L"$\frac{\partial\mathcal{F}}{\partial u_1}$")
  ax[:plot](t,favg1[6,:]*1000.,linewidth=2,label=L"$\frac{\partial\overline{\mathcal{F}}}{\partial u_1}$")
  ax[:plot](t,favg0[7,:]*1000.,linestyle="--",linewidth=2,label=L"$\frac{\partial\mathcal{F}}{\partial u_2}$")
  ax[:plot](t,favg1[7,:]*1000.,linewidth=2,label=L"$\frac{\partial\overline{\mathcal{F}}}{\partial u_2}$")
  ax[:set_xlabel]("Time")
  ax[:set_ylabel]("Flux derivative [ppt]")
  ax[:legend](loc="upper center")
end
savefig("integrate_transit_gradient.pdf", bbox_inches="tight")
