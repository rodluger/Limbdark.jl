
# Computes derivatives over the timestep.
using PyPlot
using Limbdark
import Limbdark: Transit_Struct

#include("../../../src/integrate_lightcurve.jl")

# Test it out:

include("../../../test/loglinspace.jl")

t1 = -1.5; t2 = 1.5; nt = 10000; dt = 0.3
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
favg4 = zeros(nt,5+nu)
favg5 = zeros(nt,5+nu)
favg6 = zeros(nt,5+nu)
favg7 = zeros(nt,5+nu)
neval3= zeros(Int64,nt)
neval4= zeros(Int64,nt)
neval5= zeros(Int64,nt)
neval6= zeros(Int64,nt)
neval7= zeros(Int64,nt)
depthmax3= zeros(Int64,nt)
depthmax4= zeros(Int64,nt)
depthmax5= zeros(Int64,nt)
depthmax6= zeros(Int64,nt)
depthmax7= zeros(Int64,nt)
trans = transit_init(r,b0,u_n,true)
param = [0.0,1.0,b0]   # [t_0,v,b_0]

# Now, carry out speed test:

function compute_lightcurve!(trans::Transit_Struct{T},param::Array{T,1},t::Array{T,1},favg0::Array{T,2},nt::Int64) where {T <: Real}
@inbounds for i=1:nt
  trans.b = sqrt(param[3]^2+(param[2]*(t[i]-param[1]))^2)
  favg0[i,1] =  transit_poly!(trans)-1
  favg0[i,2] =  trans.dfdrb[1]
  favg0[i,3] =  trans.dfdrb[2]/trans.b*param[2]^2*(param[1]-t[i])
  favg0[i,4] =  trans.dfdrb[2]*param[2]/trans.b*(t[i]-param[1])^2
  favg0[i,5] =  trans.dfdrb[2]*param[3]/trans.b
  favg0[i,6:5+trans.n] = trans.dfdu
end
return
end

compute_lightcurve!(trans,param,t,favg0,nt)
timing = zeros(8)
#@time compute_lightcurve!(trans,param,t,favg0,nt)
tstart=time_ns(); compute_lightcurve!(trans,param,t,favg0,nt); timing[1] = (time_ns()-tstart)*1e-9

integrate_lightcurve!(trans,param,t,dt,favg1,nt,1e-4,3,neval1,depthmax1)
#@time integrate_lightcurve!(trans,param,t,dt,favg1,nt,1e-4,8,neval1,depthmax1)
#tic(); integrate_lightcurve!(trans,param,t,dt,favg1,nt,1e-4,8,neval1,depthmax1); timing[2] = toc()
tstart=time_ns(); integrate_lightcurve!(trans,param,t,dt,favg1,nt,1e-4,8,neval1,depthmax1); timing[2] = (time_ns()-tstart)*1e-9
#@time integrate_lightcurve!(trans,param,t,dt,favg2,nt,1e-6,128,neval2,depthmax2)
tstart=time_ns(); integrate_lightcurve!(trans,param,t,dt,favg2,nt,1e-6,128,neval2,depthmax2); timing[3] =  (time_ns()-tstart)*1e-9
#@time integrate_lightcurve!(trans,param,t,dt,favg3,nt,1e-8,128,neval3,depthmax3)
tstart=time_ns(); integrate_lightcurve!(trans,param,t,dt,favg3,nt,1e-8,128,neval3,depthmax3); timing[4] =  (time_ns()-tstart)*1e-9
#@time integrate_lightcurve!(trans,param,t,dt,favg4,nt,1e-10,128,neval4,depthmax4)
tstart=time_ns(); integrate_lightcurve!(trans,param,t,dt,favg4,nt,1e-10,128,neval4,depthmax4); timing[5] =  (time_ns()-tstart)*1e-9
#@time integrate_lightcurve!(trans,param,t,dt,favg5,nt,1e-12,128,neval5,depthmax5)
tstart=time_ns(); integrate_lightcurve!(trans,param,t,dt,favg5,nt,1e-12,128,neval5,depthmax5); timing[6] =  (time_ns()-tstart)*1e-9
#@time integrate_lightcurve!(trans,param,t,dt,favg6,nt,1e-14,128,neval6,depthmax6)
tstart=time_ns(); integrate_lightcurve!(trans,param,t,dt,favg6,nt,1e-14,128,neval6,depthmax6); timing[7] =  (time_ns()-tstart)*1e-9
#@time integrate_lightcurve!(trans,param,t,dt,favg7,nt,1e-16,128,neval7,depthmax7)
tstart=time_ns(); integrate_lightcurve!(trans,param,t,dt,favg7,nt,1e-16,128,neval7,depthmax7); timing[8] =  (time_ns()-tstart)*1e-9

# Now, make some precision plots:

absolute_precision = zeros(7,6)
fractional_precision = zeros(7,6)
depthmax = [maximum(depthmax1),maximum(depthmax2),maximum(depthmax3),maximum(depthmax4),maximum(depthmax5)]
neval = [sum(neval1),sum(neval2),sum(neval3),sum(neval4),sum(neval5),sum(neval6)]
for j=1:7
  absolute_precision[j,1] = maximum(abs.(favg1[:,j]-favg7[:,j]))
  fractional_precision[j,1] = maximum(abs.(favg1[:,j]-favg7[:,j]))/maximum(abs.(favg7[:,j]))
end
for j=1:7
  absolute_precision[j,2] = maximum(abs.(favg2[:,j]-favg7[:,j]))
  fractional_precision[j,2] = maximum(abs.(favg2[:,j]-favg7[:,j]))/maximum(abs.(favg7[:,j]))
end
for j=1:7
  absolute_precision[j,3] = maximum(abs.(favg3[:,j]-favg7[:,j]))
  fractional_precision[j,3] = maximum(abs.(favg3[:,j]-favg7[:,j]))/maximum(abs.(favg7[:,j]))
end
for j=1:7
  absolute_precision[j,4] = maximum(abs.(favg4[:,j]-favg7[:,j]))
  fractional_precision[j,4] = maximum(abs.(favg4[:,j]-favg7[:,j]))/maximum(abs.(favg7[:,j]))
end
for j=1:7
  absolute_precision[j,5] = maximum(abs.(favg5[:,j]-favg7[:,j]))
  fractional_precision[j,5] = maximum(abs.(favg5[:,j]-favg7[:,j]))/maximum(abs.(favg7[:,j]))
end
for j=1:7
  absolute_precision[j,6] = maximum(abs.(favg6[:,j]-favg7[:,j]))
  fractional_precision[j,6] = maximum(abs.(favg6[:,j]-favg7[:,j]))/maximum(abs.(favg7[:,j]))
end

tol = [1e-4,1e-6,1e-8,1e-10,1e-12,1e-14]

clf()
labels = [L"$F-1$",L"$r$",L"$t_0$",L"$v$",L"$b_0$",L"$u_1$",L"$u_2$"]
for j=1:7
#  loglog(tol,absolute_precision[j,:],label=labels[j])
  loglog(neval/nt,absolute_precision[j,:],"o-",label=labels[j])
  #loglog(neval/nt,absolute_precision[j,:],"o",color="k")
#  loglog(tol,fractional_precision[j,:],label=labels[j],linestyle=":")
end
#plot(tol,tol,label="tolerance",linestyle=":")
plot(neval/nt,10*tol,label=L"$10\epsilon_{tol}$",linestyle=":")
legend(loc="lower left")
xlabel("Average number of evaluations per exposure")
ylabel("Absolute precision")
xlim(5, 1200)

twiny()
loglog(timing[2:7]/timing[1],10*tol,label=L"$10\epsilon_{tol}$",linestyle=:dashed)
legend(loc="upper right")
xlim(5, 1200)
xlabel("Time ratio per exposure")

savefig("compare_precision.pdf", bbox_inches="tight")
