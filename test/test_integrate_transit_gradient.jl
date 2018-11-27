# This is test of the code for computing a transit model integrated over a time step.
# Computes derivatives over the timestep.
using PyPlot

include("../src/integrate_transit.jl")

# Test it out:

include("../test/loglinspace.jl")

t1 = -1.2; t2 = 1.2; nt = 1000; dt = (t2-t1)/nt*100
t = linearspace(t1,t2,nt)
# The following compares different tolerances and maxdepths:
r = 0.1; b0 = 0.5; u_n = [0.5,0.5]; nu = length(u_n)
favg0 = zeros(nt,3+nu)
favg1 = zeros(nt,3+nu)
trans = transit_init(r,b0,u_n,true)
param = [0.0,1.0,r,b0]

# Now, carry out speed test:

@time for i=1:nt
        trans.b = sqrt(param[4]^2+(param[2]*(t[i]-param[1]))^2)
        favg0[i,1] = transit_poly!(trans)
        favg0[i,2:3] = trans.dfdrb
        favg0[i,4:3+nu] = trans.dfdu
      end
dtinv = inv(dt)
@time for i=1:nt
        ftmp = integrate_timestep_gradient(param,trans,t[i],dt,1e-5*r^2,8)*dtinv
        favg1[i,1:3]=ftmp[1:3]
        # Convert from d_n to u_n derivatives:
        favg1[i,4:3+nu]=BLAS.gemv('T',1.0,trans.dddu,ftmp[4:4+nu])
      end

# Now plot the results:
fig,axes = subplots(2,2)
ax = axes[1]
ax[:plot](t,favg0[:,1],linewidth=2,label=L"$\mathcal{F}$")
ax[:plot](t,favg1[:,1],linestyle="--",linewidth=2,label=L"$\overline\mathcal{F}$")
ax[:plot]([-dt/2,-dt/2,-dt/2,dt/2,dt/2,dt/2],[0.9925,0.9935,0.993,0.993,0.9925,0.9935],color="b")
ax[:annotate](L"$\Delta t$",xy=[-0.08,0.991])
ax[:legend](loc="upper center")
ax[:set_xlabel]("Time")
ax[:set_ylabel]("Flux")
ax = axes[2]
ax[:plot](t,favg0[:,2],linewidth=2,label=L"$\frac{\partial\mathcal{F}}{\partial r}$")
ax[:plot](t,favg1[:,2],linestyle="--",linewidth=2,label=L"$\frac{\partial\overline\mathcal{F}}{\partial r}$")
ax[:legend](loc="upper center")
ax[:set_xlabel]("Time")
ax[:set_ylabel]("Flux derivative")
ax = axes[3]
ax[:plot](t,favg0[:,3],linewidth=2,label=L"$\frac{\partial\mathcal{F}}{\partial b}$")
ax[:plot](t,favg1[:,3],linestyle="--",linewidth=2,label=L"$\frac{\partial\overline\mathcal{F}}{\partial b}$")
ax[:set_xlabel]("Time")
ax[:set_ylabel]("Flux derivative")
ax[:legend](loc="upper center")
ax = axes[4]
ax[:plot](t,favg0[:,4],linewidth=2,label=L"$\frac{\partial\mathcal{F}}{\partial u_1}$")
ax[:plot](t,favg1[:,4],linestyle="--",linewidth=2,label=L"$\frac{\partial\overline\mathcal{F}}{\partial u_1}$")
ax[:plot](t,favg0[:,5],linewidth=2,label=L"$\frac{\partial\mathcal{F}}{\partial u_2}$")
ax[:plot](t,favg1[:,5],linestyle="--",linewidth=2,label=L"$\frac{\partial\overline\mathcal{F}}{\partial u_2}$")
ax[:set_xlabel]("Time")
ax[:set_ylabel]("Flux derivative")
ax[:legend](loc="upper center",fontsize=10)
#b = sqrt.(param[4]^2+(param[2]*(t-param[1])).^2)
#plot(dt^2*(d2fdb2.*(param[2]*(t-param[1])./b).^2+dfdb.*param[2]^2./b)/24.0)
