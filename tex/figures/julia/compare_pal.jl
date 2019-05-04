# Compares with Pal's computation:

if VERSION >= v"0.7"
  using DelimitedFiles
end
using PyPlot
include("../../../src/transit_poly_struct.jl")

data_pal = readdlm("ntiq_test.out")

# Now, run the same with transit_poly:

u = [0.2,0.3]
nb = 1001; b=zeros(nb)
r = 0.1  # radius ratio is 0.1
flux = zeros(nb); fgrad = zeros(4)
flux_big = zeros(nb); fgrad_big = zeros(BigFloat,4)
flux_grad = zeros(nb,4)
flux_grad_big = zeros(nb,4)
t1 = time_ns()
for i=0:1000
  b[i+1] = sqrt(((i-500.0)/500.0*(1.0+2.0*r))^2)
  flux[i+1]= transit_poly!(r,b[i+1],u,fgrad)
  flux_grad[i+1,:] = fgrad
end
println("elapsed time: " * string(1e-9 * (time_ns() - t1)) * " seconds")
for i=0:1000
  flux_big[i+1]= convert(Float64,transit_poly!(big(r),big(b[i+1]),big.(u),fgrad_big))
  flux_grad_big[i+1,:] = convert(Array{Float64,1},fgrad_big)
end

# Now, plot results:

fig,axes = subplots(3,1, figsize=(8,12))

ax = axes[1]
ax[:plot](b,flux,label="ALFM (2019)",linestyle="-",lw=2,color="C0")
ax[:plot](b,1 .- data_pal[:,1],label="Pal (2008)",linestyle="--",lw=1,color="C1")
ax[:set_ylabel]("Flux", fontsize=14)
ax[:legend](loc="upper left", fontsize=8)
ax[:set_title](L"$r = 0.1$", fontsize=14)

ax = axes[2]
ax[:plot](b,flux_grad[:,1],linestyle="-",lw=2,label="ALFM (2019)", color="C0")
ax[:plot](b,flux_grad[:,2],linestyle="-",lw=2, color="C0")
ax[:plot](b,50.0*flux_grad[:,3],linestyle="-",lw=2, color="C0")
ax[:plot](b,50.0*flux_grad[:,4],linestyle="-",lw=2, color="C0")
ax[:plot](b,-data_pal[:,2],ls="--",lw=1, label="Pal (2008)", color="C1")
ax[:plot](b,-data_pal[:,3],ls="--",lw=1, color="C1")
ax[:plot](b,-50.0*data_pal[:,4],ls="--",lw=1, color="C1")
ax[:plot](b,-50.0*data_pal[:,5],ls="--",lw=1, color="C1")
ax[:legend](loc="upper left",fontsize=8)
ax[:set_ylabel]("Derivatives", fontsize=14)

ax[:annotate](L"$\partial f / \partial r$", xy=(0.7, -0.185), xycoords="data",
              xytext=(0, 0), textcoords="offset points",
              ha="center", va="center", rotation=11, fontsize=10)

ax[:annotate](L"$\partial f / \partial b$", xy=(0.1, 0.025), xycoords="data",
              xytext=(0, 0), textcoords="offset points",
              ha="center", va="center", rotation=0)

ax[:annotate](L"$50\times\partial f / \partial u_1$", xy=(0.1, -0.18), xycoords="data",
              xytext=(0, 0), textcoords="offset points",
              ha="center", va="center", rotation=4)

ax[:annotate](L"$50\times\partial f / \partial u_2$", xy=(0.1, -0.085), xycoords="data",
              xytext=(0, 0), textcoords="offset points",
              ha="center", va="center", rotation=0)


ax = axes[3]
ax[:semilogy](b,abs.(1 .- data_pal[:,1]-flux_big),label=L"$f$")
ax[:plot](b,abs.(-data_pal[:,2]-flux_grad_big[:,1]),label=L"$\frac{\partial f}{\partial r}$")
ax[:plot](b,abs.(-data_pal[:,3]-flux_grad_big[:,2]),label=L"$\frac{\partial f}{\partial b}$")
ax[:plot](b,abs.(-data_pal[:,4]-flux_grad_big[:,3]),label=L"$\frac{\partial f}{\partial u_1}$")
ax[:plot](b,abs.(-data_pal[:,5]-flux_grad_big[:,4]),label=L"$\frac{\partial f}{\partial u_2}$")
ax[:set_xlabel](L"$b$", fontsize=18)
ax[:set_ylabel]("Difference", fontsize=14)
ax[:legend](loc="upper right",fontsize=10)

savefig("compare_pal.pdf", bbox_inches="tight")
