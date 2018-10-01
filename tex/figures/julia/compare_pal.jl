# Compares with Pal's computation:

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
tic()
for i=0:1000
  b[i+1] = sqrt(((i-500.0)/500.0*(1.0+2.0*r))^2)
  flux[i+1]= transit_poly!(r,b[i+1],u,fgrad)
  flux_grad[i+1,:] = fgrad
end
toc()
for i=0:1000
  flux_big[i+1]= convert(Float64,transit_poly!(big(r),big(b[i+1]),big.(u),fgrad_big))
  flux_grad_big[i+1,:] = convert(Array{Float64,1},fgrad_big)
end

# Now, plot results:

fig,axes = subplots(3,1)

ax = axes[1]
ax[:plot](b,1.0-data_pal[:,1],label="Pal",lw=3)
ax[:plot](b,flux,label="AL 2018",linestyle="-.",lw=3,color="r")
ax[:set_xlabel]("b")
ax[:set_ylabel]("Flux")
ax[:legend](loc="upper left")
ax[:set_title]("r = 0.1")

ax = axes[2]
ax[:plot](b,-data_pal[:,2],label="df/dr Pal",lw=3)
ax[:plot](b,-data_pal[:,3],label="df/db Pal",lw=3)
ax[:plot](b,-50.0*data_pal[:,4],label="50x df/du1 Pal",lw=3)
ax[:plot](b,-50.0*data_pal[:,5],label="50x df/du2 Pal",lw=3)
ax[:plot](b,flux_grad[:,1],label="df/dr AL 2018",linestyle="-.",lw=3)
ax[:plot](b,flux_grad[:,2],label="df/db AL 2018",linestyle="-.",lw=3)
ax[:plot](b,50.0*flux_grad[:,3],label="50x df/du1 AL 2018",linestyle="-.",lw=3)
ax[:plot](b,50.0*flux_grad[:,4],label="50x df/du2 AL 2018",linestyle="-.",lw=3)
ax[:legend](loc="center left",fontsize=6)
ax[:set_xlabel]("b")
ax[:set_ylabel]("Derivatives")


ax = axes[3]
ax[:semilogy](b,abs.(1.-data_pal[:,1]-flux_big),label="diff flux Pal-AL 2018 BigFloat")
ax[:plot](b,abs.(-data_pal[:,2]-flux_grad_big[:,1]),label="diff df/dr Pal-AL 2018 BigFloat")
ax[:plot](b,abs.(-data_pal[:,3]-flux_grad_big[:,2]),label="diff df/db Pal-AL 2018 BigFloat")
ax[:plot](b,abs.(-data_pal[:,4]-flux_grad_big[:,3]),label="diff df/du1 Pal-AL 2018 BigFloat")
ax[:plot](b,abs.(-data_pal[:,5]-flux_grad_big[:,4]),label="diff df/du2 Pal-AL 2018 BigFloat")
ax[:set_xlabel]("b")
ax[:set_ylabel]("Differences of flux, derivatives")
ax[:legend](loc="lower left",fontsize=6)

savefig("compare_pal_r0pt1.pdf", bbox_inches="tight")
