# Defines the function that computes the transit:
include("../../../src/transit_poly_struct.jl")
using PyPlot
using PyCall


# System properties
npts = 1000
nu = 30
u_n = 0.5*ones(nu) / nu
r = 0.1
b = zeros(npts)
for i=1:npts
    b[i] = 3.0 * (i / float(npts) - 0.5)
end

# Julia flux
flux_julia = zeros(npts)
#flux_big = zeros(npts)
trans = transit_init(r, 0.0, u_n, false)
dfdx = zeros(nu+2,npts)
#trans_big = transit_init(big(r), big(0.0), big.(u_n), false)
for i=1:npts
    trans.b = abs(b[i])
#    trans = transit_init(r, abs(b[i]), u_n, true)
    flux_julia[i] = transit_poly!(trans)
    dfdx[:,i] = trans.dfdrbu
#    trans_big.b = big(abs(b[i]))
#    flux_big[i] = convert(Float64,transit_poly!(trans_big))
#    flux_julia[i] = transit_poly(r,abs(b[i]),u_n)
end

# Starry flux
@pyimport starry
map = starry.Map(nu)
for i=1:nu
    map[:__setitem__](i, u_n[i])
end
flux_starry, grad = map[:flux](xo=b, ro=r, gradient = true)

# Plot
fig,axes = subplots(2,2)
ax = axes[1]
ax[:plot](b, flux_julia, linewidth=2, label="julia")
ax[:plot](b, flux_starry, linewidth=1, label="starry")
ax[:set_xlabel]("Impact parameter")
ax[:set_ylabel]("Flux")
ax[:legend]()
ax[:axis]([-1.5,1.5,0.985, 1.01])

ax = axes[2]
#ax[:semilogy](b, abs.(flux_julia-flux_big), linewidth=2, label="julia-big")
ax[:semilogy](b, abs.(flux_julia-flux_starry), linewidth=2, label="julia-starry")
#ax[:plot](b, abs.(flux_starry-flux_big), linewidth=2, label="starry-big")
ax[:set_xlabel]("abs(Difference)")
ax[:legend]()
ax[:axis]([-1.5,1.5,1e-10,1e-5])

ax = axes[3]
ax[:plot](b, grad["ro"],label=L"$df/dr$ (starry)")
ax[:plot](b, dfdx[1,:],label=L"$df/dr$ (limbdark)")
ax[:plot](b, grad["xo"],label=L"$df/db$ (starry)")
ax[:plot](b, dfdx[2,:],label=L"$df/db$ (limbdark)")
for iu=1:nu
  ax[:plot](b, grad["u"][iu,:],label=L"$df/du$ (starry)")
  ax[:plot](b, dfdx[2+iu,:],label=L"$df/du$ (limbdark)")
end
ax[:set_xlabel]("Derivative")
ax[:set_ylabel]("Flux")
#ax[:legend](fontsize=6,loc="lower right")

ax = axes[4]
ax[:semilogy](b, abs.(grad["ro"]-dfdx[1,:]),label=L"$df/dr$ diff")
ax[:plot](b, abs.(grad["xo"]-dfdx[2,:]),label=L"$df/db$ diff")
#for iu=1:nu
#  ax[:plot](b, abs.(grad["u"][iu,:]-dfdx[2+iu,:]),label=L"$df/du$ diff")
#end
ax[:set_xlabel]("Derivative differences")
ax[:set_ylabel]("Flux")
ax[:legend](fontsize=1,loc="lower right")
#ax[:axis]([-1.5,1.5,1e-10,1e-5])

savefig("compare_to_starry_" * string(nu) * ".pdf", bbox_inches="tight")
