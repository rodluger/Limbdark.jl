# Defines the function that computes the transit:
include("../../../src/transit_poly_struct.jl")
using PyPlot
using PyCall


# System properties
npts = 1000
iu = 30
u_n = 0.5*ones(iu) / iu
r = 0.1
b = zeros(npts)
for i=1:npts
    b[i] = 3.0 * (i / float(npts) - 0.5)
end

# Julia flux
flux_julia = zeros(npts)
flux_big = zeros(npts)
trans = transit_init(r, 0.0, u_n, false)
trans_big = transit_init(big(r), big(0.0), big.(u_n), false)
for i=1:npts
    trans.b = abs(b[i])
#    trans = transit_init(r, abs(b[i]), u_n, true)
    flux_julia[i] = transit_poly!(trans)
    trans_big.b = big(abs(b[i]))
    flux_big[i] = convert(Float64,transit_poly!(trans_big))
#    flux_julia[i] = transit_poly(r,abs(b[i]),u_n)
end

# Starry flux
@pyimport starry
map = starry.Map(iu)
for i=1:iu
    map[:__setitem__](i, u_n[i])
end
flux_starry = map[:flux](xo=b, ro=r)

# Plot
fig,axes = subplots(2,1)
ax = axes[1]
ax[:plot](b, flux_julia, linewidth=2, label="julia")
ax[:plot](b, flux_starry, linewidth=1, label="starry")
ax[:set_xlabel]("Impact parameter")
ax[:set_ylabel]("Flux")
ax[:legend]()
ax[:axis]([-1.5,1.5,0.985, 1.01])

ax = axes[2]
ax[:semilogy](b, abs.(flux_julia-flux_big), linewidth=2, label="julia-big")
ax[:plot](b, abs.(flux_julia-flux_starry), linewidth=2, label="julia-starry")
ax[:plot](b, abs.(flux_starry-flux_big), linewidth=2, label="starry-big")
ax[:legend]()
ax[:axis]([-1.5,1.5,1e-10,1e-5])

savefig("compare_to_starry_" * string(iu) * ".pdf", bbox_inches="tight")


