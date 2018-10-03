# Defines the function that computes the transit:
include("../../../src/transit_poly_struct.jl")
using PyPlot
using PyCall


# System properties
npts = 1000
iu = 3
u_n = 0.5*ones(iu) / iu
r = 0.1
b = zeros(npts)
for i=1:npts
    b[i] = 3.0 * (i / float(npts) - 0.5)
end

# Julia flux
flux_julia = zeros(npts)
trans = transit_init(r, 0.0, u_n, true)
for i=1:npts
    trans.b = abs(b[i])
    flux_julia[i] = transit_poly!(trans)
end

# Starry flux
@pyimport starry
map = starry.Map(iu)
for i=1:iu
    map[:__setitem__](i, u_n[i])
end
flux_starry = map[:flux](xo=b, ro=r)

# Plot
plot(b, flux_julia, linewidth=2, label="julia")
plot(b, flux_starry, linewidth=1, label="starry")
xlabel("Impact parameter")
ylabel("Flux")
legend()
ylim(0.985, 1.01)
savefig("compare_to_starry_" * string(iu) * ".pdf", bbox_inches="tight")
