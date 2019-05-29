using PyPlot
using BenchmarkTools
using PyCall
exo = pyimport("exoplanet")
include("../../../src/integrate_lightcurve.jl")
include("../../../test/loglinspace.jl")
fig, axes = subplots(3,1, figsize=(8,14))
p1 = axes[1]

period_vals = collect(logarithmspace(0.4, 3.0, 25)) # orbital period (~2.5 to 10000 days)
r = 0.1 # radius of planet in units of stellar radius
b = 0.3 # impact parameter

a = zeros(length(period_vals))
aonr = zeros(length(period_vals))

for i = 1:length(period_vals)

    orbit = exo.orbits.KeplerianOrbit(period=period_vals[i], b=0.3)
    ai = orbit.a.eval()
    t0 = orbit.t0.eval()
    #transit = orbit.in_transit
    println("$ai")
    #a[i] = ai
    #aonr[i] = a[i] / r

end

#p1.plot(period_vals, a, label = "a", linestyle = "-", lw=2, color=:green)
#p1.plot(period_vals, aonr, label = "aonr", linestyle = "-", lw=2, color=:blue)
#p1.set_xscale(:log)
#p1.set_yscale(:log)
p1.legend(loc="upper left",fontsize=10)
