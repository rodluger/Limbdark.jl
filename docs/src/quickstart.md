# Quick start

## Compute a limb-darkened transit light curve

```@example 1
using Plots # hide
using Plots.PlotMeasures # hide
push!(LOAD_PATH, "../../src/") # hide
import Limbdark

# Array of impact parameters
npts = 100
b = zeros(npts)
for i = 1:npts
    b[i] = 3.0 * ((i - 0.5) / float(npts) - 0.5)
end

# 5th order limb darkening coefficients
nu = 5
u_n = 0.5 * ones(nu) / nu

# Occultor radius
r = 0.1

# Initialize the transit struct
trans = Limbdark.transit_init(r, 0.0, u_n, false)

# Compute the flux
flux = zeros(npts)
for i = 1:npts
    trans.b = abs(b[i])
    flux[i] = Limbdark.transit_poly!(trans)
end

# Plot it
plot(b, flux, linewidth=2, 
     xlabel="Impact parameter", 
     ylabel="Flux",
     label="light curve");
savefig("example1.svg"); nothing # hide
```

![](example1.svg)