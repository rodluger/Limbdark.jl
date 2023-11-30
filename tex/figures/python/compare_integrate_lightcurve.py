"""Time integration comparison tests"""
# TODO - this is a work in progress. It's purpose is to see if calling python
# from Julia has overhead that isn't present when running that code directly
import exoplanet as exo
import matplotlib.pyplot as pl
import numpy as np

# Input params
u_1 = 0.2
u_2 = 0.2
rplanet = 0.1   # fraction of stellar radius
b0 = 0.2        # impact parameter
P = 1000          # orbital period in days
texp = 2.0 / 60.0

# Timing params
#nN = 8
#Nmax = 5
#Narr = np.array(np.logspace(1, Nmax, nN), dtype=int)
runs = 1

num_datapoints = np.zeros(runs)
limbdark_time = np.zeros(runs)
exoplanet_time = np.zeros(runs)

# Plot
fig = pl.figure(figsize=(7, 4))
ax = pl.subplot2grid((2, 5), (0, 0), colspan=4, rowspan=2)
axleg1 = pl.subplot2grid((2, 5), (0, 4))
axleg1.axis('off')

# Loop over number of cadences
for i in range(runs):

    dataPoints = 1000
    t = np.linspace(-1.5, 1.5, dataPoints)

    orbit = exo.orbits.KeplerianOrbit(period=P, b=b0, ecc=0.0, omega=0.0)

    lc = exo.StarryLightCurve([u_1, u_2]).get_light_curve(orbit=orbit, r=rplanet, t=t, texp=texp, order=2).eval()

    ax.plot(t, lc, label="exoplanet (n=$dataPoints)", linestyle = "-", lw=1)

# TODO get plotting working
pl.show(fig)
