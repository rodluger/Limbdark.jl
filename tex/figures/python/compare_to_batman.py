"""Starry speed tests."""
from starry.kepler import Primary, Secondary, System
from starry import Map
import time
import matplotlib.pyplot as pl
import numpy as np
import batman
import pytransit
import subprocess

# Marker size is proportional to log error
def ms(error):
    return 18 + np.log10(error)

# Input params
u1 = 0.4
u2 = 0.26
mstar = 1       # solar masses
rstar = 1       # solar radii
rplanet = 0.1   # fraction of stellar radius
b0 = 0.5        # impact parameter
P = 50          # orbital period in days

# Compute the semi-major axis from Kepler's
# third law in units of rstar (for batman)
a = ((P * 86400) ** 2 * (1.32712440018e20 * mstar) /
     (4 * np.pi ** 2)) ** (1. / 3.) / (6.957e8 * rstar)

# Get the inclination in degrees
inc = np.arccos(b0 / a) * 180 / np.pi

# Timing params
number = 10
nN = 8
Nmax = 5
Narr = np.logspace(1, Nmax, nN)
agol_time = np.zeros(nN)
agol_grad_time = np.zeros(nN)
pytransit_time = np.zeros(nN)
batman_time = np.zeros(nN)

# Loop over number of cadences
for i, N in enumerate(Narr):

    # Time array
    t = np.linspace(-0.15, 0.15, N)

    # Use starry to compute b(t)
    star = Primary()
    star[1] = u1
    star[2] = u2
    planet = Secondary()
    planet.r = rplanet
    planet.inc = inc
    planet.porb = P
    planet.a = a
    planet.lambda0 = 90
    system = System(star, planet)
    system.compute(t)
    b = np.sqrt(planet.X ** 2 + planet.Y ** 2)

    # Feed b(t) to julia
    # HACK: PyJulia is currently broken, so this is how we have to do this...
    np.savetxt("b.txt", X=b)
    np.savetxt("u.txt", X=[u1, u2])
    foo = subprocess.check_output(['julia', "compare_to_batman.jl"])
    agol_time[i] = float(foo.decode('utf-8'))
    agol_flux = np.loadtxt("flux.txt")
    foo = subprocess.check_output(['julia', "compare_to_batman_grad.jl"])
    agol_grad_time[i] = float(foo.decode('utf-8'))

    # batman
    params = batman.TransitParams()
    params.limb_dark = "quadratic"
    params.u = [u1, u2]
    params.t0 = 0.
    params.ecc = 0
    params.w = 90.
    params.rp = rplanet
    params.a = a
    params.per = P
    params.inc = inc
    m = batman.TransitModel(params, t, nthreads=1)
    tstart = time.time()
    for k in range(number):
        batman_flux = m.light_curve(params)
    batman_time[i] = (time.time() - tstart) / number

    # pytransit
    m = pytransit.MandelAgol(interpolate=False, nthr=0)
    tstart = time.time()
    for k in range(number):
        pytransit_flux = m(b, rplanet, [u1, u2])
    pytransit_time[i] = (time.time() - tstart) / number

    # Multiprecision
    if i == 4:
        star_multi = Primary(multi=True)
        star_multi[1] = u1
        star_multi[2] = u2
        planet_multi = Secondary(multi=True)
        planet_multi.r = rplanet
        planet_multi.inc = inc
        planet_multi.porb = P
        planet_multi.a = a
        planet_multi.lambda0 = 90
        system_multi = System(star_multi, planet_multi)
        system_multi.compute(t)
        flux_multi = np.array(system_multi.lightcurve)
        err_agol = np.nanmedian(np.abs(agol_flux - flux_multi))
        err_pytransit = np.nanmedian(np.abs(pytransit_flux - flux_multi))
        err_batman = np.nanmedian(np.abs(batman_flux - flux_multi))


# Plot
fig = pl.figure(figsize=(7, 4))
ax = pl.subplot2grid((2, 5), (0, 0), colspan=4, rowspan=2)
axleg1 = pl.subplot2grid((2, 5), (0, 4))
axleg2 = pl.subplot2grid((2, 5), (1, 4))
axleg1.axis('off')
axleg2.axis('off')

ax.plot(Narr, batman_time, 'o', ms=ms(err_batman), color='C1')
ax.plot(Narr, batman_time, '-', lw=0.75, color='C1')
ax.plot(Narr, agol_time, 'o', ms=ms(err_agol), color='C0')
ax.plot(Narr, agol_time, '-', lw=0.75, color='C0')
ax.plot(Narr, agol_grad_time, 'o', ms=ms(err_agol), color='C0')
ax.plot(Narr, agol_grad_time, '--', lw=0.75, color='C0')
ax.plot(Narr, pytransit_time, 'o', ms=ms(err_pytransit), color='C4')
ax.plot(Narr, pytransit_time, '-', lw=0.75, color='C4')

# Tweak and save
ax.set_ylabel("Time [s]", fontsize=10)
ax.set_xlabel("Number of points", fontsize=10)
ax.set_xscale('log')
ax.set_yscale('log')

# Legend
axleg1.plot([0, 1], [0, 1], color='C0', label='this work', lw=1.5)
axleg1.plot([0, 1], [0, 1], '--', color='C0', label='this work\n(+ gradients)', lw=1.5)
axleg1.plot([0, 1], [0, 1], color='C4', label='PyTransit', lw=1.5)
axleg1.plot([0, 1], [0, 1], color='C1', label='batman', lw=1.5)
axleg1.set_xlim(2, 3)
leg = axleg1.legend(loc='center', frameon=False, fontsize=8)
leg.set_title('method', prop={'weight': 'bold'})
for logerr in [-16, -12, -8, -4, 0]:
    axleg2.plot([0, 1], [0, 1], 'o', color='gray',
                ms=ms(10 ** logerr),
                label=r'$%3d$' % logerr)
axleg2.set_xlim(2, 3)
leg = axleg2.legend(loc='center', labelspacing=1, frameon=False)
leg.set_title('log error', prop={'weight': 'bold'})

# Print average ratio
print(np.nanmedian(agol_grad_time / batman_time))
print(np.nanmedian(agol_time / batman_time))
fig.savefig("compare_to_batman.pdf", bbox_inches='tight')
