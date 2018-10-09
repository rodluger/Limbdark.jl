"""High order limb-darkening with starry."""
import starry
import numpy as np
import matplotlib.pyplot as pl
from scipy.optimize import curve_fit


def IofMu(mu, *u):
    """Return the specific intensity as a function of `mu`."""
    return (1 - np.sum([u[l] * (1 - mu) ** l for l in range(1, len(u))],
        axis=0))


# Say we have the following specific intensity values from a stellar model:
mu = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]
I = [1.0, 0.95904139, 0.91381385, 0.86332286, 0.80617997, 0.74036269,
     0.66275783, 0.56820172, 0.44715803, 0.2787536,  0.0]

# Second-order model
guess = [1, 0, 0]
u, _ = curve_fit(IofMu, mu, I, guess)

# Tenth-order model
guess = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
u10, _ = curve_fit(IofMu, mu, I, guess)

# Compute transit light curves
npts = 10000
r = 0.1
b = np.linspace(-1.5, 1.5, npts)

# Second-order
map = starry.Map()
map[1] = u[1]
map[2] = u[2]
flux = map.flux(xo=b, yo=0, ro=r)

# Tenth-order
map10 = starry.Map(10)
for l in range(1, len(u10)):
    map10[l] = u10[l]
flux10 = map10.flux(xo=b, yo=0, ro=r)

# Plot
fig = pl.figure(figsize=(10, 4))
fig.subplots_adjust(wspace=0.35, hspace=0.1)
axI = pl.subplot2grid((2, 2), (0, 0), rowspan=2)
axI.plot(mu, I, 'ko', label='Data')
mu_hires = np.linspace(0, 1, 100)
axI.plot(mu_hires, IofMu(mu_hires, *u), label='Quadratic LD')
axI.plot(mu_hires, IofMu(mu_hires, *u10), label='10th order LD')
axI.set_xlim(1, 0)
axI.set_xlabel(r'$\mu$', fontsize=16)
axI.set_ylabel('Specific Intensity', fontsize=16)
axI.legend(loc='lower left')
axF = pl.subplot2grid((2, 2), (0, 1))
axR = pl.subplot2grid((2, 2), (1, 1))
axF.plot(b, flux, '-', label='Quadratic LD')
axF.plot(b, flux10, '-', label='10th order LD')
axR.plot(b, flux - flux10, '-')
axR.set_xlim(-1.5, 1.5)
axF.set_xlim(-1.5, 1.5)
axF.set_ylabel('Flux', fontsize=16, labelpad=16)
axR.set_ylabel('Relative error', fontsize=16)
axR.set_xlabel('Impact parameter', fontsize=16)
axF.set_xticklabels(())
fig.savefig('high_order_ld.pdf', bbox_inches='tight')
