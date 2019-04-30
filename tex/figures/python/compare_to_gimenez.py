"""Gimenez comparison."""
import time
import matplotlib.pyplot as pl
import numpy as np
import subprocess
import pytransit
import os
import sys
sys.path.insert(0, os.path.join(os.getenv("HOME", "~"), "starry"))
import starry
import exoplanet
from scipy.optimize import curve_fit
np.random.seed(43)


# Marker size is proportional to log error
def ms(error):
    return 18 + np.log10(max(1.e-16, error))


def Polynomial(mu, *u):
    """The polynomial limb darkening model."""
    return 1 - np.sum([u[l] * (1 - mu) ** (l + 1) for l in range(len(u))], axis=0)


def PolynomialJac(mu, *u):
    """The derivative matrix of the polynomial model."""
    jac = -np.array([(1 - mu) ** (l + 1) for l in range(len(u))]).transpose()
    return jac


def GetPolynomialCoeffs(g):
    """Get the polynomial coefficents."""
    N = 100
    mu = np.linspace(0, 1, N)
    I = GimenezPolynomial(mu, *g)
    guess = g
    u, _ = curve_fit(Polynomial, mu, I, guess, jac=PolynomialJac)
    IPoly = Polynomial(mu, *u)
    err = np.max(np.abs(I - IPoly))
    return u, err


def GimenezPolynomial(mu, *u):
    """The Gimenez polynomial limb darkening model."""
    return 1 - np.sum([u[l] * (1 - mu ** (l + 1)) for l in range(len(u))], axis=0)


def GimenezPolynomialJac(mu, *u):
    """The derivative matrix of the Gimenez polynomial model."""
    jac = -np.array([(1 - mu ** (l + 1)) for l in range(len(u))]).transpose()
    return jac


def GetGimenezPolynomialCoeffs(u):
    """Get the polynomial coefficents that approximate the nonlinear model."""
    N = 100
    mu = np.linspace(0, 1, N)
    I = Polynomial(mu, *u)
    guess = u
    g, _ = curve_fit(GimenezPolynomial, mu, I, guess, jac=GimenezPolynomialJac)
    IPoly = GimenezPolynomial(mu, *g)
    err = np.max(np.abs(I - IPoly))
    return g, err


Narr = [1, 3, 5, 10, 15, 20, 30, 40, 50]
b = np.linspace(0.0, 1.1, 1000)

agol_time = np.empty(len(Narr))
agol_grad_time = np.empty(len(Narr))
starry_ylm_time = np.empty(len(Narr))
starry_time = np.empty(len(Narr))
starry_grad_time = np.empty(len(Narr))
pytransit_time = np.empty(len(Narr))
err_pytransit = np.empty(len(Narr))
err_agol = np.empty(len(Narr))
err_starry = np.empty(len(Narr))
err_starry_ylm = np.empty(len(Narr))

# Loop over polynomial degree
for i, N in enumerate(Narr):

    u = np.random.randn(N) * np.exp(-np.arange(N) / 5)
    u_g, err = GetGimenezPolynomialCoeffs(u)
    print(err)

    # Feed b(t) to julia
    # HACK: PyJulia is currently broken, so this is how we have to do this...
    np.savetxt("b.txt", X=b)
    np.savetxt("u.txt", X=u)
    foo = subprocess.check_output(['julia', "compare_to_batman.jl"])
    agol_time[i] = float(foo.decode('utf-8'))
    agol_flux = np.loadtxt("flux.txt")
    flux_multi = np.loadtxt("flux_multi.txt")
    foo = subprocess.check_output(['julia', "compare_to_batman_grad.jl"])
    agol_grad_time[i] = float(foo.decode('utf-8'))

    # pytransit
    m = pytransit.Gimenez(nldc=len(u_g), interpolate=False)
    tstart = time.time()
    for k in range(10):
        pytransit_flux = m(b, 0.1, u_g)
    pytransit_time[i] = (time.time() - tstart) / 10

    # starry (sph)
    # Using the dense spherical harmonic
    # integration algorithm (slow, since we don't
    # actually need the majority of the terms!)
    if N < 30:
        map = starry.Map(ydeg=1, udeg=N)
        map[1:] = u
        tstart = time.time()
        for k in range(10):
            starry_ylm_flux = map.flux(xo=b, ro=0.1)
        starry_ylm_time[i] = (time.time() - tstart) / 10
        tstart = time.time()
    else:
        # Let's not even bother computing these!
        starry_ylm_flux = np.zeros_like(b) * np.nan
        starry_ylm_time[i] = np.nan

    # starry (ld)
    map = starry.Map(ydeg=0, udeg=N)
    map[1:] = u
    tstart = time.time()
    for k in range(10):
        starry_flux = map.flux(b=b, ro=0.1)
    starry_time[i] = (time.time() - tstart) / 10
    tstart = time.time()
    for k in range(10):
        map.flux(b=b, ro=0.1, bf=np.ones_like(b))
    starry_grad_time[i] = (time.time() - tstart) / 10

    # Multiprecision
    err_agol[i] = np.nanmedian(np.abs(agol_flux - flux_multi))
    err_pytransit[i] = np.nanmedian(np.abs(pytransit_flux - flux_multi))
    err_starry[i] = np.nanmedian(np.abs(starry_flux - flux_multi))
    err_starry_ylm[i] = np.nanmedian(np.abs(starry_ylm_flux - flux_multi))

# Plot
fig = pl.figure(figsize=(7, 4))
ax = pl.subplot2grid((2, 5), (0, 0), colspan=4, rowspan=2)
axleg1 = pl.subplot2grid((2, 5), (0, 4))
axleg2 = pl.subplot2grid((2, 5), (1, 4))
axleg1.axis('off')
axleg2.axis('off')

for i in range(len(Narr)):
    ax.plot(Narr[i], agol_time[i], 'o', ms=ms(err_agol[i]), color='C0')
ax.plot(Narr, agol_time, '-', lw=0.75, color='C0')
for i in range(len(Narr)):
    ax.plot(Narr[i], agol_grad_time[i], 'o', ms=ms(err_agol[i]), color='C0')
ax.plot(Narr, agol_grad_time, '--', lw=0.75, color='C0')

for i in range(len(Narr)):
    ax.plot(Narr[i], starry_time[i], 'o', ms=ms(err_starry[i]), color='C2')
ax.plot(Narr, starry_time, '-', lw=0.75, color='C2')
for i in range(len(Narr)):
    ax.plot(Narr[i], starry_grad_time[i], 'o', ms=ms(err_starry[i]), color='C2')
ax.plot(Narr, starry_grad_time, '--', lw=0.75, color='C2')

for i in range(len(Narr)):
    ax.plot(Narr[i], starry_ylm_time[i], 'o', ms=ms(err_starry_ylm[i]), color='C3')
ax.plot(Narr, starry_ylm_time, '-', lw=0.75, color='C3')

for i in range(len(Narr)):
    ax.plot(Narr[i], pytransit_time[i], 'o', ms=ms(err_pytransit[i]), color='C4')
ax.plot(Narr, pytransit_time, '-', lw=0.75, color='C4')

# Tweak and save
ax.set_ylabel("Time [s]", fontsize=10)
ax.set_xlabel("Degree of limb darkening", fontsize=10)
ax.set_yscale('log')
#ax.set_ylim(5e-5, 1e-1)

# Legend
axleg1.plot([0, 1], [0, 1], color='C0', label='limbdark', lw=1.5)
axleg1.plot([0, 1], [0, 1], '--', color='C0', label='limbdark\n(+ gradients)', lw=1.5)
axleg1.plot([0, 1], [0, 1], color='C2', label='starry', lw=1.5)
axleg1.plot([0, 1], [0, 1], '--', color='C2', label='starry\n(+ gradients)', lw=1.5)
axleg1.plot([0, 1], [0, 1], color='C3', label='starry\n(dense)', lw=1.5)
axleg1.plot([0, 1], [0, 1], color='C4', label='PyTransit', lw=1.5)
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

# Save
fig.savefig("compare_to_gimenez.pdf", bbox_inches='tight')