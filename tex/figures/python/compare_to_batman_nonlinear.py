"""Starry speed tests."""
from starry.kepler import Primary, Secondary, System
import time
import matplotlib.pyplot as pl
import numpy as np
import batman
import subprocess
from scipy.optimize import curve_fit
from scipy.special import gamma
from scipy.integrate import dblquad
import pytransit


def ms(error):
    """Marker size, proportional to log error."""
    return 18 + np.log10(error)


def NonLinear(mu, *c):
    """The nonlinear limb darkening model."""
    return 1 - c[0] * (1 - mu ** 0.5) \
             - c[1] * (1 - mu ** 1.0) \
             - c[2] * (1 - mu ** 1.5) \
             - c[3] * (1 - mu ** 2.0)


def Polynomial(mu, *u):
    """The polynomial limb darkening model."""
    return 1 - np.sum([u[l] * (1 - mu) ** (l + 1) for l in range(len(u))], axis=0)


def PolynomialJac(mu, *u):
    """The derivative matrix of the polynomial model."""
    jac = -np.array([(1 - mu) ** (l + 1) for l in range(len(u))]).transpose()
    return jac


def GetPolynomialCoeffs(c, order):
    """Get the polynomial coefficents that approximate the nonlinear model."""
    N = 1000
    mu = np.linspace(0, 1, N)
    I = NonLinear(mu, *c)
    X = np.vander((1 - mu), N=order + 1, increasing=True)
    guess = -np.linalg.solve(np.dot(X.transpose(), X), np.dot(X.transpose(), I))[1:]
    u, _ = curve_fit(Polynomial, mu, I, guess, jac=PolynomialJac)
    IPoly = Polynomial(mu, *u)
    err = np.sum((I - IPoly) ** 2) / N
    return u, err


def GimenezPolynomial(mu, *u):
    """The Gimenez polynomial limb darkening model."""
    return 1 - np.sum([u[l] * (1 - mu ** (l + 1)) for l in range(len(u))], axis=0)


def GimenezPolynomialJac(mu, *u):
    """The derivative matrix of the Gimenez polynomial model."""
    jac = -np.array([(1 - mu ** (l + 1)) for l in range(len(u))]).transpose()
    return jac


def GetGimenezPolynomialCoeffs(c, order):
    """Get the polynomial coefficents that approximate the nonlinear model."""
    N = 1000
    mu = np.linspace(0, 1, N)
    I = NonLinear(mu, *c)
    X = np.vander((1 - mu), N=order + 1, increasing=True)
    guess = -np.linalg.solve(np.dot(X.transpose(), X), np.dot(X.transpose(), I))[1:]
    u, _ = curve_fit(GimenezPolynomial, mu, I, guess, jac=GimenezPolynomialJac)
    IPoly = GimenezPolynomial(mu, *u)
    err = np.sum((I - IPoly) ** 2) / N
    return u, err


def NumericalFlux(b, r, c):
    """Compute the flux by numerical integration of the surface integral."""
    # I'm only coding up a specific case here
    assert r <= 1, "Invalid range."
    if b < 0:
        b = np.abs(b)

    # No occ
    if b >= 1 + r:
        return 1

    # Get points of intersection
    if b > 1 - r:
        yi = (1. + b ** 2 - r ** 2) / (2. * b)
        xi = (1. / (2. * b)) * np.sqrt(4 * b ** 2 - (1 + b ** 2 - r ** 2) ** 2)
    else:
        yi = np.inf
        xi = r

    # Specific intensity map
    def I(y, x):
        mu = np.sqrt(1 - x ** 2 - y ** 2)
        return 1 - c[0] * (1 - mu ** 0.5) - c[1] * (1 - mu) - c[2] * (1 - mu ** 1.5) - c[3] * (1 - mu ** 2)

    # Total flux
    total, _ = dblquad(I, -1, 1, lambda x: 0, lambda x: np.sqrt(1 - x ** 2), epsabs=1e-12, epsrel=1e-12)
    total *= 2

    # Lower integration limit
    def y1(x):
        if yi <= b:
            # Lower occultor boundary
            return b - np.sqrt(r ** 2 - x ** 2)
        elif b <= 1 - r:
            # Lower occultor boundary
            return b - np.sqrt(r ** 2 - x ** 2)
        else:
            # Tricky: we need to do this in two parts
            return b - np.sqrt(r ** 2 - x ** 2)

    # Upper integration limit
    def y2(x):
        if yi <= b:
            # Upper occulted boundary
            return np.sqrt(1 - x ** 2)
        elif b <= 1 - r:
            # Upper occultor boundary
            return b + np.sqrt(r ** 2 - x ** 2)
        else:
            # Tricky: we need to do this in two parts
            return np.sqrt(1 - x ** 2)

    # Compute the total flux
    flux, _ = dblquad(I, -xi, xi, y1, y2, epsabs=1e-12, epsrel=1e-12)

    # Do we need to solve an additional integral?
    if not (yi <= b) and not (b <= 1 - r):

        def y1(x):
            return b - np.sqrt(r ** 2 - x ** 2)

        def y2(x):
            return b + np.sqrt(r ** 2 - x ** 2)

        additional_flux, _ = dblquad(I, -r, -xi, y1, y2,
                                     epsabs=1e-12, epsrel=1e-12)

        flux += 2 * additional_flux

    return (total - flux) / total


# Input params
c = [0.2, 0.2, 0.2, 0.2]
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

# Get the polynomial coeffs for l = 15
order = 15
u, err = GetPolynomialCoeffs(c, order)
print("Polynomial fit error: %.3e" % err)

# Get the Gimenez polynomial coeffs
u_g, err = GetGimenezPolynomialCoeffs(c, order)
print("Polynomial fit error (Gimenez): %.3e" % err)

# Timing params
number = 30
nN = 8
Nmax = 5
Narr = np.logspace(1, Nmax, nN)
agol_time = np.zeros(nN) * np.nan
agol_grad_time = np.zeros(nN) * np.nan
pytransit_time = np.zeros(nN)
batman_time = np.zeros(nN)

# Loop over number of cadences
for i, N in enumerate(Narr):

    # Time array
    t = np.linspace(-0.15, 0.15, N)

    # Use starry to compute b(t)
    star = Primary()
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
    np.savetxt("u.txt", X=u)
    foo = subprocess.check_output(['julia', "compare_to_batman.jl"])
    agol_time[i] = float(foo.decode('utf-8'))
    agol_flux = np.loadtxt("flux.txt")
    foo = subprocess.check_output(['julia', "compare_to_batman_grad.jl"])
    agol_grad_time[i] = float(foo.decode('utf-8'))

    # batman
    params = batman.TransitParams()
    params.limb_dark = "nonlinear"
    params.u = c
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
    if (N < 100000):
        m = pytransit.Gimenez(nldc=len(u_g), interpolate=False, nthr=0)
        tstart = time.time()
        for k in range(number):
            pytransit_flux = m(b, 0.1, u_g)
        pytransit_time[i] = (time.time() - tstart) / number
    else:
        # HACK: PyTransit segfaults on Travis for N = 100000
        # Let's be generous and linearly extrapolate
        pytransit_flux = b * np.nan
        pytransit_time[i] = pytransit_time[i - 1] * Narr[i] / Narr[i - 1]

    # Multiprecision
    if i == 1:
        flux_multi = [NumericalFlux(bi, 0.1, c) for bi in b]
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

ax.plot(Narr, agol_time, 'o', ms=ms(err_agol), color='C0')
ax.plot(Narr, agol_time, '-', lw=0.75, color='C0')
ax.plot(Narr, agol_grad_time, 'o', ms=ms(err_agol), color='C0')
ax.plot(Narr, agol_grad_time, '--', lw=0.75, color='C0')
ax.plot(Narr, pytransit_time, 'o', ms=ms(err_pytransit), color='C4')
ax.plot(Narr, pytransit_time, '-', lw=0.75, color='C4')
ax.plot(Narr, batman_time, 'o', ms=ms(err_batman), color='C1')
ax.plot(Narr, batman_time, '-', lw=0.75, color='C1')

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

# Print average time and error ratios
print(np.nanmedian(agol_time / batman_time))
print(err_agol / err_batman)
fig.savefig("compare_to_batman_nonlinear.pdf", bbox_inches='tight')
