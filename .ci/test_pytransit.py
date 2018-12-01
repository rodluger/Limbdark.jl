"""Test PyTransit."""
from starry.kepler import Primary, Secondary, System
import time
import numpy as np
import pytransit
from scipy.optimize import curve_fit


def NonLinear(mu, *c):
    """The nonlinear limb darkening model."""
    return 1 - c[0] * (1 - mu ** 0.5) \
             - c[1] * (1 - mu ** 1.0) \
             - c[2] * (1 - mu ** 1.5) \
             - c[3] * (1 - mu ** 2.0)


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


def test(order=15, N=100, number=30):
    print("Running %d, %d, %d... " % (order, N, number), end="")
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
    u_g, err = GetGimenezPolynomialCoeffs(c, order)

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

    # pytransit
    m = pytransit.Gimenez(nldc=order, interpolate=False, nthr=0)
    tstart = time.time()
    for k in range(number):
        pytransit_flux = m(b, rplanet, u_g)
    pytransit_time = (time.time() - tstart) / number
    print("OK!")