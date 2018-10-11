"""Fitting a polynomial to the nonlinear limb darkening model."""
import matplotlib.pyplot as pl
import numpy as np
from scipy.optimize import curve_fit
cmap = pl.get_cmap("plasma")

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
    X = np.vander((1 - mu), N=order, increasing=True)
    guess = -np.linalg.solve(np.dot(X.transpose(), X), np.dot(X.transpose(), I))[1:]
    u, _ = curve_fit(Polynomial, mu, I, guess, jac=PolynomialJac)
    IPoly = Polynomial(mu, *u)
    err = np.sum((I - IPoly) ** 2) / N
    return u, err


# Non-linear coefficients
c = [0.2, 0.2, 0.2, 0.2]

# Mu array and non-linear intensity
mu = np.linspace(0, 1, 1000)
INL = NonLinear(mu, *c)

# Plot the error for each polynomial order from 2 - 200
orders = np.logspace(np.log10(2), np.log10(500), 10)
orders = np.array(orders, dtype=int)
fig, ax = pl.subplots(1)
for i, order in enumerate(orders):
    u, err = GetPolynomialCoeffs(c, order)
    IPoly = Polynomial(mu, *u)
    pl.plot(mu, np.abs(INL - IPoly), color=cmap(float(i) / len(orders)), label=order)

pl.legend(loc="lower left", ncol=2, fontsize=8)
pl.yscale("log")
pl.xlabel(r"$\mu$", fontsize=16)
pl.ylabel(r"Error", fontsize=16)
fig.savefig("nonlinear_fit.pdf", bbox_inches="tight")
