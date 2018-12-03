import occultquad
import numpy as np


def s2(r, b):
    if (b >= 1 + r):
        return (2 * np.pi / 3)
    elif (b <= r - 1):
        return 0
    r_ = np.array([r], dtype=float)
    b_ = np.array([b], dtype=float)
    flux = np.empty(1)
    mu0 = np.empty(1)
    try:
        occultquad.occultquad(b_, 1.0, 0.0, r_, flux, mu0)
    except:
        print("FOO")
    return (2 * np.pi / 3) * flux[0]


if __name__ == "__main__":
    # Sanity checks
    import starry
    import numpy as np
    import matplotlib.pyplot as pl
    map = starry.Map()
    map[0, 0] = 0
    map[1, 0] = np.pi / np.sqrt(3)
    xo = np.linspace(-1.5, 1.5, 100)
    ro = 0.1
    flux_starry = map.flux(xo=xo, ro=ro)
    flux_ma = [s2(ro, np.abs(b)) for b in xo]
    pl.plot(xo, flux_starry)
    pl.plot(xo, flux_ma, '--')
    pl.show()