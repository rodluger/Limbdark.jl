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
    nb = 1001
    nr = 1001
    b = np.linspace(0.0, 2.0, nb, endpoint=True)
    r = np.linspace(0.0, 2.0, nr, endpoint=True)
    fgrid = np.zeros((nr, nb))
    for ib in range(nb):
        for ir in range(nr):
            fgrid[ir, ib] = s2(r[ir], b[ib])
    np.savetxt("s2_MA2002.txt", X=fgrid)