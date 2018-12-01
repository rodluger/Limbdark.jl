from pytransit import Gimenez
import numpy as np
np.random.seed(43)


def test(npts, nldc):
    print("Running %d points, %d coeffs... " % (npts, nldc), end="")
    m = Gimenez(nldc=nldc, interpolate=False, nthr=0)
    b = np.linspace(0, 1.5, npts)
    u = np.random.randn(nldc)
    flux = m(b, 0.1, u)
    print("OK!")