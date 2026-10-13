# Evidence only (floating point): real zeros u near 0 of H(e^{iu}-1) for the M8' cubic
# H = x2x3x4 + x2x3x5 - x2x4x5 - x3x4x5 - x2x3x6 + x4x5x6, and a crude check for small integer relations.
import numpy as np, itertools
from scipy.optimize import least_squares
def H(x):
    x2, x3, x4, x5, x6 = x
    return x2*x3*x4 + x2*x3*x5 - x2*x4*x5 - x3*x4*x5 - x2*x3*x6 + x4*x5*x6
def res(u):
    v = H(np.exp(1j*u) - 1) / (np.linalg.norm(u)**3)
    return [v.real, v.imag]
rng = np.random.default_rng(3)
vecs = [np.array(v) for v in itertools.product(range(-3, 4), repeat=5) if any(v)]
good = 0
for trial in range(40):
    u0 = rng.normal(size=5) * 0.05
    sol = least_squares(res, u0, xtol=1e-15, ftol=1e-15, gtol=1e-15)
    u = sol.x
    r = abs(H(np.exp(1j*u) - 1))
    if r < 1e-14 and np.linalg.norm(u) > 1e-3:
        mind = min(abs(v @ u) for v in vecs) / np.linalg.norm(u)
        if mind > 1e-4:
            good += 1
        if trial < 6:
            print(f"|u|={np.linalg.norm(u):.3e} |H|={r:.1e} min_(|v|_inf<=3) |<v,u>|/|u| = {mind:.2e}")
print("solutions with |H|<1e-14 and no relation |v|_inf<=3 (to 1e-4):", good, "of 40 starts")
