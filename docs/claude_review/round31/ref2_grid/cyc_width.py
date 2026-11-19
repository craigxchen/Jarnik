"""Cyclic Fourier factors: random check of nonvanishing at width < pi/3.
(A search for the smallest vanishing width was dropped: it degenerates to coincident points.)"""
import numpy as np
from scipy.optimize import minimize
rng = np.random.default_rng(0)
xi = np.exp(2j*np.pi/6)
def facs(th):
    z = np.exp(1j*th)
    return np.array([abs(sum(xi**(k*j)*z[j] for j in range(6))) for k in range(1,6)])
mn = np.inf
for t in range(20000):
    w = rng.uniform(1e-3, np.pi/3*0.999)
    th = np.sort(rng.uniform(0, w, 6)); th[0]=0; th[-1]=w
    th = np.sort(th)
    f = facs(th) / np.prod(np.diff(th))**(1/5)
    mn = min(mn, f.min()/w)
print('random ordered sextuples with width<pi/3: min |factor| / width > 0 :', mn > 0, mn)
