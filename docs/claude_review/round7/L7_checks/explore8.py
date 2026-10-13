"""NUMERICAL: zeros of Phi in each variable (Newton from many starts), genuineness test."""
import numpy as np, sys
from fractions import Fraction as Fr
import param2 as P2
from explore7_lib import genuine
rng = np.random.default_rng(5)
base = dict(b=3.0, tau=-2.0, mu=2.5, sp=7/3)
def phi(var, x):
    a = dict(base); a[var] = x
    d = P2.stage1(complex(a['b']), complex(a['tau']), complex(a['mu']), complex(a['sp']))
    return P2.stage2(d), d
for var in ['sp', 'b', 'tau', 'mu']:
    zs = []
    for k in range(400):
        x = complex(rng.normal() * 4, rng.normal() * 4)
        okk = False
        for it in range(100):
            try:
                f, _ = phi(var, x); h = 1e-7 * (1 + abs(x))
                df = (phi(var, x + h)[0] - phi(var, x - h)[0]) / (2 * h)
                st = -f / df
            except ZeroDivisionError:
                break
            x = x + st
            if abs(st) < 1e-14 * (1 + abs(x)): okk = True; break
        if not okk or not np.isfinite(x): continue
        if all(abs(x - z) > 1e-7 * (1 + abs(x)) for z in zs): zs.append(x)
    gen = []
    for z in zs:
        try:
            _, d = phi(var, z)
            if genuine(d): gen.append(z)
        except Exception: pass
    print(var, 'zeros found', len(zs), 'genuine', len(gen), np.round(gen, 5), flush=True)
