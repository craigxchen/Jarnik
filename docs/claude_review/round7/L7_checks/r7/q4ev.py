"""NUMERICAL: degree of ev_{s0}: Q4 -> X (members through a point x with x at intrinsic parameter s0)."""
import numpy as np, sys, warnings, glob
warnings.filterwarnings('ignore')
from q4c import solve_all, coords, FREE, EXPECT
from levelsets import analyse
P = np.polynomial.polynomial
f0v = np.load(sys.argv[1]); f0 = dict(zip(FREE, f0v))
rng = np.random.default_rng(int(sys.argv[2]))
s0 = 0.31 + 0.27j
def point(f, s):
    _, r, lam = solve_all(f); cs = coords(r, lam)
    v = np.array([P.polyval(s, c) for c in cs])
    return v[:4] / v[4]
x = point(f0, s0)
def G(u):
    f = dict(zip(FREE, u))
    e, r, lam = solve_all(f)
    return np.concatenate([e, point(f, s0) - x])
def newton(u):
    u = np.array(u, complex)
    for _ in range(150):
        g = G(u); J = np.zeros((7, 7), complex)
        for k in range(7):
            du = np.zeros(7, complex); du[k] = 1e-7 * (1 + abs(u[k]))
            J[:, k] = (G(u + du) - G(u - du)) / (2 * du[k])
        st = np.linalg.solve(J, -g)
        if not np.all(np.isfinite(st)): return None
        n = np.linalg.norm(st)
        if n > 1: st /= n
        u = u + st
        if n < 1e-13 * (1 + np.linalg.norm(u)): break
    return u if np.linalg.norm(G(u)) < 1e-10 else None
sols = [f0v]; conv = 0; gen = 0
for trial in range(int(sys.argv[3])):
    scale = rng.choice([0.1, 0.5, 1.5, 4.0])
    u0 = f0v + scale * (rng.normal(size=7) + 1j * rng.normal(size=7))
    try: u = newton(u0)
    except Exception: u = None
    if u is None: continue
    conv += 1
    try:
        _, r, lam = solve_all(dict(zip(FREE, u))); cnt, bad, md, infb = analyse(coords(r, lam))
        ok = dict(cnt) == EXPECT and not bad and md > 1e-5 and not infb
    except Exception:
        ok = False
    if not ok: continue
    gen += 1
    if all(np.linalg.norm(u - w) > 1e-6 * (1 + np.linalg.norm(w)) for w in sols):
        sols.append(u); print('new genuine member through x at s0:', np.round(u, 4), flush=True)
print('converged', conv, 'genuine', gen, 'distinct genuine (incl. base)', len(sols))
