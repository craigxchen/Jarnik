"""NUMERICAL: Newton-homotopy path tracking for Q4 members (gauge of q4b.py)."""
import numpy as np, sys, warnings
warnings.filterwarnings('ignore')
from q4b import eqs, genuine, coords, NAMES
rng = np.random.default_rng(int(sys.argv[1]))
def jac(u, var):
    J = np.zeros((12, len(var)), complex)
    for a, k in enumerate(var):
        du = np.zeros(16, complex); du[k] = 1e-7 * (1 + abs(u[k]))
        J[:, a] = (eqs(u + du) - eqs(u - du)) / (2 * du[k])
    return J
def track(u0, fixed, nsteps=400):
    var = [k for k in range(16) if k not in fixed]
    F0 = eqs(u0); u = np.array(u0, complex)
    gamma = np.exp(2j * np.pi * rng.random())
    s = 0.0; ds = 1.0 / nsteps
    while s < 1.0:
        s1 = min(1.0, s + ds)
        # predictor: Euler  dH/du du + dH/ds ds = 0, H = F(u) - (1-s) gamma F0
        J = jac(u, var)
        try: du = np.linalg.solve(J, -gamma * F0 * (s1 - s))
        except np.linalg.LinAlgError: return None
        v = u.copy(); v[var] += du
        ok = False
        for _ in range(8):
            H = eqs(v) - (1 - s1) * gamma * F0
            try: st = np.linalg.solve(jac(v, var), -H)
            except np.linalg.LinAlgError: break
            v[var] += st
            if np.linalg.norm(st) < 1e-10 * (1 + np.linalg.norm(v)): ok = True; break
        if not ok or not np.all(np.isfinite(v)) or np.linalg.norm(v) > 1e6:
            ds /= 2
            if ds < 1e-7: return None
            continue
        u, s = v, s1
        ds = min(ds * 1.5, 0.05)
    # final polish
    for _ in range(20):
        try: st = np.linalg.solve(jac(u, var), -eqs(u))
        except np.linalg.LinAlgError: return None
        u[var] += st
        if np.linalg.norm(st) < 1e-13 * (1 + np.linalg.norm(u)): break
    return u if np.linalg.norm(eqs(u)) < 1e-9 * (1 + np.linalg.norm(u)) ** 4 else None
ngen = 0; nend = 0
for trial in range(int(sys.argv[2])):
    u0 = rng.normal(size=16) + 1j * rng.normal(size=16)
    fixed = [0, 7, 10, 12]   # r12, r45, t134, t135 fixed
    u = track(u0, fixed)
    if u is None: continue
    nend += 1
    ok, cnt, bad, md = genuine(u)
    if ok:
        ngen += 1; np.save('q4h_member_%s_%d.npy' % (sys.argv[1], ngen), u)
        print('GENUINE', np.round(u, 4), 'min contact distance %.3g' % md, flush=True)
print('paths ended', nend, 'genuine', ngen)
