"""NUMERICAL: degree of the contact map W -> E_126 (contact point with Kapranov plane 126 = {X3=X4=X5},
at t = u = infinity in the gauge).  Data: lam2, h3 (leading coeffs) and normal ratio rho."""
import numpy as np, warnings, sys
warnings.filterwarnings('ignore')
from kapcurve import curve, f_r23, f_lam2
from kapcheck import summary
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
def data(p):
    c, res, tv = curve(p)
    h = c[0]
    R25, R34, R35, W = p
    R23 = f_r23(*p)
    rho = (1 + 0 - R23 - R35) / (R25 + R35 - 1 - R34)
    return np.array([res, f_lam2(*p), h[3], rho]), c
def newton(p, target, it=80):
    try:
        return _newton(p, target, it)
    except (OverflowError, ZeroDivisionError, FloatingPointError, ValueError):
        return None
def _newton(p, target, it=80):
    p = np.array(p, complex)
    for _ in range(it):
        f = data(p)[0] - target
        J = np.zeros((4, 4), complex); hh = 1e-7
        for k in range(4):
            dp = np.zeros(4, complex); dp[k] = hh*(1 + abs(p[k]))
            J[:, k] = (data(p + dp)[0] - data(p - dp)[0]) / (2*dp[k])
        try: st = np.linalg.solve(J, -f)
        except np.linalg.LinAlgError: return None
        if not np.all(np.isfinite(st)): return None
        p = p + st
        if np.linalg.norm(st) < 1e-12*(1 + np.linalg.norm(p)): break
    f = data(p)[0] - target
    if np.linalg.norm(f) > 1e-8*(1 + np.linalg.norm(target)): return None
    return p
def genuine(p):
    try:
        c = data(p)[1]; s = summary(c)
        return s['nplanes'] == 13 and not s['bad'] and s['ncontacts'] == 19 and s['mindist'] > 1e-5
    except Exception: return False
# find a genuine base point: solve DD=0 in w with others random complex
while True:
    p0 = rng.normal(size=4) + 1j*rng.normal(size=4)
    tg = np.zeros(4, complex)
    # only impose resid=0: Newton in w alone
    p = p0.copy()
    ok = False
    for _ in range(60):
        f = data(p)[0][0]; hh = 1e-7
        dp = np.zeros(4, complex); dp[3] = hh
        df = (data(p + dp)[0][0] - data(p - dp)[0][0])/(2*hh)
        p[3] -= f/df
        if abs(f) < 1e-12: ok = True; break
    if ok and genuine(p): break
target = data(p)[0]; target[0] = 0
print('base point', np.round(p, 6), 'target', np.round(target, 6))
sols = [p]
nconv = 0; nbase = 0; ndeg = 0
for k in range(int(sys.argv[2]) if len(sys.argv) > 2 else 300):
    q = newton(rng.normal(size=4)*2 + 1j*rng.normal(size=4)*2, target)
    if q is None: continue
    nconv += 1
    if np.linalg.norm(q - p) < 1e-6*(1+np.linalg.norm(p)): nbase += 1
    if not genuine(q): ndeg += 1; continue
    if all(np.linalg.norm(q - s) > 1e-6*(1 + np.linalg.norm(s)) for s in sols):
        sols.append(q); print('new genuine solution', np.round(q, 6), flush=True)
print('converged', nconv, 'to base', nbase, 'nongenuine', ndeg)
print('number of genuine points in the fibre:', len(sols))
