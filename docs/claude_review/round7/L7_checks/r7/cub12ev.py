"""NUMERICAL: degree of ev_{s0}: V -> X for the Q-rational cub12 family V (params r23,r25,r34,r35;
gauge t126=inf, t123=0 (r45), t135=1 (r24)).  Count v with c(s0; v) proportional to a random x."""
import numpy as np, sys, warnings
warnings.filterwarnings('ignore')
import cub13exact as C
rng = np.random.default_rng(int(sys.argv[1]))
def build12(v):
    r23, r25, r34, r35 = v
    r24, r45 = 1.0 + 0j, 0j
    vv = (r25 * r35 - r24 * r34) / (r25 + r35 - r24 - r34)
    X4 = C.prod_lin([r24, r34, r45], 1 + 0j); X5 = C.prod_lin([r25, r35, r45], 1 + 0j)
    X3 = C.prod_lin([r23, r34, r35], 1 + 0j)
    lam2 = C.ev(X4, vv) / ((vv - r23) * (vv - r24) * (vv - r25))
    X2 = C.prod_lin([r23, r24, r25], lam2)
    def other(P, known, deg):
        s = -P[deg - 1] / P[deg]
        for k in known: s = s - k
        return s
    D35 = C.padd(X3, C.pscale(-1, X5)); D34 = C.padd(X3, C.pscale(-1, X4))
    D25 = C.padd(X2, C.pscale(-1, X5)); D24 = C.padd(X2, C.pscale(-1, X4))
    t246 = other(D35, [r35], 2); t256 = other(D34, [r34], 2)
    t346 = other(D25, [r25, vv], 3); t356 = other(D24, [r24, vv], 3)
    pts = [(t246, C.ev(X3, t246)), (t256, C.ev(X3, t256)), (t346, C.ev(X2, t346)), (t356, C.ev(X2, t356))]
    xs = np.array([p[0] for p in pts]); ys = np.array([p[1] for p in pts])
    h = np.linalg.solve(np.vander(xs, 4, increasing=True), ys)
    return [h, X2, X3, X4, X5]
def F(v, s0, x):
    cs = build12(v)
    val = np.array([np.polyval(np.array(c)[::-1], s0) for c in cs])
    return val[:4] / val[4] - x
def newton(v, s0, x):
    v = np.array(v, complex)
    for _ in range(100):
        f = F(v, s0, x); J = np.zeros((4, 4), complex)
        for k in range(4):
            dv = np.zeros(4, complex); dv[k] = 1e-7 * (1 + abs(v[k]))
            J[:, k] = (F(v + dv, s0, x) - F(v - dv, s0, x)) / (2 * dv[k])
        st = np.linalg.solve(J, -f)
        if not np.all(np.isfinite(st)): return None
        v = v + st
        if np.linalg.norm(st) < 1e-12 * (1 + np.linalg.norm(v)): break
    if np.linalg.norm(F(v, s0, x)) > 1e-9: return None
    return v
v0 = rng.normal(size=4) + 1j * rng.normal(size=4)
s0 = 0.37 - 0.21j
x = F(v0, s0, np.zeros(4))   # a point on the curve of v0 at s0 (so a solution exists)
sols = [v0]; conv = 0
for k in range(int(sys.argv[2])):
    try:
        v = newton(v0 + (rng.normal(size=4) + 1j * rng.normal(size=4)) * rng.choice([0.3, 1, 3]), s0, x)
    except Exception:
        v = None
    if v is None: continue
    conv += 1
    if all(np.linalg.norm(v - w) > 1e-6 * (1 + np.linalg.norm(w)) for w in sols):
        sols.append(v); print('new solution', np.round(v, 5), flush=True)
print('converged', conv, 'distinct solutions (incl. base)', len(sols))
