"""curve of the Kapranov cub13 parametrisation (kap1/kap2) at given parameter values.
Works numerically (complex) via lambdify; exact version via sympy subs."""
import sympy as sp, pickle, numpy as np
t, r25, r34, r35, w = sp.symbols('t r25 r34 r35 w')
d = pickle.load(open('/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round7/L7_checks/r7/kap2.pkl', 'rb'))
vals, r23e, lam2e = d['vals'], d['r23'], d['lam2']
r24, r45 = 1, 0
SY = (r25, r34, r35, w)
f_r23 = sp.lambdify(SY, r23e, 'numpy')
f_lam2 = sp.lambdify(SY, lam2e, 'numpy')
f_vals = [(sp.lambdify(SY, a, 'numpy'), sp.lambdify(SY, b, 'numpy')) for a, b in vals]
Pn = np.polynomial.polynomial
def L(r): return np.array([-r, 1], complex)
def curve(p, which=(0, 1, 2, 3)):
    """returns c = [h, X2, X3, X4, X5] (low-first coefficient arrays) and the residual of the 5th condition"""
    p = [complex(x) for x in p]
    R23 = f_r23(*p); lam2 = f_lam2(*p)
    R25, R34, R35, W = p
    X2 = lam2*Pn.polymul(Pn.polymul(L(R23), L(r24)), L(R25))
    X3 = Pn.polymul(Pn.polymul(L(R23), L(R34)), L(R35))
    X4 = Pn.polymul(Pn.polymul(L(r24), L(R34)), L(r45))
    X5 = Pn.polymul(Pn.polymul(L(R25), L(R35)), L(r45))
    tv = [(fa(*p), fb(*p)) for fa, fb in f_vals]
    xs = np.array([tv[k][0] for k in which]); ys = np.array([tv[k][1] for k in which])
    V = np.vander(xs, 4, increasing=True)
    h = np.linalg.solve(V, ys)
    k5 = [k for k in range(5) if k not in which][0]
    resid = Pn.polyval(tv[k5][0], h) - tv[k5][1]
    return [h, X2, X3, X4, X5], resid, tv
