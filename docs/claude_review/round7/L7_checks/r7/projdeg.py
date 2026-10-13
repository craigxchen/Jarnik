"""Degree of projections W -> A^3 of the form (r25, r34, phi) on the slice r25=3, r34=-2 (exact elimination
+ numerical genuineness test of the solutions)."""
import sympy as sp, pickle, numpy as np, warnings, random
warnings.filterwarnings('ignore')
from kapcurve import curve
from kapcheck import summary
t, r25, r34, r35, w = sp.symbols('t r25 r34 r35 w')
S = pickle.load(open('slice_r25_3_r34_-2.pkl', 'rb'))
G = S['G']
d = pickle.load(open('kap2.pkl', 'rb'))
sub = {r25: 3, r34: -2}
vals = d['vals']
r23 = sp.cancel(d['r23'].subs(sub)); lam2 = sp.cancel(d['lam2'].subs(sub))
v = sp.cancel(d['v'].subs(sub))
tD23 = sp.cancel(((lam2 * 1 * 3 - (-2) * r35) / (lam2 - 1)) / w)
phis = dict(r23=r23, v=v, lam2=lam2, tD23=tD23,
            t246=sp.cancel(vals[0][0].subs(sub)), t256=sp.cancel(vals[1][0].subs(sub)),
            t346=sp.cancel(vals[2][0].subs(sub)), t356=sp.cancel(vals[3][0].subs(sub)),
            r35=r35, w=w)
def genuine(r35v, wv):
    try:
        c, res, tv = curve([3, -2, complex(r35v), complex(wv)])
        if abs(res) > 1e-6 * (1 + max(abs(x) for x in np.concatenate(c))): return False
        sm = summary(c)
        return sm['nplanes'] == 13 and not sm['bad'] and sm['ncontacts'] == 19 and sm['mindist'] > 1e-6
    except Exception:
        return False
random.seed(5)
for name, phi in phis.items():
    c = sp.Rational(random.randint(-40, 40), random.randint(1, 13))
    num, den = sp.fraction(sp.together(phi - c))
    R = sp.resultant(sp.Poly(G, w), sp.Poly(sp.expand(num), w))
    R = sp.Poly(R, r35)
    sols = 0; cand = 0
    for f, m in sp.factor_list(R.as_expr())[1]:
        fp = sp.Poly(f, r35)
        if fp.degree() < 1: continue
        xs = np.roots([complex(cc) for cc in fp.all_coeffs()])
        GP = sp.Poly(G, w); NP = sp.Poly(sp.expand(num), w)
        Gc = [sp.lambdify(r35, cc) for cc in GP.all_coeffs()]
        Nc = [sp.lambdify(r35, cc) for cc in NP.all_coeffs()]
        for x in xs:
            gco = [complex(fn(x)) for fn in Gc]; nco = [complex(fn(x)) for fn in Nc]
            for y in np.roots(gco):
                if abs(np.polyval(nco, y)) < 1e-6 * (1 + sum(abs(cc) for cc in nco)) * (1 + abs(y))**(len(nco) - 1):
                    cand += 1
                    if genuine(x, y): sols += 1
    print('projection (r25, r34, %s): genuine fibre points = %d  (candidates %d, resultant deg %d)' % (name, sols, cand, R.degree()), flush=True)
