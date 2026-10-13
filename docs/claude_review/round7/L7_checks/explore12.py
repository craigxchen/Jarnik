"""param2 (b, tau, mu, sp): exact numerator in each variable, rational roots stripped, genuine roots counted."""
import numpy as np, warnings
warnings.filterwarnings('ignore')
from fractions import Fraction as Fr
import qpoly as Q, param2 as P2
from modp import M, ratdeg
from explore10 import exact
from explore7_lib import genuine
base = dict(b=Fr(3), tau=Fr(-2), mu=Fr(5, 2), sp=Fr(7, 3))
order = ['b', 'tau', 'mu', 'sp']
for v in order:
    fm = lambda x: P2.stage2(P2.stage1(*[x if w == v else M(base[w]) for w in order]))
    dn, dd, _, _ = ratdeg(fm, 40)
    f = lambda x: P2.stage2(P2.stage1(*[x if w == v else base[w] for w in order]))
    num, den = exact(f, dn, dd)
    rr = Q.rational_roots(num)
    rest = num; mult = {}
    for r in rr:
        while True:
            qq, rem = Q.divmod_(rest, Q.lin(r))
            if rem: break
            rest = qq; mult[str(r)] = mult.get(str(r), 0) + 1
    roots = np.polynomial.polynomial.polyroots(np.array([float(x) for x in rest]))
    gen = []
    for r in roots:
        args = [complex(r) if w == v else complex(base[w]) for w in order]
        try:
            d = P2.stage1(*args); P2.stage2(d)
            if genuine(d): gen.append(r)
        except Exception: pass
    print(v, 'num', dn, 'den', dd, 'rational roots', mult, 'rest deg', Q.deg(rest), 'genuine', len(gen), np.round(gen, 4), flush=True)
