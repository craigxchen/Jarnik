"""Exact search for Q-rational points of W = {G = 0}: fix three coordinates to small rationals,
factor DD numerator in the 4th, report rational roots of the genuine (top-degree) factor."""
import sympy as sp, pickle, random, sys, itertools, warnings
warnings.filterwarnings('ignore')
from kapcurve import curve
from kapcheck import summary
t, r25, r34, r35, w = sp.symbols('t r25 r34 r35 w')
SY = [r25, r34, r35, w]
d = pickle.load(open('kap2.pkl', 'rb'))
vals = d['vals']
def DDnum(sub):
    tv = [(sp.cancel(a.subs(sub)), sp.cancel(b.subs(sub))) for a, b in vals]
    s = 0
    for k, (tk, yk) in enumerate(tv):
        den = 1
        for j, (tj, _) in enumerate(tv):
            if j != k: den *= (tk - tj)
        s += yk/den
    return sp.fraction(sp.factor(sp.together(s)))[0]
random.seed(int(sys.argv[1]))
var = SY[int(sys.argv[2])]
others = [x for x in SY if x != var]
smalls = sorted(set(sp.Rational(a, b) for a in range(-6, 7) for b in range(1, 5)) - {0, 1})
found = 0
for trial in range(int(sys.argv[3])):
    sub = {x: random.choice(smalls) for x in others}
    try:
        n = DDnum(sub)
    except Exception as e:
        continue
    if n == 0: continue
    for f, m in sp.factor_list(n)[1]:
        if sp.degree(f, var) < 1: continue
        for r in sp.roots(sp.Poly(f, var), filter='Q'):
            p = [complex(sub[x]) if x != var else complex(r) for x in SY]
            try:
                c, res, tv = curve(p); sm = summary(c)
                ok = sm['nplanes'] == 13 and not sm['bad'] and sm['ncontacts'] == 19 and sm['mindist'] > 1e-6
            except Exception:
                ok = False
            if ok:
                print('GENUINE RATIONAL POINT', {str(k): str(v) for k, v in sub.items()}, str(var), '=', r, flush=True)
                found += 1
print('done, found', found)
