"""NUMERICAL: on coordinate lines, roots of each numerator factor of DD -> genuine cub13?"""
import sympy as sp, numpy as np, pickle, warnings
warnings.filterwarnings('ignore')
from kapcurve import curve, SY
from kapcheck import summary
r25, r34, r35, w = SY
d = pickle.load(open('kap2.pkl', 'rb'))
vals = d['vals']
def DD(sub):
    tv = [(sp.cancel(a.subs(sub)), sp.cancel(b.subs(sub))) for a, b in vals]
    s = 0
    for k, (tk, yk) in enumerate(tv):
        den = 1
        for j, (tj, _) in enumerate(tv):
            if j != k: den *= (tk - tj)
        s += yk/den
    return sp.factor(sp.together(s))
base = {r25: sp.Rational(3), r34: sp.Rational(-2), r35: sp.Rational(5, 3), w: sp.Rational(7, 2)}
for var in SY:
    sub = {k: vv for k, vv in base.items() if k != var}
    n, dn = sp.fraction(DD(sub))
    for f, m in sp.factor_list(n)[1]:
        rts = sp.Poly(f, var).nroots(n=30)
        gen = 0
        for rt in rts:
            p = [complex(base[x]) if x != var else complex(rt) for x in SY]
            try:
                c, res, tv = curve(p)
                s = summary(c)
                ok = (s['nplanes'] == 13 and not s['bad'] and s['ncontacts'] == 19 and s['mindist'] > 1e-4)
            except Exception as e:
                ok = False; s = {'err': str(e)}
            gen += ok
        print(var, 'factor deg', sp.degree(f, var), 'mult', m, 'genuine roots', gen, '/', len(rts), flush=True)
