"""DD restricted to lines: factor structure in each variable (exact)."""
import sympy as sp, pickle, random
t, r25, r34, r35, w = sp.symbols('t r25 r34 r35 w')
d = pickle.load(open('kap2.pkl', 'rb'))
vals = d['vals']
def DD(vs, sub):
    tv = [(sp.cancel(a.subs(sub)), sp.cancel(b.subs(sub))) for a, b in vs]
    s = 0
    for k, (tk, yk) in enumerate(tv):
        den = 1
        for j, (tj, _) in enumerate(tv):
            if j != k: den *= (tk - tj)
        s += yk/den
    return sp.factor(sp.together(s))
random.seed(1)
base = {r25: sp.Rational(3), r34: sp.Rational(-2), r35: sp.Rational(5, 3), w: sp.Rational(7, 2)}
for var in [r25, r34, r35, w]:
    sub = {k: vv for k, vv in base.items() if k != var}
    e = DD(vals, sub)
    n, dn = sp.fraction(e)
    fl = sp.factor_list(n)
    print(var, 'numerator factors (deg):', [(sp.degree(f, var), m) for f, m in fl[1]])
    print('   denominator factors:', [(sp.degree(f, var), m) for f, m in sp.factor_list(dn)[1]], flush=True)
