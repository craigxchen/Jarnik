"""Exact: genuine factor of DD in w at (r25, r34, r35) = given rationals; prints integer coefficients."""
import sympy as sp, pickle, sys, json
t, r25, r34, r35, w = sp.symbols('t r25 r34 r35 w')
d = pickle.load(open('kap2.pkl', 'rb'))
vals = d['vals']
a = [sp.Rational(x) for x in sys.argv[1:4]]
sub = {r25: a[0], r34: a[1], r35: a[2]}
tv = [(sp.cancel(x.subs(sub)), sp.cancel(y.subs(sub))) for x, y in vals]
s = 0
for k, (tk, yk) in enumerate(tv):
    den = 1
    for j, (tj, _) in enumerate(tv):
        if j != k: den *= (tk - tj)
    s += yk/den
n = sp.fraction(sp.factor(sp.together(s)))[0]
facs = sp.factor_list(n)[1]
for f, m in facs:
    print('deg', sp.degree(f, w), 'mult', m)
f = max(facs, key=lambda fm: sp.degree(fm[0], w))[0]
P = sp.Poly(f, w)
c = [int(x) for x in P.all_coeffs()]
print(json.dumps(c))
json.dump(dict(point=[str(x) for x in a], coeffs_leading_first=c), open('minpoly_%s_%s_%s.json' % tuple(x.replace('/', 'o') for x in sys.argv[1:4]), 'w'))
