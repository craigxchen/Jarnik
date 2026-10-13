"""Exact: genuine factor of DD restricted to a 2-plane (two of r25,r34,r35,w fixed to rationals).
usage: python3 slice2.py fixedvar1 val1 fixedvar2 val2  -> pickles bivariate genuine factor."""
import sympy as sp, pickle, sys
t, r25, r34, r35, w = sp.symbols('t r25 r34 r35 w')
S = dict(r25=r25, r34=r34, r35=r35, w=w)
d = pickle.load(open('kap2.pkl', 'rb'))
vals = d['vals']
def DD(sub):
    tv = [(sp.cancel(a.subs(sub)), sp.cancel(b.subs(sub))) for a, b in vals]
    s = 0
    for k, (tk, yk) in enumerate(tv):
        den = 1
        for j, (tj, _) in enumerate(tv):
            if j != k: den *= sp.factor(tk - tj)
        s += yk/den
    return sp.together(s)
v1, a1, v2, a2 = sys.argv[1], sp.Rational(sys.argv[2]), sys.argv[3], sp.Rational(sys.argv[4])
sub = {S[v1]: a1, S[v2]: a2}
e = DD(sub)
n, dn = sp.fraction(e)
fl = sp.factor_list(sp.expand(n))
free = [S[k] for k in S if k not in (v1, v2)]
best = max(fl[1], key=lambda fm: sp.Poly(fm[0], *free).total_degree())
for f, m in fl[1]:
    print('factor degrees', [sp.degree(f, x) for x in free], 'total', sp.Poly(f, *free).total_degree(), 'mult', m)
G = best[0]
pickle.dump(dict(G=G, free=[str(x) for x in free], sub={v1: str(a1), v2: str(a2)}), open('slice_%s_%s_%s_%s.pkl' % (v1, sys.argv[2].replace('/', 'o'), v2, sys.argv[4].replace('/', 'o')), 'wb'))
print('genuine candidate: degrees', [sp.degree(G, x) for x in free], 'nterms', len(sp.Poly(G, *free).terms()))
