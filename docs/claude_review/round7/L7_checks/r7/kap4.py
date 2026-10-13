"""Full multivariate factorisation of the DD numerator (exact, sympy)."""
import sympy as sp, pickle, time
t, r25, r34, r35, w = sp.symbols('t r25 r34 r35 w')
d = pickle.load(open('kap2.pkl', 'rb'))
vals = d['vals']
T0 = time.time()
tv = vals
# common-denominator approach: DD = sum_k y_k / prod_{j!=k} (t_k - t_j)
s = 0
for k, (tk, yk) in enumerate(tv):
    den = 1
    for j, (tj, _) in enumerate(tv):
        if j != k: den *= sp.factor(sp.together(tk - tj))
    s += yk/den
print('summing', time.time()-T0, flush=True)
s = sp.together(s)
n, dn = sp.fraction(s)
print('together', time.time()-T0, flush=True)
fl = sp.factor_list(sp.expand(n))
print('factored', time.time()-T0, flush=True)
facs = []
for f, m in fl[1]:
    degs = [sp.degree(f, x) for x in (r25, r34, r35, w)]
    print('factor degs (r25,r34,r35,w)=', degs, 'mult', m, 'nterms', len(sp.Poly(f, r25, r34, r35, w).terms()))
    facs.append((f, m))
pickle.dump(dict(const=fl[0], facs=facs), open('kap4.pkl', 'wb'))
