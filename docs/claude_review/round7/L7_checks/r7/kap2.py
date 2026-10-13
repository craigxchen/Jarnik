"""continue kap1: contact parameters and final divided-difference condition (exact, sympy)."""
import sympy as sp, time, pickle
t, r25, r34, r35, w = sp.symbols('t r25 r34 r35 w')
r24, r45 = sp.Integer(1), sp.Integer(0)
v = (r25*r35 - r24*r34)/(r25 + r35 - r24 - r34)
r23 = (r25*r35 - r34)*(r25*r34*r35 - r25*r34 - 2*r25*r35*w + r25*r35 + r25*w**2 - r34*r35 - r34*w**2 + 2*r34*w + r35*w**2 - w**2)/((r25 - 1)*(r34 - w)*(r35 - w)*(r25 - r34 + r35 - 1))
r23 = sp.factor(r23)
X4 = (t - r24)*(t - r34)*(t - r45)
X5 = (t - r25)*(t - r35)*(t - r45)
lam2 = sp.factor(X4.subs(t, v) / ((v - r23)*(v - r24)*(v - r25)))
X2 = lam2*(t - r23)*(t - r24)*(t - r25)
X3 = (t - r23)*(t - r34)*(t - r35)
def other_root(P, known):
    P = sp.Poly(sp.expand(P), t)
    c = P.all_coeffs()  # leading first
    d = P.degree()
    # sum of roots = -c1/c0
    return sp.factor(-c[1]/c[0] - sum(known))
T0 = time.time()
# X3 - X5 is quadratic (u = inf), roots r35 and t246
t246 = other_root(sp.expand(X3 - X5), [r35])
t256 = other_root(sp.expand(X3 - X4), [r34])
print('t246 =', t246); print('t256 =', t256, flush=True)
P25 = sp.together(X2 - X5); P25n = sp.numer(P25)
t346 = other_root(P25n, [r25, v])
P24 = sp.together(X2 - X4); P24n = sp.numer(P24)
t356 = other_root(P24n, [r24, v])
print('t346 =', t346); print('t356 =', t356, flush=True)
pts = [(t246, X3), (t256, X3), (t346, X2), (t356, X2), (w, X2)]
vals = [(sp.factor(tk), sp.factor(Xk.subs(t, tk))) for tk, Xk in pts]
print('values done', time.time() - T0, flush=True)
pickle.dump(dict(vals=vals, r23=r23, lam2=lam2, v=v), open('kap2.pkl', 'wb'))
