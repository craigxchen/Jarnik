"""cub13 family in the Kapranov psi_7 model (exact, sympy).
P^4, p1..p5 = e0..e4, p6 = (1,1,1,1,1); curve c(t) = (h, X2, X3, X4, X5), cubics in t.
Gauge: u = t_126 = infinity  (so lam3 = lam4 = lam5 = 1), r45 = 0, r24 = 1 (affine normalisation).
Free: r25, r34, r35, w (= t_456).  r23 solved (Moebius), v = t_136 rational.
Output: the final condition DD (4th divided difference), factored."""
import sympy as sp
t, r25, r34, r35, w = sp.symbols('t r25 r34 r35 w')
r24, r45 = sp.Integer(1), sp.Integer(0)
r23 = sp.Symbol('r23')
lam2 = sp.Symbol('lam2')
# u = infinity: lam3 = lam4 = lam5 = 1 (leading coefficients equal)
# v: X2(v)=X4(v)=X5(v); X4=X5 at v: (v-r24)(v-r34) = (v-r25)(v-r35)
v = (r25*r35 - r24*r34)/(r25 + r35 - r24 - r34)
X3 = (t - r23)*(t - r34)*(t - r35)
X4 = (t - r24)*(t - r34)*(t - r45)
X5 = (t - r25)*(t - r35)*(t - r45)
X2 = lam2*(t - r23)*(t - r24)*(t - r25)
# lam2 from X2(v) = X4(v)
lam2v = sp.cancel(X4.subs(t, v) / ((v - r23)*(v - r24)*(v - r25)))
# w: X2(w) = X3(w)  ->  lam2 (w-r24)(w-r25) = (w-r34)(w-r35)   [r23 cancels]
eq = sp.together(lam2v*(w - r24)*(w - r25) - (w - r34)*(w - r35))
num = sp.factor(sp.numer(eq))
print('r23 equation numerator factored:', num)
sol = sp.solve(sp.numer(eq), r23)
print('r23 solutions:', len(sol))
r23s = sp.factor(sol[0])
print('r23 =', r23s)
