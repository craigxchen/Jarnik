import sympy as sp
from check2 import family_sym, cayley
t, L = sp.symbols('t lam')
# eigen identity: omega=(lam-1)+i t, lam^2-(t^2+2)lam+1=0  =>  omega^2 = t(lam-1)(t+2i)
w2 = sp.expand(((L-1)+sp.I*t)**2 - t*(L-1)*(t+2*sp.I))
print('omega^2 - t(lam-1)beta mod minpoly:', sp.rem(w2, L**2-(t**2+2)*L+1, L))
# check omega is eigenvector of A with eigenvalue lam (i.e. (lam-1, t))
A = sp.Matrix([[t**2+1, t],[t,1]]); v = sp.Matrix([L-1, t])
print('A v - lam v mod minpoly:', [sp.rem(sp.expand(e), L**2-(t**2+2)*L+1, L) for e in (A*v-L*v)])
S = sp.symbols('s1:5')
for n in (1,2,3):
    marks = family_sym(n, list(S), sp.Integer(2))
    us = [cayley(m, marks[0]) for m in marks[1:]]
    Jm = sp.Matrix([[sp.diff(u, s) for s in S] for u in us])
    pt = {S[0]: 1, S[1]: 2, S[2]: 3, S[3]: 5}
    d = sp.nsimplify(Jm.subs(pt).det())
    print('T=2 n=%d det='%n, sp.N(d,6), d!=0)
