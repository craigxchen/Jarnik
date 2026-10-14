"""Prop 2.3 check: the restriction of L_8 (width fat-point system of M8) to a Fano span U_l is zero.
psi(t1,t2,t3) = Phi(t1 p1 + t2 p2 + t3 p3): conditions D_v^{tau+1} psi = 0 for v = (a1,a2,a3) in the p-basis with tau_A(v) < 8."""
import itertools, pickle, math
import numpy as np
import sympy as sp
from m8lib import *
A, cols, H, JA, W = pickle.load(open('m8_step1.pkl','rb'))
Aa = np.array(A)
print("all k' with |k'|_1 <= 2 in A:", all(tuple(k) in set(A) for k in itertools.product((-1,0,1), repeat=7) if sum(map(abs,k)) <= 2))
a, b = 1, 2
Ta = set(i for i in range(7) if cols[a][i]); Tb = set(i for i in range(7) if cols[b][i])
ind = lambda S: tuple(1 if i in S else 0 for i in range(7))
p = [ind(Ta & Tb), ind(Ta - Tb), ind(Tb - Ta)]
def tau_h(v):
    vals = Aa @ np.array(v); return len(set(vals.tolist())) - 1, int(vals.max() - vals.min())
for name, cf in [("p1",(1,0,0)),("p2",(0,1,0)),("p3",(0,0,1)),("p1+p2",(1,1,0)),("p1+p3",(1,0,1)),("p2+p3",(0,1,1)),("z=p1+p2+p3",(1,1,1))]:
    v = tuple(sum(c*pp[i] for c, pp in zip(cf, p)) for i in range(7))
    print(name, v, "(tau,h) =", tau_h(v))
T = sp.symbols('t1 t2 t3')
mons = [e for e in itertools.product(range(9), repeat=3) if sum(e) == 8]
C = sp.symbols('c0:%d' % len(mons))
psi = sum(c * T[0]**e[0] * T[1]**e[1] * T[2]**e[2] for c, e in zip(C, mons))
def conds(points):
    eqs = []
    for cf in points:
        v = tuple(sum(c*pp[i] for c, pp in zip(cf, p)) for i in range(7))
        tau, h = tau_h(v)
        if tau >= 8: continue
        Dv = lambda f: sum(cf[i] * sp.diff(f, T[i]) for i in range(3))
        f = psi
        for _ in range(tau + 1): f = Dv(f)
        f = sp.Poly(sp.expand(f), *T)
        eqs += list(f.coeffs())
    M = sp.Matrix([[sp.Poly(e, *C).coeff_monomial(c) for c in C] for e in eqs])
    return len(mons) - M.rank()
seven = [(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1),(1,1,1)]
print("dim of ternary octics with the multiplicities at the 7 points P_i, Q_ij, Z:", conds(seven))
allpts = [v for v in itertools.product(range(-2,3), repeat=3) if any(v) and [t for t in v if t][0] > 0 and math.gcd(*[abs(t) for t in v if t]) == 1]
print("dim with all integer points (coefficients in [-2,2]) of U_l with tau < 8:", conds(allpts))
