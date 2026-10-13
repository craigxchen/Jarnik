"""LP upper bound for zeta(A) (rigorous via a rational dual certificate).
Every zonotope k0 + sum [0,L_i] v_i in A satisfies, for each constraint functional lam with |lam| <= b on A,
sum_i L_i |lam(v_i)| <= 2b.  Dual: if y_lam >= 0 and sum_lam y_lam |lam(v)| >= 1 for every candidate primitive v
(|v_x| <= 2T), then sum L_i <= sum_lam 2 b_lam y_lam.  We solve the LP, rationalise y, and verify exactly."""
import sys, itertools, math, ast
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog
from zono import functionals, prim_vectors, dot

def lp_bound(n, T, prof, L, d=None, lams=None):
    if lams is None:
        d = n - 1
        lams = functionals(n, T, prof, L)
    vecs = prim_vectors(d, 2 * T)
    caps = [2 * b for l, b in lams]
    vecs = [v for v in vecs if all(abs(dot(l, v)) <= c for (l, b), c in zip(lams, caps))]
    Mv = np.array([[abs(dot(l, v)) for l, b in lams] for v in vecs], dtype=float)
    # dual LP: min sum caps*y s.t. Mv y >= 1, y >= 0
    res = linprog(np.array(caps, float), A_ub=-Mv, b_ub=-np.ones(len(vecs)), bounds=[(0, None)] * len(lams), method='highs')
    y = [Fraction(x).limit_denominator(10 ** 6) for x in res.x]
    # make feasible exactly: scale up by the worst violation
    worst = min(sum(Fraction(abs(dot(l, v))) * yy for (l, b), yy in zip(lams, y)) for v in vecs)
    assert worst > 0
    y = [yy / worst for yy in y] if worst < 1 else y
    assert all(sum(Fraction(abs(dot(l, v))) * yy for (l, b), yy in zip(lams, y)) >= 1 for v in vecs)
    return sum(Fraction(c) * yy for c, yy in zip(caps, y)), res.fun

if __name__ == "__main__":
    n = int(sys.argv[1]); T = int(sys.argv[2]); L = int(sys.argv[3]); prof = ast.literal_eval(sys.argv[4])
    b, f = lp_bound(n, T, prof, L)
    print("LP bound zeta <= %s = %.4f (LP value %.4f)" % (b, float(b), f))
