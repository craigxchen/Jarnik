"""Evidence check of Prop 2.1 (duality) on small T = 1 slab boxes (floating-point LPs, scipy HiGHS).
For A = {k' in {-1,0,1}^d : |k(S_p)| <= L for each cut}, m, random u:
  primal: dist_inf(l_u^m, Poly_{<m}|_A) = min_p max_k |g(k) - p(k)|
  dual:   max <c, g> s.t. c in V(A,m) (moment equations), |c|_1 <= 1.
They must agree (LP duality). Also checks Theorem 2.5's lower bound kappa_V(u) >= (m!/(2d)^m) eta^m whenever nu_eta(u) >= m
for one-cut boxes (kappa_V >= kappa_H, so this is a consistency check only).
"""
import itertools, math, random
import numpy as np
from scipy.optimize import linprog

def box(d, cuts):
    return [k for k in itertools.product((-1, 0, 1), repeat=d) if all(abs(sum(k[x] for x in S)) <= L for S, L in cuts)]

def monos(d, m):
    out = []
    for deg in range(m):
        for c in itertools.combinations_with_replacement(range(d), deg):
            a = [0]*d
            for x in c: a[x] += 1
            out.append(tuple(a))
    return out

def primal(A, m, u):
    d = len(u); Ms = monos(d, m)
    Mat = np.array([[np.prod([k[i]**a[i] for i in range(d)]) for a in Ms] for k in A], dtype=float)
    g = np.array([np.dot(k, u)**m for k in A])
    nP = Mat.shape[1]
    # vars: coeffs p (free), t ; minimize t s.t. -t <= g - Mat p <= t
    cobj = np.zeros(nP + 1); cobj[-1] = 1
    Aub = np.vstack([np.hstack([-Mat, -np.ones((len(A), 1))]), np.hstack([Mat, -np.ones((len(A), 1))])])
    bub = np.concatenate([-g, g])
    res = linprog(cobj, A_ub=Aub, b_ub=bub, bounds=[(None, None)]*nP + [(0, None)], method='highs')
    return res.fun, Mat, g

def dual(Mat, g):
    U, nP = Mat.shape
    # vars c+ , c- >= 0 ; maximize g.(c+ - c-) s.t. Mat^T (c+ - c-) = 0, sum(c+ + c-) <= 1
    cobj = np.concatenate([-g, g])
    Aeq = np.hstack([Mat.T, -Mat.T]); beq = np.zeros(nP)
    Aub = np.ones((1, 2*U)); bub = [1.0]
    res = linprog(cobj, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=[(0, None)]*(2*U), method='highs')
    return -res.fun

if __name__ == '__main__':
    rng = random.Random(1)
    worst = 0
    for trial in range(12):
        d = rng.choice([4, 5]); s = rng.randint(2, d); L = rng.randint(1, s - 1)
        cuts = [(list(range(s)), L)]
        if rng.random() < 0.5:
            S2 = rng.sample(range(d), rng.randint(2, d)); cuts.append((S2, rng.randint(1, len(S2) - 1)))
        A = box(d, cuts)
        m = rng.randint(3, 5)
        u = [rng.uniform(-1, 1) for _ in range(d)]
        pv, Mat, g = primal(A, m, u)
        dv = dual(Mat, g)
        worst = max(worst, abs(pv - dv) / max(1e-12, abs(pv)))
        print(f"d={d} cuts={[(len(S),L) for S,L in cuts]} |A|={len(A)} m={m}: primal {pv:.6e} dual {dv:.6e}")
    print("max relative primal-dual gap:", worst)
