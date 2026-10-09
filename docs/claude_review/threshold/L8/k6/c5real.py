"""NUMERICAL check that the k=6 C5 ansatz has nondegenerate (real) balanced cubic curves:
random rational (rho,sigma) with rho12 free -> roots of the cubic -> c,tau,z by nullspace -> check all
25 boundary contacts are simple and distinct (pair roots computed numerically)."""
import sys, random, itertools
sys.path.insert(0, '..')
import numpy as np
from fractions import Fraction as Fr
from polylib import *
from c5ansatz import LABELS, EDGES, TRIPLES, nbr_edges, triples_of, matrix
from c5deg import interp

random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
found = 0
for trial in range(30):
    vals = random.sample(range(-30, 31), 9)
    rho = {e: Fr(v) for e, v in zip(EDGES[1:], vals[:4])}
    sigma = {T: Fr(v) for T, v in zip(TRIPLES, vals[4:])}
    def N(x):
        r = dict(rho); r[(1, 2)] = x
        D = prod([linpoly(-1, sigma[T]) for T in triples_of(1)] + [linpoly(-1, sigma[T]) for T in triples_of(2)])
        return det_frac(matrix(r, sigma)) * ev(D, x)
    xs = [Fr(10 ** 6 + 17 * k) for k in range(5)]
    p = interp(xs, [N(x) for x in xs])
    q, rem = divmod_(p, linpoly(1, -sigma[(1, 2, 4)]))
    roots = np.roots([float(c) for c in reversed(q)])
    for r12 in roots:
        if abs(r12.imag) > 1e-9: continue
        r12 = r12.real
        rh = {e: float(v) for e, v in rho.items()}; rh[(1, 2)] = r12
        sg = {T: float(v) for T, v in sigma.items()}
        M = np.array([[float(0)] * 5 for _ in range(5)])
        for ii, i in enumerate(LABELS):
            e1, e2 = nbr_edges(i); Ts = triples_of(i)
            for kk in range(3):
                T1, T2, T3 = Ts[kk], Ts[(kk + 1) % 3], Ts[(kk + 2) % 3]
                M[ii, TRIPLES.index(T1)] = (sg[T2] - sg[T3]) / ((sg[T1] - rh[e1]) * (sg[T1] - rh[e2]))
        U, S, Vt = np.linalg.svd(M)
        z = dict(zip(TRIPLES, Vt[-1]))
        u = {}
        okc = True
        for i in LABELS:
            e1, e2 = nbr_edges(i); Ts = triples_of(i)
            A = {T: (sg[T] - rh[e1]) * (sg[T] - rh[e2]) for T in Ts}
            # z_T / A = c (sigma_T - tau): fit line through 3 points
            X = np.array([[sg[T], -1.0] for T in Ts]); Y = np.array([z[T] / A[T] for T in Ts])
            sol, res, _, _ = np.linalg.lstsq(X, Y, rcond=None)
            c, ct = sol
            if abs(c) < 1e-9: okc = False; break
            tau = ct / c
            u[i] = np.poly1d([c]) * np.poly1d([1, -rh[e1]]) * np.poly1d([1, -rh[e2]]) * np.poly1d([1, -tau])
        if not okc: continue
        # contacts: roots of u_i (5x3 = 15, contains 5 shared rho -> 10 distinct) and u_i - u_j (10 pairs x 3)
        special = []
        for i in LABELS:
            special += list(np.roots(u[i].coeffs))
        pairroots = []
        for i, j in itertools.combinations(LABELS, 2):
            pairroots += list(np.roots((u[i] - u[j]).coeffs))
        allpts = special + pairroots
        # cluster with tolerance; expected 25 distinct contact parameters
        cl = []
        for x in allpts:
            if not any(abs(x - y) < 1e-6 for y in cl): cl.append(x)
        lead = [u[i].coeffs[0] for i in LABELS]
        mlead = min(abs(a - b) for a, b in itertools.combinations(lead, 2))
        if len(cl) == 25 and mlead > 1e-6:
            found += 1
print('k=6 C5 ansatz: real nondegenerate balanced cubic curves found:', found)
