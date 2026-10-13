#!/usr/bin/env python3
"""Referee check of the claim (p3lower.md, after Cor. A.1 and Sec. 0 item 4):
 'the strongest consequence of "W I - 2 mu is PSD and mu >= D - c_C" is exactly the Perron bound'.
Counterexample when c_C > 0 (entries of B = D - c_C (J-I) negative), M = 5, via the Horn matrix.
 D = 2 * adjacency of the 5-cycle, c_C = 1, so B = D - (J - I) = I - Horn.
Exact arithmetic in Q(sqrt5) for the PSD certificate; LP for W_cut."""
import itertools, math
from fractions import Fraction as F
import sympy as sp
import numpy as np
from scipy.optimize import linprog

M = 5
D = [[0] * M for _ in range(M)]
for i in range(M):
    D[i][(i + 1) % M] = D[(i + 1) % M][i] = 2
cC = 1
B = [[(D[x][y] - cC) if x != y else 0 for y in range(M)] for x in range(M)]
Horn = [[(1 if x == y else -B[x][y]) for y in range(M)] for x in range(M)]
print("B =", B)

# 1) Perron-type bound sup_{u>=0,|u|=1} u^T B u = 1 - min u^T Horn u = 1 (Horn copositive).
#    Exact: u^T Horn u >= 0 on u >= 0 (Hall-Newman). Check by fine grid on the simplex (evidence) and
#    exact attainment u = (e1+e2)/sqrt2 -> u^T B u = 1.
Bs = sp.Matrix(B)
u = sp.Matrix([1, 1, 0, 0, 0]) / sp.sqrt(2)
print("u^T B u at (e1+e2)/sqrt2 =", sp.nsimplify((u.T * Bs * u)[0]))
best = -1e9
N = 24
for c in itertools.product(range(N + 1), repeat=4):
    if sum(c) > N: continue
    v = np.array(list(c) + [N - sum(c)], dtype=float)
    v /= np.linalg.norm(v)
    val = v @ np.array(B, dtype=float) @ v
    best = max(best, val)
print("grid max of u^T B u over u>=0 (evidence):", round(best, 6))

# Cor A.1 formula: 2 lmax(D) - 2 (M-1) c_C^+
lmD = max(np.linalg.eigvalsh(np.array(D, dtype=float)))
print("Cor A.1 bound 2 lmax(D) - 2(M-1)c_C =", 2 * lmD - 2 * (M - 1) * cC)

# 2) PSD consequence: for Z doubly nonnegative with tr Z = 1, PSD(W I - 2 mu) and mu >= B give
#    W >= 2 <Z, mu> >= 2 <Z, B>.  Z = circ(1, 3/5, 0, 0, 3/5)/5.
a = sp.Rational(3, 5)
Z = sp.zeros(M, M)
for i in range(M):
    Z[i, i] = sp.Rational(1, 5)
    Z[i, (i + 1) % M] = Z[(i + 1) % M, i] = a / 5
eigs = [sp.nsimplify(sp.Rational(1, 5) * (1 + 2 * a * sp.cos(2 * sp.pi * k / 5))) for k in range(5)]
eigs = [sp.radsimp(sp.simplify(e)) for e in eigs]
print("eigenvalues of Z (exact):", eigs, " all >= 0:", all(sp.N(e) >= 0 for e in eigs))
print("min eigenvalue numeric:", min(float(sp.N(e)) for e in eigs))
inner = sum(Z[x, y] * B[x][y] for x in range(M) for y in range(M) if x != y)
print("<Z,B> =", inner, " => PSD consequence W >= 2<Z,B> =", 2 * inner)

# 3) exact W_cut(B) via LP
S = []
for bits in itertools.product([1, -1], repeat=M - 1):
    if all(b == 1 for b in bits): continue
    S.append((1,) + bits)
pairs = [(x, y) for x in range(M) for y in range(x + 1, M)]
A = np.array([[0.5 if s[x] != s[y] else -0.5 for s in S] for (x, y) in pairs])
b = np.array([B[x][y] for (x, y) in pairs], dtype=float)
res = linprog(np.ones(len(S)), A_ub=-A, b_ub=-b, bounds=(0, None), method="highs")
print("W_cut(B) =", res.fun)
ok = (2 * inner > 2 * 1) and all(sp.N(e) >= 0 for e in eigs) and best <= 1 + 1e-9
print("COUNTEREXAMPLE CONFIRMED: PSD bound (>= %s) exceeds Perron-type bound 2*sup u^TBu = 2" % (2 * inner) if ok else "NOT CONFIRMED")
