"""Independent recomputation of the II.3 table (Walsh model on F_2^5, 31 cuts, e_p=1, equal weights).
Assignment optima by LP (HiGHS, assignment polytope is integral), cross-checked by brute force for d=1
and by scipy LSA. Collected exponent by the factorization (min of a linear functional over a Minkowski
sum is additive), cross-checked for d=1 against the support of the expanded determinant."""
from fractions import Fraction
from itertools import permutations, product
from math import comb
import numpy as np
from scipy.optimize import linprog, linear_sum_assignment
import sympy as sp

def total_degree(k, d):
    if k == 0:
        yield (); return
    for a in range(d + 1):
        for t in total_degree(k - 1, d - a):
            yield (a,) + t

def lp_assign(C):
    n = C.shape[0]
    Aeq = np.zeros((2 * n, n * n))
    for i in range(n):
        Aeq[i, i * n:(i + 1) * n] = 1
        Aeq[n + i, i::n] = 1
    r = linprog(C.ravel(), A_eq=Aeq, b_eq=np.ones(2 * n), bounds=(0, 1), method='highs')
    val = round(r.fun)
    assert abs(r.fun - val) < 1e-7
    r2, c2 = linear_sum_assignment(C)
    assert C[r2, c2].sum() == val
    return val

def walsh():
    return [[bin(x & v).count('1') % 2 for x in range(32)] for v in range(1, 32)]

def build(d):
    A = list(total_degree(5, d))
    node = lambda h, j: 5 * j + h
    sp_ = [(d, 0, 0, 0, 0), (d - 1, 1, 0, 0, 0)]
    rows, fac = [], []
    for b in A:
        if b in sp_: continue
        rows.append(tuple(node(h, b[h]) for h in range(5)))
        for h in range(5):
            for j in range(b[h]):
                fac.append((node(h, b[h]), node(h, j)))
    u, v = 5 * (d + 1), 5 * (d + 1) + 1
    tail = tuple(node(h, 0) for h in (2, 3, 4))
    rows += [(u, v) + tail, (v, u) + tail]
    for j in range(d - 1):
        fac += [(u, node(0, j)), (v, node(0, j))]
    fac += [(u, v), (u, v, node(0, d - 1), node(1, 0))]
    return A, rows, fac

def raw(A, rows, bits):
    E = np.array(A)
    Val = np.array([[bits[j] for j in r] for r in rows])
    C = Val @ E.T
    V = int(E.sum())
    return V + lp_assign(C) - (-lp_assign(-C))

def collected(fac, bits):
    return sum(1 for f in fac if len({bits[j] for j in f}) == 1)

cuts = walsh()
print('d  n   V   raw kappa   collected kappa   fair-cut collected ((V-1)/2+1/8)/V')
for d in (1, 2, 3, 4):
    A, rows, fac = build(d)
    V = sum(map(sum, A)); assert V == len(fac) == 5 * comb(d + 5, 6)
    R = sum(raw(A, rows, b) for b in cuts)
    Cl = sum(collected(fac, b) for b in cuts)
    print(d, len(A), V, Fraction(R, 31 * V), Fraction(Cl, 31 * V), Fraction(2 * (V - 1) + 1, 8 * V) * 2 / 2 if False else Fraction(4 * (V - 1) + 1, 8 * V))

# d = 1 brute-force check of raw and of collected vs expanded support
A, rows, fac = build(1)
x = sp.symbols('x0:32')
Mx = sp.Matrix([[sp.prod([x[r[h]] ** al[h] for h in range(5)]) for al in A] for r in rows])
F = sp.Poly(sp.expand(Mx.det()), *x)
supp = [m for m, c in F.terms()]
tot_raw = tot_true = 0
for bits in cuts:
    E = np.array(A); Val = np.array([[bits[j] for j in r] for r in rows]); C = Val @ E.T
    costs = [sum(C[i, p[i]] for i in range(6)) for p in permutations(range(6))]
    tot_raw += 5 + min(costs) - max(costs)
    s = [sum(m[j] * bits[j] for j in range(32)) for m in supp]
    tot_true += 5 + min(s) - max(s)
print('d=1 brute force raw kappa', Fraction(tot_raw, 155), ' true-support kappa', Fraction(tot_true, 155))

# cyclic design: all 64 patterns, and Walsh labels 0,1,2,4,8,16, and random labels
z = sp.symbols('z0:6')
Mc = sp.Matrix([[1] + [z[(i + h) % 6] for h in range(5)] for i in range(6)])
Fc = sp.Poly(sp.expand(Mc.det()), *z)
suppc = [m for m, c in Fc.terms()]
ok = True
for pat in product((0, 1), repeat=6):
    s = [sum(m[j] * pat[j] for j in range(6)) for m in suppc]
    val = 5 + min(s) - max(s)
    ok &= (val == (5 if len(set(pat)) == 1 else 0))
print('cyclic: true-support exponent = 5*[constant pattern] on all 64 patterns:', ok)
import random
random.seed(3)
vals = []
for t in range(200):
    lab = random.sample(range(32), 6)
    tot = 0
    for bits in cuts:
        pat = [bits[j] for j in lab]
        tot += 5 if len(set(pat)) == 1 else 0
    vals.append(Fraction(tot, 155))
print('cyclic on 200 random 6-label sets in Walsh model: max kappa', max(vals), ' (=', float(max(vals)), ')')
