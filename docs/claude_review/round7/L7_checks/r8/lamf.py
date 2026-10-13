# lambda_{f_i}(gamma) = -psi_i.gamma + sum_j D_ij.gamma, computed from D_S-degree vectors (exact).
import itertools
from fractions import Fraction
LAB = ['O','a1','a2','b1','b2','c1','c2']; n = 7
def key(T):
    T = frozenset(T); C = frozenset(range(n)) - T
    return min(tuple(sorted(T)), tuple(sorted(C)), key=lambda s: (len(s), s)) if len(T) != len(C) else min(tuple(sorted(T)), tuple(sorted(C)))
def split_key(T):
    T = frozenset(T); C = frozenset(range(n)) - T
    return frozenset([T, C])
def psi_dot(i, cls):
    j, k = [a for a in range(n) if a != i][:2]
    tot = 0
    others = [a for a in range(n) if a not in (i, j, k)]
    for r in range(0, len(others) + 1):
        for S in itertools.combinations(others, r):
            T = (i,) + S
            if 2 <= len(T) <= n - 2: tot += cls.get(split_key(T), 0)
    return tot
def lam_f(i, cls):
    return -psi_dot(i, cls) + sum(cls.get(split_key((i, j)), 0) for j in range(n) if j != i)
idx = {l: a for a, l in enumerate(LAB)}
gA = {('O','a1','a2'): 2, ('O','b1','b2'): 2, ('O','c1','c2'): 2, ('a1','b1','c1'): 1, ('a1','b1','c2'): 1,
      ('a1','b2','c1'): 1, ('a1','b2','c2'): 1, ('a2','b1','c1'): 1, ('a2','b2','c1'): 1, ('a2','b1'): 1,
      ('a2','b2'): 1, ('a2','c2'): 2, ('b1','c2'): 1, ('b2','c2'): 1}
cls = {split_key([idx[l] for l in T]): v for T, v in gA.items()}
beta = {split_key(T): 1 for r in (2, 3) for T in itertools.combinations(range(n), r)}
print('check psi.beta =', psi_dot(0, beta), ' lambda_f(beta) =', lam_f(0, beta))
f3 = {split_key((3, j)): 1 for j in range(n) if j != 3}
print('check lambda_{f_O}(f_b1) =', lam_f(0, f3), ' lambda_{f_b1}(f_b1) =', lam_f(3, f3))
vals = {LAB[i]: lam_f(i, cls) for i in range(n)}
print('lambda_{f_i}(gamma_A) by label i:', vals, ' psi_i.gamma_A:', {LAB[i]: psi_dot(i, cls) for i in range(n)})
mu = min(vals.values()); print('min =', mu, ' (35/12)*min =', Fraction(35, 12)*mu, ' vs -9')
