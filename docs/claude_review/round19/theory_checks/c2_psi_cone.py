"""c2: beta is a nonnegative combination of psi-monomial curve classes psi^a (|a| = 4) on Mbar_{0,8}.
Each psi_i is nef (Kapranov: |psi_i| is base-point-free), so every psi^a is nonnegative on every
effective divisor; hence so is beta.  Independent of Keel-McKernan's symmetric effective cone theorem.
D_T . psi^a = int_{Mbar_{0,|T|+1}} prod_{i in T} psi_i^{a_i} * int_{Mbar_{0,|T^c|+1}} prod_{i in T^c} psi_i^{a_i},
with int_{Mbar_{0,m}} prod psi^{b} = (m-3)!/prod b! if sum b = m-3, else 0.
Exact arithmetic.  We symmetrize over S_8: s_lambda = sum of psi^a over a of exponent type lambda.
"""
from fractions import Fraction as Fr
from itertools import combinations, product
from math import factorial
import sys
def splits(n):
    S = []
    for s in range(2, n//2+1):
        for T in combinations(range(n), s):
            T = frozenset(T); Tc = frozenset(range(n)) - T
            if s == n-s and min(Tc) < min(T): continue
            S.append(T)
    return S
DIV = splits(8)
def psi_int(bs, m):
    if sum(bs) != m-3: return 0
    r = factorial(m-3)
    for b in bs: r //= 1 if b == 0 else 1  # placeholder
    num = factorial(m-3); den = 1
    for b in bs: den *= factorial(b)
    assert num % den == 0
    return num // den
def DT_dot(T, a):
    Tc = [i for i in range(8) if i not in T]
    return psi_int([a[i] for i in T], len(T)+1) * psi_int([a[i] for i in Tc], len(Tc)+1)
# all exponent vectors with |a| = 4
A = [a for a in product(range(5), repeat=8) if sum(a) == 4]
def typ(a): return tuple(sorted([x for x in a if x > 0], reverse=True))
types = sorted(set(typ(a) for a in A), reverse=True)
print("exponent types:", types, " #monomials:", len(A))
# symmetric pairings: for each type and each split size s in {2,3,4}, D_T.s_lambda (same for all T of that size)
M = {}
for lam in types:
    row = {}
    for T in DIV:
        v = sum(DT_dot(T, a) for a in A if typ(a) == lam)
        s = len(T)
        if s in row: assert row[s] == v, (lam, s)
        row[s] = v
    M[lam] = row
    print("type", lam, ": D_T.s_lambda for |T|=2,3,4 :", row[2], row[3], row[4])
# also list individual monomials psi^a: check 'beta = sum c_a psi^a' via symmetric types; solve LP by enumeration of
# 3-subsets of types (basic solutions), exact.
import itertools
sols = []
for basis in itertools.combinations(types, 3):
    # solve 3x3 system  sum_l c_l M[l][s] = 1, s = 2,3,4
    mat = [[Fr(M[l][s]) for l in basis] for s in (2,3,4)]
    rhs = [Fr(1)]*3
    # gaussian elimination
    m = [row[:] + [r] for row, r in zip(mat, rhs)]
    ok = True
    for c in range(3):
        piv = next((r for r in range(c, 3) if m[r][c] != 0), None)
        if piv is None: ok = False; break
        m[c], m[piv] = m[piv], m[c]
        for r in range(3):
            if r != c and m[r][c] != 0:
                f = m[r][c] / m[c][c]
                m[r] = [x - f*y for x, y in zip(m[r], m[c])]
    if not ok: continue
    c = [m[i][3]/m[i][i] for i in range(3)]
    if all(x >= 0 for x in c):
        sols.append((basis, c))
print("nonnegative basic solutions of beta = sum c_lambda s_lambda:")
for b, c in sols: print("  ", dict(zip(b, [str(x) for x in c])))
if sols:
    b, c = sols[0]
    # verify on every one of the 119 divisors with the full (unsymmetrized) sum
    for T in DIV:
        v = sum(cl * sum(DT_dot(T, a) for a in A if typ(a) == l) for l, cl in zip(b, c))
        assert v == 1, (T, v)
    print("verified: D_T.(sum c s_lambda) = 1 for all 119 D_T (exact)")
else:
    print("NO nonnegative solution with psi-monomials"); sys.exit(1)
