"""Exact checks for the canonical-class computation used in Sections 4 and 6 of ptolemy.md.

(K1) In Kapranov's model (Mbar_{0,n} = iterated blow-up of P^{n-3} at q_1..q_{n-1} and the linear
     spans of subsets I with 1 <= |I| <= n-4, in increasing dimension), the canonical class is
        K = -(n-2) H + sum_I (n-3-|I|) E_I .
     The boundary divisors are D_{I+n} = E_I (1<=|I|<=n-4) and
        D_{I+n} = H - sum_{0 != K subset I, K != I} E_K   (|I| = n-3).
     We verify the Keel-McKernan/Pandharipande formula
        K = sum_S ( s(n-s)/(n-1) - 2 ) D_S          (s = |S|, S up to complement)
     coefficient-by-coefficient in this basis, for n = 5..10, with exact fractions.
(K2) For the 'balanced' curve class beta (beta . D_S = 1 for every boundary S), compute
        K.beta = sum_S ( s(n-s)/(n-1) - 2 ) = 2^(n-3)(n-8) + n + 2,
     (K+B).beta, #boundary divisors, expected dimension of rational curves of class beta,
     and the integer heuristic exponent 1 - K.beta (per shift class, in units of w),
     and check the closed forms.
"""
from fractions import Fraction as Fr
from itertools import combinations
from math import comb

def kapranov_check(n):
    pts = list(range(1, n))          # q_1..q_{n-1}; the n-th marked point is the psi point
    Is = [frozenset(I) for r in range(1, n - 3) for I in combinations(pts, r)]
    basis = ['H'] + Is
    def vec():
        return {b: Fr(0) for b in basis}
    K = vec(); K['H'] = Fr(-(n - 2))
    for I in Is:
        K[I] = Fr(n - 3 - len(I))
    # boundary divisors: S up to complement; represent each by the side containing n
    total = vec()
    nb = 0
    for r in range(1, n - 2):                  # |I| = r, S = I + {n}, |S| = r+1 in [2, n-2]
        for I in combinations(pts, r):
            I = frozenset(I)
            s = len(I) + 1
            coeff = Fr(s * (n - s), n - 1) - 2
            D = vec()
            if len(I) <= n - 4:
                D[I] = Fr(1)
            else:                               # |I| = n-3: proper transform of a hyperplane
                D['H'] = Fr(1)
                for K2 in Is:
                    if K2 < I:
                        D[K2] -= 1
            for b in basis:
                total[b] += coeff * D[b]
            nb += 1
    assert nb == 2 ** (n - 1) - n - 1
    assert total == K, (n, {b: (total[b], K[b]) for b in basis if total[b] != K[b]})
    return True

for n in range(5, 11):
    kapranov_check(n)
print("(K1) K = sum_S (s(n-s)/(n-1) - 2) D_S verified in Kapranov's basis for n = 5..10")

print("(K2)  n  #bdry   (K+B).beta     K.beta  closed-form  exp.dim(g=0)  heuristic exponent per shift class")
for n in range(4, 13):
    KB = Fr(0); Kb = Fr(0); nb = 0
    for s in range(2, n // 2 + 1):
        cnt = comb(n, s) if 2 * s != n else comb(n, s) // 2
        KB += cnt * (Fr(s * (n - s), n - 1) - 1)
        Kb += cnt * (Fr(s * (n - s), n - 1) - 2)
        nb += cnt
    assert nb == 2 ** (n - 1) - n - 1
    closed = 2 ** (n - 3) * (n - 8) + n + 2 if n >= 3 else None
    assert Kb == closed, (n, Kb, closed)
    assert KB == 2 ** (n - 3) * (n - 4) + 1
    expdim = -Kb + (n - 3) - 3          # rational curves of class beta in Mbar_{0,n}, modulo Aut(P^1)
    m = n - 1
    heur = 2 ** (m - 2) * (7 - m) - m - 2   # angular-only integer count per shift class (check_system_rank.py)
    assert heur == 1 - Kb, (n, heur, Kb)
    print(f"      {n:2d} {nb:5d}   {str(KB):>10}   {str(Kb):>8}   {closed:>8}   {str(expdim):>10}       {heur:>6}")
print("(K2) closed forms (K+B).beta = 2^(n-3)(n-4)+1, K.beta = 2^(n-3)(n-8)+n+2, heuristic = 1 - K.beta: verified")
print("     K.beta > 0  <=>  n >= 8.")
