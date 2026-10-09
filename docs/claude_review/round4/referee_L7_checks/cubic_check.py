"""EXACT: psi_i-degrees of the 'cubic meeting 13 planes' class of L7 Prop 3, the cone decomposition
beta = (35/13) sym(cubic) + (21/26) sym(f), and the linear tests (B+8K), 7psi_i+15K.  Also recomputes
psi_i two ways (two choices of the auxiliary pair j,k) as a consistency check of psi_coeffs."""
import sys, itertools
from fractions import Fraction as Fr
sys.path.insert(0, '../L7')
import classes7 as C
planes = list(itertools.combinations(range(1, 7), 3))
best = None
for choice in itertools.combinations(planes, 13):
    vec = C.kapranov(3, {}, {}, {p: 1 for p in choice})
    if min(vec.values()) >= 0: best = choice; break
cub = C.kapranov(3, {}, {}, {p: 1 for p in best})
print('planes met:', best)
print('unmet planes:', [p for p in planes if p not in best])
def psi(i, vec, alt=False):
    others = [x for x in C.P if x != i]
    j, k = (others[:2] if not alt else others[-2:])
    s = 0
    for S in C.SPL:
        side = S if i in S else frozenset(C.P) - S
        if i in side and j not in side and k not in side: s += vec[S]
    return s
print('psi degrees     :', [psi(i, cub) for i in C.P])
print('psi degrees(alt):', [psi(i, cub, True) for i in C.P])
pair = sum(cub[S] for S in C.SPL if C.is_pair(S)); trip = sum(cub[S] for S in C.SPL if C.is_triple(S))
mK = -C.dot(C.Kc, cub)
print('pair', pair, 'triple', trip, '-K', mK, '(B+8K)', trip - Fr(5, 3)*pair,
      '7psi_i+15K:', [7*psi(i, cub) - 15*mK for i in C.P])
# cone decomposition on S7-invariant parts: averages over pair/triple divisors
a = Fr(35, 13); b = Fr(21, 26)
print('pair-average check :', a*Fr(pair, 21) + b*Fr(6, 21), ' triple-average check:', a*Fr(trip, 35) + b*0)
# Kapranov 3-planes condition: every 4-subset A of [6] contains <= 3 met planes
print('max met planes inside a 4-set:', max(sum(1 for p in best if set(p) <= set(A)) for A in itertools.combinations(range(1, 7), 4)))
