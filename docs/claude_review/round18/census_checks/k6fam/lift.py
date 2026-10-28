"""Find interior rational points of the source curve E (a^2 = f_b(b) = f_c(g)) by taking points of
E' of the form Q + 4R and solving the three quadratics; exact."""
from fractions import Fraction as Fr
from math import isqrt
from fam import family, disc_T
from ell import quotient_curve, is_sq, ev

def roots_rational(PQ, T):
    P, Q = PQ
    A = [P[k] - T * Q[k] for k in range(3)]
    if A[2] == 0:
        return None
    D = A[1] ** 2 - 4 * A[2] * A[0]
    r = is_sq(D)
    if r is None:
        return None
    return [(-A[1] + r) / (2 * A[2]), (-A[1] - r) / (2 * A[2])]

def lift_T(F, T):
    a = is_sq(T)
    if a is None:
        return None
    rb = roots_rational(F['maps']['b'], T)
    rc = roots_rational(F['maps']['c'], T)
    if rb is None or rc is None:
        return None
    return [(s * a, b, c) for s in (1, -1) for b in rb for c in rc]

if __name__ == '__main__':
    import sys
    xs = [Fr(2), Fr(3), Fr(-5), Fr(1)]
    F, C, g, pts = quotient_curve(xs)
    k = g[3]
    found = []
    Q = pts
    # candidate combos
    import itertools
    cands = []
    for n1, n2, n3, n4 in itertools.product(range(-4, 5), repeat=4):
        if (n1, n2, n3, n4) == (0, 0, 0, 0):
            continue
        cands.append((n1, n2, n3, n4))
    cands.sort(key=lambda v: sum(abs(x) for x in v))
    seen = set()
    for v in cands[:400]:
        R = None
        for n, P in zip(v, Q):
            if n < 0:
                P = (P[0], -P[1]); n = -n
            for _ in range(n):
                R = C.add(R, P)
        if R is None:
            continue
        T = R[0] / k
        if T in seen:
            continue
        seen.add(T)
        L = lift_T(F, T)
        if L:
            found.append((v, T, L))
    print('lifted', len(found), 'of', len(seen), 'distinct T tried')
    for v, T, L in found[:8]:
        print(v, 'height(T) digits', len(str(T.numerator)), 'point', [str(c) for c in L[0]][:3] if len(str(L[0][0])) < 80 else '(large)')
