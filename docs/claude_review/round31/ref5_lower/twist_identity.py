"""Referee's identity for monic-normalised (unit-twisted) points at T=0:
if 0 in S_i, 1 in S_j (distinct forms), A = prod_{v in S_i cap S_j}(theta+v), P_J(X) = sum_k c_k X^{2^k},
O = sum_{k odd} c_k, then F_i/theta - F_j/theta^2 = A * O.  Check, and compare deg O with deg P_J(1)."""
import itertools, random
from f4poly import *

def subspace_coeffs(J):
    """coefficients c_k in F_2[t] (ints) of P_J(X) = prod_{v in J}(X+v) = sum c_k X^{2^k}"""
    # multiply out over F_2[t] as polynomial in X with int (F_2[t]) coefficients
    poly = [1]  # coefficient list in X
    for v in J:
        new = [0] * (len(poly) + 1)
        for i, c in enumerate(poly):
            new[i + 1] ^= c
            new[i] ^= f2mul(c, v)
        poly = new
    n = len(poly) - 1
    out = {}
    for i, c in enumerate(poly):
        if c:
            assert i & (i - 1) == 0, "not additive"
            out[i.bit_length() - 1] = c
    return out

for d in range(3, 9):
    Q = 1 << d
    fs = forms_and_sides(d)
    worst = None
    checked = 0
    pairs = list(itertools.combinations(range(len(fs)), 2))
    for a, b in pairs:
        J = sorted(fs[a][1] & fs[b][1])
        c = subspace_coeffs(J)
        O = 0
        for k, ck in c.items():
            if k % 2 == 1:
                O ^= ck
        PJ1 = 0
        for k, ck in c.items():
            PJ1 ^= ck
        assert f2deg(PJ1) == sum(f2deg(1 ^ v) for v in J)
        loss = f2deg(PJ1) - (f2deg(O) if O else -10**9)
        worst = loss if worst is None else max(worst, loss)
        if d <= 6 and checked < 40:
            # exact check of the identity with S_i = ker_a (contains 0), S_j = 1 + ker_b (contains 1)
            Si, Sj = fs[a][1], fs[b][2]
            Fi = prod(linfac(0, v) for v in sorted(Si))
            Fj = prod(linfac(0, v) for v in sorted(Sj))
            A = prod(linfac(0, v) for v in sorted(Si & Sj))
            lhs = add(mul((INV[2],), Fi), mul((INV[3],), Fj))
            assert lhs == mul(A, from_f2(O)), (d, a, b)
            checked += 1
    print(f"d={d}: identity checked on {checked} pairs; max over pairs of deg P_J(1) - deg O = {worst}")
