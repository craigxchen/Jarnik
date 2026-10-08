"""Exact eleven-prefix reciprocal-Gram obstruction and tail budget."""

from fractions import Fraction as Q
from itertools import product


MASKS = (128, 64, 206, 205, 35, 19, 240, 8, 4, 60, 195)
DELTA = (0, 0, 0, 0, -1, -1, 1, 1)


if __name__ == "__main__":
    S = [[1 if mask >> i & 1 else -1 for mask in MASKS] for i in range(8)]
    dA = [(x - y) // 2 for x, y in zip(S[0], S[1])]
    dB = [(x - y) // 2 for x, y in zip(S[4], S[5])]
    assert [(j + 1, x) for j, x in enumerate(dA) if x] == [(3, -1), (4, 1)]
    assert [(j + 1, x) for j, x in enumerate(dB) if x] == [(5, -1), (6, 1)]
    assert all(x * y == 0 for x, y in zip(dA, dB))
    heights, total = [], [0] * 8
    for t in range(5):
        r = [(total[i] - total[0]) // 2 for i in range(8)]
        heights.append(sum(abs(r[i] - DELTA[i]) for i in range(1, 8)))
        total = [v + S[i][t] for i, v in enumerate(total)]
    assert heights == [4, 3, 2, 7, 4]
    assert Q(1, 6) - Q(1, 8) == Q(1, 24)
    assert Q(1, 10) - Q(1, 12) == Q(1, 60)
    critical_prefix = [(j + 1, S[1][j]) for j in range(11) if S[1][j] != S[6][j]]
    assert critical_prefix == [(2, -1), (4, -1), (5, 1), (6, 1), (7, -1)]
    assert -Q(1, 4) - Q(1, 8) + Q(1, 10) == -Q(11, 40)
    assert Q(1, 4)**2 + Q(1, 8)**2 + Q(1, 10)**2 == Q(141, 1600)
    # The sign tail is a repetition of (-,+,+,-); each block sums to zero.
    partial, tail = 0, (-1, 1, 1, -1)
    for value in tail:
        partial += value
        assert abs(partial) <= 1
    assert partial == 0
    assert sum(tail[:2]) == 0  # handles either parity of the critical degree
    assert 141 - 11**2 == 20  # rationalization of 40/(sqrt(141)-11)
    assert 141 < 12**2        # 2(sqrt(141)+11)<46

    # A degree-two polynomial identity; three values per variable suffice.
    for b0, b1, b4, b5 in product(range(3), repeat=4):
        def cross(x, y):
            return -Q((x - y) ** 2, 2)
        rectangle = cross(b0, b4) - cross(b0, b5) - cross(b1, b4) + cross(b1, b5)
        assert rectangle == (b0 - b1) * (b4 - b5)

    # Exact coefficient checks in N: both sides are polynomials of degree <=2.
    for N in (0, 1, 2):
        H = 60 * N + 40
        assert (Q(H, 24) - N) * (Q(H, 60) - N) == N + Q(10, 9)
    # 128 n^2-[60(n-11)+40] at n=12+t has all positive coefficients.
    coefficients = (128 * 12**2 - 60 * 12 + 620, 256 * 12 - 60, 128)
    assert all(x > 0 for x in coefficients)

    def polynomial(B):
        return B**4 - 26 * B**3 - 124 * B**2 + 6120 * B + 28000

    # After clearing denominators, these are polynomial identities of
    # degree at most six, so seven distinct nodes certify the identities.
    for B in range(14, 21):
        A = B - 2
        upper = Q(11, 40) - Q(1, A) + Q(1, B) + Q(1, B + 10)
        gap = Q(141, 1600) + Q(1, A**2) + Q(1, B**2) - upper**2
        assert gap == Q(polynomial(B), 80 * B * (B - 2) * (B + 10)**2)
    # Each translated numerator is degree four: five nodes suffice.
    for t in range(5):
        assert polynomial(14 + t) == 56448 - 1664*t - 40*t**2 + 30*t**3 + t**4
        assert polynomial(20 + t) == 52800 + 1960*t + 716*t**2 + 54*t**3 + t**4
    assert 56448 - 1664*6 - 40*36 == 45024 > 0
    assert all(c > 0 for c in (52800, 1960, 716, 54, 1))
    # Cancellation of u^2 gives the claimed affine, increasing dependence.
    for u, v, w in product((Q(0), Q(1), Q(2)), repeat=3):
        P, c = Q(11, 40), Q(141, 1600)
        assert c + u*u + v*v - (P-u+v+w)**2 == c + v*v - (P+v+w)**2 + 2*u*(P+v+w)

    # Scope fixture: all sixteen prefix deficits reverse at a_j=2^j.
    reciprocal = [Q(1, 2**j) for j in range(1, 12)]
    gaps, ratios = [], []
    for i in range(4):
        for k in range(4, 8):
            support = [j for j in range(11) if S[i][j] != S[k][j]]
            p1 = sum(S[i][j] * reciprocal[j] for j in support)
            p2 = sum(reciprocal[j]**2 for j in support)
            gaps.append((p1*p1-p2, (i, k)))
            ratios.append((p1*p1/p2, (i, k)))
    assert min(gaps) == (Q(695, 1048576), (0, 5))
    assert min(ratios) == (Q(1062961, 1048901), (3, 7))
    beta = [sum(S[i][j]*reciprocal[j] for j in range(11)) for i in range(8)]
    assert sorted(range(8), key=beta.__getitem__) == [4, 5, 0, 1, 2, 3, 6, 7]
    for v in range(-4, 5):
        assert v - max(v, 0) + max(-v, 0) == 0
    print("Passed prefix supports, axial heights, rectangle, critical-tail bounds, and every-shift spacing exclusion.")
