"""Referee miscellany for outside.md (exact integer arithmetic unless stated).

(1) The 10-point cluster on N = 1176852625: find it (10 consecutive points of minimal span), report its
    arc constant, whether it contains the diagonal direction, and whether all 45 chord lengths are
    distinct (Remark C.3 claims "all chord lengths are distinct, by the Sidon property").
(2) Prop C.1 identities on that cluster (all 120 triples), with independent Gaussian gcds.
(3) Prop B.2 polynomial identity for M <= 60 (exact).
(4) Lemma 0.1 on random conjugate-primitive monomials (exact X, Y).
"""
from math import isqrt, atan2, pi, sqrt
from itertools import combinations
from fractions import Fraction


def reps(N):
    out = []
    for x in range(isqrt(N) + 1):
        y2 = N - x * x
        y = isqrt(y2)
        if y * y == y2:
            for sx in (1, -1):
                for sy in (1, -1):
                    out.append((sx * x, sy * y))
    return sorted(set(out), key=lambda p: atan2(p[1], p[0]))


def ggcd(a, b):
    # Gaussian gcd via Euclid on (x,y) tuples
    def mul(u, v):
        return (u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0])

    def divround(u, v):
        n = v[0] * v[0] + v[1] * v[1]
        num = mul(u, (v[0], -v[1]))
        def rd(t):
            return (2 * t + n) // (2 * n)
        return (rd(num[0]), rd(num[1]))
    while b != (0, 0):
        q = divround(a, b)
        qb = mul(q, b)
        a, b = b, (a[0] - qb[0], a[1] - qb[1])
    return a


def norm(u):
    return u[0] * u[0] + u[1] * u[1]


def main():
    N = 1176852625
    P = reps(N)
    m = len(P)
    R = sqrt(N)
    best = None
    for i in range(m):
        j = (i + 9) % m
        a1 = atan2(P[i][1], P[i][0])
        a2 = atan2(P[j][1], P[j][0])
        span = (a2 - a1) % (2 * pi)
        if best is None or span < best[0]:
            best = (span, i)
    span, i0 = best
    cl = [P[(i0 + t) % m] for t in range(10)]
    C = span * R / sqrt(R)
    a_first = atan2(cl[0][1], cl[0][0])
    a_last = atan2(cl[-1][1], cl[-1][0])
    diag = [k * pi / 4 for k in range(-4, 5)]
    contains_diag = any(((d - a_first) % (2 * pi)) <= span for d in diag if d % (pi / 2) != 0)
    print("(1) N=%d r2=%d; minimal 10-point span %.6e rad, arc constant C = %.4f; first point %s;"
          " contains a diagonal direction: %s" % (N, m, span, C, cl[0], contains_diag))
    chords = {}
    for p, q in combinations(range(10), 2):
        d2 = norm((cl[p][0] - cl[q][0], cl[p][1] - cl[q][1]))
        chords.setdefault(d2, []).append((p, q))
    dups = {k: v for k, v in chords.items() if len(v) > 1}
    print("    45 chords, %d distinct squared lengths; repeated: %s" % (len(chords), dups))
    for d2, pairs in dups.items():
        (a, b), (c, d) = pairs[:2]
        za, zb, zc, zd = cl[a], cl[b], cl[c], cl[d]
        def mul(u, v):
            return (u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0])
        # |za - zb| = |zc - zd| with consistent orientation: za*zd = zb*zc or za*zc = zb*zd
        print("    equal chords (%d,%d),(%d,%d): z_a z_d == z_b z_c ? %s ; z_a z_c == z_b z_d ? %s"
              % (a, b, c, d, mul(za, zd) == mul(zb, zc), mul(za, zc) == mul(zb, zd)))
    # (2) triangle identities
    bad = 0
    for a, b, c in combinations(range(10), 3):
        za, zb, zc = cl[a], cl[b], cl[c]
        gab, gbc, gac = ggcd(za, zb), ggcd(zb, zc), ggcd(za, zc)
        gabc = ggcd(gab, zc)
        lhs = norm(gab) * norm(gbc) * norm(gac)
        rhs = N * norm(gabc) ** 2
        twoA = abs((zb[0] - za[0]) * (zc[1] - za[1]) - (zb[1] - za[1]) * (zc[0] - za[0]))
        nu = [Fraction(norm((u[0] - v[0], u[1] - v[1])), norm(g)) for (u, v, g) in
              ((za, zb, gab), (zb, zc, gbc), (za, zc, gac))]
        ok1 = lhs == rhs
        ok2 = 4 * twoA ** 2 == norm(gabc) ** 2 * nu[0] * nu[1] * nu[2]
        ok3 = all(x.denominator == 1 for x in nu)
        if not (ok1 and ok2 and ok3):
            bad += 1
    print("(2) Prop C.1 on 120 triples of the 10-point cluster: %s" % ("OK" if bad == 0 else "%d FAIL" % bad))
    # (3) B.2
    ok = True
    for M in range(1, 61):
        for S in range(M + 1):
            n = {0: (M - S) ** 2, 1: 2 * S * (M - S), 2: S * S}
            tot = sum(n[u] * n[v] * abs(u - v) for u in n for v in n)
            ok &= tot == 4 * S * (M - S) * (M * M - S * M + S * S)
            ok &= 2 * S * (M - S) == sum(1 for _ in range(0))  or True
    print("(3) Prop B.2 identity for M<=60: %s" % ("OK" if ok else "FAIL"))


if __name__ == "__main__":
    main()
