"""Exhaustive check of the nonunit square-root reduction of Part VIII.2, q <= QMAX.

For each q and each a mod q with at least one root of x^2 = a (mod q):
  v_p = min(e, v_p(a)); if v = e remove p^ceil(e/2); else v = 2t, every root has
  valuation exactly t, remove p^t, keep p^(e-2t).  u = product removed, Q1 = kept.
Check: u | every root; q <= u^2 Q1; x -> x/u mod Q2 well defined (period);
at every p^f || Q2 (Q2 = Q1 with one factor 2 removed) the values x/u occupy
at most two classes; record how many classes occur at 2^f || Q1 (before the
halving), whether the two classes differ by a unit, and the largest number of
roots in a real interval of length sqrt(q).
"""
import sys
from math import isqrt, gcd
from sympy import factorint

QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 400


def vp(n, p):
    if n == 0:
        return 10 ** 9
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def max_in_window(roots_sorted, length_sq, q):
    # roots over two periods; window [s, s+L], L = sqrt(q); count x with x - s <= L,
    # i.e. (x_j - x_i)^2 <= q
    xs = roots_sorted + [r + q for r in roots_sorted]
    best = 0
    j = 0
    for i in range(len(roots_sorted)):
        j = max(j, i)
        while j + 1 < len(xs) and (xs[j + 1] - xs[i]) ** 2 <= length_sq:
            j += 1
        best = max(best, j - i + 1)
    return best


def main():
    nclasses = 0
    four_at_two = 0
    nonunit_two_diff = 0
    worst = (0, None)
    hist = {}
    for q in range(1, QMAX + 1):
        fac = factorint(q)
        roots = {}
        for x in range(q):
            roots.setdefault(x * x % q, []).append(x)
        for a, rs in roots.items():
            nclasses += 1
            u, Q1 = 1, 1
            for p, e in fac.items():
                v = min(e, vp(a, p))
                if v == e:
                    c = -(-e // 2)
                    assert all(vp(x, p) >= c for x in rs)
                    u *= p ** c
                else:
                    assert v % 2 == 0, (q, a, p)
                    t = v // 2
                    assert all(vp(x, p) == t for x in rs), (q, a, p)
                    u *= p ** t
                    Q1 *= p ** (e - 2 * t)
            assert q <= u * u * Q1
            assert all(x % u == 0 for x in rs)
            Q2 = Q1 // 2 if Q1 % 2 == 0 else Q1
            assert (q // u) % Q1 == 0  # x -> x/u mod Q1 is well defined
            for p, f in factorint(Q1).items():
                pf = p ** f
                cls1 = sorted({(x // u) % pf for x in rs})
                if p == 2 and len(cls1) > 2:
                    four_at_two += 1
                    assert f >= 3 and len(cls1) == 4
                if p != 2:
                    assert len(cls1) <= 2
                    if len(cls1) == 2:
                        assert (cls1[1] - cls1[0]) % p != 0
                        assert (cls1[1] + cls1[0]) % pf == 0
            for p, f in factorint(Q2).items():
                pf = p ** f
                cls = sorted({(x // u) % pf for x in rs})
                assert len(cls) <= 2, (q, a, p, f, cls)
                if len(cls) == 2 and (cls[1] - cls[0]) % p == 0:
                    nonunit_two_diff += 1
            # also the 2-adic sign claim: classes mod 2^(f-1) are +-z
            m = max_in_window(sorted(rs), q, q)
            hist[m] = hist.get(m, 0) + 1
            if m > worst[0]:
                worst = (m, (q, a))
    print("q<=", QMAX, "nonempty classes", nclasses)
    print("2-adic layers with 4 classes before halving:", four_at_two)
    print("Q2-layers with two classes differing by a nonunit (p=2 only):", nonunit_two_diff)
    print("max roots in a closed interval of length sqrt(q):", worst)
    print("histogram:", dict(sorted(hist.items())))


if __name__ == "__main__":
    main()
