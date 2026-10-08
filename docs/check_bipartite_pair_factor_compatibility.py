"""Check the exact bipartite pair-factor profile and its four-cycle budget."""

from decimal import Decimal, getcontext
from itertools import combinations


getcontext().prec = 70


def cuts_for(m):
    k = m // 2
    return [frozenset(u + v)
            for u in combinations(range(k), k // 2)
            for v in combinations(range(k, m), k // 2)]


def signatures(m, cuts):
    return {pair: tuple((pair[0] in s) != (pair[1] in s) for s in cuts)
            for pair in combinations(range(m), 2)}


def split_primes(n):
    out = []
    p = 5
    while len(out) < n:
        if p % 4 == 1 and all(p % d for d in range(2, int(p ** .5) + 1)):
            out.append(p)
        p += 2
    return out


def least_odd_weight(p, t):
    lp = Decimal(p).ln()
    s = int(t / lp)
    if s % 2 == 0:
        s += 1
    while Decimal(s) * lp < t:
        s += 2
    return s, Decimal(s) * lp


def check():
    m = 12
    k = m // 2
    for size in (8, 12):
        cs = cuts_for(size)
        distinct = len(set(signatures(size, cs).values()))
        assert distinct == (22 if size == 8 else size * (size - 1) // 2)

    cuts = cuts_for(m)
    l = len(cuts)
    assert l == 400
    ps = split_primes(l + 1)
    t = Decimal(100000)
    exponents, weights = zip(*(least_odd_weight(p, t) for p in ps[:-1]))
    err_bound = 2 * sum((Decimal(p).ln() for p in ps[:-1]), Decimal(0))
    lam = Decimal(16).ln()  # C=1/2: log(4/C^2)=log 16.
    star = ps[-1]
    star_log = Decimal(star).ln()
    sstar = 1
    while Decimal(sstar) * star_log <= err_bound + 2 * lam + 10:
        sstar += 2
    v = Decimal(sstar) * star_log
    assert v > err_bound + 2 * lam
    all_cuts = cuts + [frozenset(range(k))]
    all_exponents = exponents + (sstar,)
    all_weights = weights + (v,)
    assert all(len(cut) == k for cut in all_cuts)
    w = sum(all_weights, Decimal(0))
    upper = w / 2 + w / m
    lower = w / 2 + lam
    sig = signatures(m, all_cuts)
    assert all(sum(bits[j] for bits in sig.values()) == m * m // 4
               for j in range(len(all_cuts)))
    good = []
    for pair, bits in sig.items():
        d = sum((wp for wp, separated in zip(all_weights, bits) if separated), Decimal(0))
        assert d > lower
        is_cross = (pair[0] < k) != (pair[1] < k)
        assert (d < upper) == is_cross
        if d < upper:
            good.append(pair)
    assert len(good) == k * k

    # Signed local exponent: positive selects pi, negative selects bar(pi).
    def delta(i, j, cut, s):
        return s * (int(i in cut) - int(j in cut))

    triangles = 0
    for a, b, c in combinations(range(m), 3):
        for cut, s in zip(all_cuts, all_exponents):
            x, y, z = delta(a, b, cut, s), delta(b, c, cut, s), delta(a, c, cut, s)
            q = min(abs(x), abs(y)) if x * y < 0 else 0
            assert x + y == z
            assert abs(x) + abs(y) == abs(z) + 2 * q
        triangles += 1

    cycles = 0
    beta_count = 72  # alpha^2/2 * 400 for M=12.
    for a, c in combinations(range(k), 2):
        for b, d in combinations(range(k, m), 2):
            support1 = support2 = 0
            odd_reduced = False
            for cut, s in zip(all_cuts, all_exponents):
                u1, u2 = delta(a, b, cut, s), delta(c, d, cut, s)
                v1, v2 = delta(a, d, cut, s), delta(c, b, cut, s)
                pos1 = max(u1, 0) + max(u2, 0)
                neg1 = max(-u1, 0) + max(-u2, 0)
                pos2 = max(v1, 0) + max(v2, 0)
                neg2 = max(-v1, 0) + max(-v2, 0)
                y1, y2 = pos1 - neg1, pos2 - neg2
                q1, q2 = min(pos1, neg1), min(pos2, neg2)
                assert y1 == y2
                assert q1 * q2 == 0
                assert q1 in (0, s) and q2 in (0, s)
                support1 += q1 > 0
                support2 += q2 > 0
                odd_reduced |= bool(y1 % 2)
            assert support1 == support2 == beta_count
            assert odd_reduced
            cycles += 1
    assert cycles == 225
    print(f"PASS: M={m}, {l} balanced source cuts, {len(good)} bipartite good pairs, "
          f"{triangles} triangles and {cycles} exact four-cycle reductions; "
          "distinct nonsquare pair norms and a finite C=1/2 norm-window instance.")


if __name__ == "__main__":
    check()
