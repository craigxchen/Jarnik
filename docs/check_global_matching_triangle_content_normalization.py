"""Exact global matching / triangle content and center-normalization audit."""

from itertools import combinations
from math import comb, gcd, isqrt, prod
from random import Random

from check_quartet_matching_gcd_cut_budget import ggcd, matching_content, mul, norm, sub


def factors(n):
    answer = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            answer[p] = answer.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        answer[n] = answer.get(n, 0) + 1
    return answer


def vp(n, p):
    assert n
    value = 0
    while n % p == 0:
        n //= p
        value += 1
    return value


def gaussian_prime(p):
    if p == 2:
        return 1, 1
    if p % 4 == 3:
        return p, 0
    for a in range(1, isqrt(p) + 1):
        b = isqrt(p - a * a)
        if a * a + b * b == p:
            return a, b
    raise AssertionError(p)


def vg(z, pi):
    assert z != (0, 0)
    d, value = norm(pi), 0
    while True:
        numerator = mul(z, (pi[0], -pi[1]))
        if numerator[0] % d or numerator[1] % d:
            return value
        z = numerator[0] // d, numerator[1] // d
        value += 1


def list_gcd(rows):
    answer = (0, 0)
    for row in rows:
        answer = ggcd(answer, row)
    return answer


def audit_eight_source_contents(points, source=None):
    """Check the exact source-gcd target dictionary, not its open bound."""
    assert len(points) == len(set(points)) == 8
    N = norm(points[0])
    assert all(norm(z) == N for z in points)
    assert norm(list_gcd(points)) == 1
    source = factors(N) if source is None else source
    assert prod(p ** e for p, e in source.items()) == N
    H = [list_gcd(z for k, z in enumerate(points) if k != i) for i in range(8)]
    A, F, G2 = prod(norm(h) for h in H), 1, 1
    pair_blocks = {}
    for i, j in combinations(range(8), 2):
        deletion = list_gcd(z for k, z in enumerate(points) if k not in (i, j))
        denominator = mul(H[i], H[j])
        numerator = mul(deletion, (denominator[0], -denominator[1]))
        nd = norm(denominator)
        assert numerator[0] % nd == numerator[1] % nd == 0
        L = numerator[0] // nd, numerator[1] // nd
        pair_blocks[i, j] = L
        F *= norm(L)
        G2 *= norm(deletion)
    layer_products = [1] * 5
    for p, e in source.items():
        assert p % 4 == 1
        levels = [vg(z, gaussian_prime(p)) for z in points]
        assert min(levels) == 0 and max(levels) == e
        for h in range(1, e + 1):
            r = sum(level >= h for level in levels)
            layer_products[min(r, 8 - r)] *= p
    p1, p2, p3, p4 = layer_products[1:]
    assert N == p1 * p2 * p3 * p4
    assert A == p1 and F == p2 and G2 == A ** 7 * F
    assert all(norm(ggcd(a, b)) == 1 for a, b in combinations(H, 2))
    assert all(norm(ggcd(a, b)) == 1
               for a, b in combinations(pair_blocks.values(), 2))
    assert all(norm(ggcd(H[k], block)) == 1
               for (i, j), block in pair_blocks.items()
               for k in range(8) if k not in (i, j))
    assert all(gcd(norm(a), norm(b)) == 1
               for (edge_a, a), (edge_b, b) in combinations(pair_blocks.items(), 2)
               if set(edge_a) & set(edge_b))
    for p in source:
        edges = [edge for edge, block in pair_blocks.items() if norm(block) % p == 0]
        assert len(edges) <= 2
        assert all(not (set(a) & set(b)) for a, b in combinations(edges, 2))
    # Cross multiplication of both exact expressions for exp(D+W4-2W).
    assert A ** 8 * F ** 3 * p3 * p4 == N * p1 ** 7 * p2 ** 2
    assert G2 ** 3 == A ** 13 * A ** 8 * F ** 3


def audit(points):
    m, N = len(points), norm(points[0])
    assert all(norm(z) == N for z in points)
    assert norm(list_gcd(points)) == 1
    triangles = []
    for i, j, k in combinations(range(m), 3):
        a, b = sub(points[j], points[i]), sub(points[k], points[i])
        triangles.append(abs(a[0] * b[1] - a[1] * b[0]))
    assert all(triangles)
    g = gcd(*triangles)
    D = list_gcd(sub(z, points[0]) for z in points[1:])
    G = list_gcd(matching_content(q) for q in combinations(points, 4))
    H = [list_gcd(z for j, z in enumerate(points) if i != j) for i in range(m)]
    E = prod(norm(h) for h in H)
    ND, NG = norm(D), norm(G)
    assert N % E == 0 and g % ND == 0
    assert (4 * E * g * g) % (NG * ND) == 0
    assert (2 * E * g * g) % NG == 0
    assert gcd(N, ND) == 1 and ND % 2 == 0
    assert all(norm(ggcd(a, b)) == 1 for a, b in combinations(H, 2))
    K = comb(m, 3)
    b = comb(m // 2, 3) + comb((m + 1) // 2, 3)
    d = comb(m - 1, 3) - b
    assert prod(t // g for t in triangles) % (N ** b * E ** d) == 0

    source = factors(N)
    if m == 8:
        audit_eight_source_contents(points, source)
    primes = set(source) | set(factors(g)) | set(factors(ND)) | set(factors(NG)) | {2}
    for p in primes:
        pi = gaussian_prime(p)
        pibar = pi[0], -pi[1]
        deltas = {(i, j): vg(sub(points[j], points[i]), pi)
                  for i, j in combinations(range(m), 2)}
        if p in source:
            e = source[p]
            levels = [vg(z, pi) for z in points]
            ordered = sorted(levels)
            assert ordered[0] == 0 and ordered[-1] == e
            a, bb = ordered[1], ordered[-2]
            assert vp(E, p) == a + e - bb
            if len(set(levels)) >= 3:
                assert vp(g, p) == 0
            else:
                residue = min(value - levels[i] for (i, j), value in deltas.items()
                              if levels[i] == levels[j])
                assert vp(g, p) == residue
            if a == e:
                expected = e + residue, residue
            elif bb == 0:
                expected = residue, e + residue
            else:
                expected = a, e - bb
            assert (vg(G, pi), vg(G, pibar)) == expected
            assert vg(D, pi) == vg(D, pibar) == 0
        else:
            r = min(deltas.values())
            tau = min(sum(deltas[pair] for pair in combinations(triple, 2))
                      for triple in combinations(range(m), 3))
            h = vg(G, pi)
            classes = []
            for i in range(m):
                for group in classes:
                    if deltas[tuple(sorted((i, group[0])))] > r:
                        group.append(i)
                        break
                else:
                    classes.append([i])
            if len(classes) >= 3:
                assert (tau, h) == (3 * r, 2 * r)
            else:
                u = min(deltas[tuple(sorted(pair))] - r for group in classes
                        for pair in combinations(group, 2))
                assert tau == 3 * r + u
                assert h == 2 * r + (u if min(map(len, classes)) == 1 else 0)
            if p == 2:
                assert r >= 1 and len(classes) == 2
                assert (vp(g, 2), vp(ND, 2), vp(NG, 2)) == ((tau - 2) // 2, r, h)
                assert tau % 2 == 0
            else:
                assert (vp(g, p), vp(ND, p), vp(NG, p)) == (tau, 2 * r, 2 * h)

    assert norm(ggcd(points[0], D)) == 1
    denominator = ND // gcd(abs(D[0]), abs(D[1]))
    numerator = mul(points[0], (D[0], -D[1]))
    assert denominator == ND // gcd(ND, *numerator)
    assert denominator * denominator >= ND


def main():
    circles = {}
    for x in range(-31, 32):
        for y in range(-31, 32):
            N = x * x + y * y
            if 0 < N <= 1000:
                circles.setdefault(N, []).append((x, y))
    circles = [points for points in circles.values() if len(points) >= 5]
    rng, checked = Random(20260915), 0
    while checked < 1200:
        circle = rng.choice(circles)
        points = rng.sample(circle, rng.randrange(5, min(9, len(circle)) + 1))
        if norm(list_gcd(points)) != 1:
            continue
        audit(points)
        checked += 1
    print(f"PASS: {checked:,} primitive circle tuples; exact local formulas,")
    print("content and strengthened triangle divisibilities, reduced center denominators.")
    # Larger repeated-prime fixtures check that singleton and pair contents
    # may use different layers of the same Gaussian prime.
    primes = [(2, 1), (3, 2), (4, 1), (5, 2), (6, 1)]
    checked = 0
    while checked < 512:
        points, source = [(1, 0)] * 8, {}
        for pi in primes:
            e = rng.randrange(1, 7)
            source[norm(pi)] = e
            levels = [rng.randrange(e + 1) for _ in points]
            a, b = rng.sample(range(8), 2)
            levels[a], levels[b] = 0, e
            for i, level in enumerate(levels):
                for j in range(e):
                    factor = pi if j < level else (pi[0], -pi[1])
                    points[i] = mul(points[i], factor)
        points = [mul(z, rng.choice([(1, 0), (0, 1), (-1, 0), (0, -1)]))
                  for z in points]
        if len(set(points)) != 8:
            continue
        audit_eight_source_contents(points, source)
        checked += 1
    print(f"PASS: {checked} additional primitive eight-point source deletion dictionaries;")
    print("mixed units, prime exponents through six, exact intrinsic-target translation.")


if __name__ == "__main__":
    main()
