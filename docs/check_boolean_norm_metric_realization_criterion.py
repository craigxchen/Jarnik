"""Exact arithmetic checks for boolean_norm_metric_realization_criterion.md."""

from fractions import Fraction
from itertools import combinations
from math import gcd, isqrt, prod


def add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def mul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def bar(a):
    return (a[0], -a[1])


def norm(a):
    return a[0] ** 2 + a[1] ** 2


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1]


def det(a, b):
    return a[0] * b[1] - a[1] * b[0]


def exact_div(a, b):
    z = mul(a, bar(b))
    d = norm(b)
    if z[0] % d or z[1] % d:
        return None
    return (z[0] // d, z[1] // d)


def norm_candidates(n):
    for x in range(-isqrt(n), isqrt(n) + 1):
        y2 = n - x * x
        y = isqrt(y2)
        if y * y == y2:
            yield (x, y)
            if y:
                yield (x, -y)


def metric_data(n, delta):
    """Return exact rational Q and the C_k, or None if compatibility fails."""
    m = len(n)
    d = delta[1, 2]
    c = {
        k: n[1] * delta[2, k] ** 2 + n[2] * delta[1, k] ** 2
        - n[k] * d * d
        for k in range(3, m + 1)
    }
    for k in range(4, m + 1):
        if c[k] * delta[1, 3] * delta[2, 3] != c[3] * delta[1, k] * delta[2, k]:
            return None
    if c[3] ** 2 != 4 * delta[1, 3] ** 2 * delta[2, 3] ** 2 * (n[1] * n[2] - d * d):
        return None
    s = Fraction(c[3], 2 * delta[1, 3] * delta[2, 3])
    return ((Fraction(n[1]), s / d, Fraction(n[2], d * d)), c)


def integer_lift(n, delta):
    """Find the note's common divisor G by exact finite norm enumeration."""
    metric = metric_data(n, delta)
    assert metric is not None
    q11, q12, q22 = metric[0]
    d = delta[1, 2]
    v = {1: (Fraction(1), Fraction(0)), 2: (Fraction(0), Fraction(d))}
    for k in range(3, len(n) + 1):
        v[k] = (Fraction(-delta[2, k], d), Fraction(delta[1, k]))
    z = {}
    for j in range(1, len(n) + 1):
        a, b = v[1]
        c, e = v[j]
        sij = q11 * a * c + q12 * (a * e + b * c) + q22 * b * e
        assert sij.denominator == 1
        z[j] = (int(sij), 0 if j == 1 else delta[1, j])
    for g in norm_candidates(n[1]):
        if all(exact_div(w, g) is not None for w in z.values()):
            p = {1: bar(g)}
            for j in range(2, len(n) + 1):
                p[j] = exact_div(z[j], g)
            assert all(norm(p[j]) == n[j] for j in p)
            assert all(det(p[i], p[j]) == delta[i, j]
                       for i, j in combinations(p, 2))
            return p
    return None


def check_abstract_integer_lifts():
    # Rational positive forms of determinant one, including the ramified
    # prime 2 and the inert prime 3, whose first coordinates need rotation.
    for n, ds in [
        ({1: 5, 2: 5, 3: 1}, (5, 1, 2)),
        ({1: 2, 2: 2, 3: 1}, (2, 1, 1)),
        ({1: 9, 2: 9, 3: 2}, (9, 3, 3)),
    ]:
        delta = {(1, 2): ds[0], (1, 3): ds[1], (2, 3): ds[2]}
        assert integer_lift(n, delta) is not None

    # Positive real metric with N_1=3: the Gaussian norm hypothesis matters.
    n = {1: 3, 2: 3, 3: 6}
    delta = {(1, 2): 3, (1, 3): 3, (2, 3): 3}
    assert metric_data(n, delta) is not None
    assert integer_lift(n, delta) is None

    # The integer lift need not be conjugate-primitive.
    p = {1: (5, 0), 2: (3, 4), 3: (2, 1)}
    n = {i: norm(z) for i, z in p.items()}
    delta = {(i, j): det(p[i], p[j]) for i, j in combinations(p, 2)}
    assert integer_lift(n, delta) is not None
    prim25 = [z for z in norm_candidates(25)
              if gcd(abs(z[0]), abs(z[1])) == 1]
    assert all(det(a, b) != 20 for a in prim25 for b in prim25)


def check_boolean_triple():
    # Fifteen disjoint split-prime norms on all nonempty four-row subsets.
    split_primes = [5, 13, 17, 29, 37, 41, 53, 61, 73, 89,
                    97, 101, 109, 113, 137]
    h = {}
    for subset, prime in zip(
        (frozenset(s) for r in range(1, 5) for s in combinations(range(1, 5), r)),
        split_primes,
    ):
        factor = next((a, b) for a in range(1, isqrt(prime) + 1)
                      for b in range(1, isqrt(prime) + 1)
                      if a * a + b * b == prime)
        h[subset] = factor
    core_norm = {s: norm(z) for s, z in h.items()}
    p = {i: (1, 0) for i in range(1, 5)}
    for subset, factor in h.items():
        for i in subset:
            p[i] = mul(p[i], factor)
    n = {i: norm(z) for i, z in p.items()}
    delta = {(i, j): det(p[i], p[j]) for i, j in combinations(p, 2)}
    assert all(x != 0 for x in delta.values())
    assert metric_data(n, delta) is not None

    pair_core = {(i, j): prod(v for s, v in core_norm.items() if i in s and j in s)
                 for i, j in combinations(p, 2)}
    assert all(delta[e] % pair_core[e] == 0 for e in delta)
    t = {e: delta[e] // pair_core[e] for e in delta}
    u = {i: prod(v for s, v in core_norm.items() if s & {1, 2, 3} == {i})
         for i in (1, 2, 3)}
    v = {(i, j): prod(val for s, val in core_norm.items()
                      if s & {1, 2, 3} == {i, j})
         for i, j in combinations((1, 2, 3), 2)}
    common = prod(val for s, val in core_norm.items() if {1, 2, 3} <= s)
    f = v[1, 2] * v[1, 3] * v[2, 3] * common ** 3
    assert pair_core[1, 2] * pair_core[1, 3] * pair_core[2, 3] == f
    c3 = n[1] * delta[2, 3] ** 2 + n[2] * delta[1, 3] ** 2 - n[3] * delta[1, 2] ** 2
    a = t[2, 3] ** 2 * u[1] * v[2, 3]
    b = t[1, 3] ** 2 * u[2] * v[1, 3]
    c = t[1, 2] ** 2 * u[3] * v[1, 2]
    assert c3 == f * (a + b - c)
    s12 = dot(p[1], p[2]) // pair_core[1, 2]
    assert a + b - c == 2 * t[1, 3] * t[2, 3] * s12
    assert s12 ** 2 + t[1, 2] ** 2 == u[1] * u[2] * v[1, 3] * v[2, 3]
    assert (a + b - c) ** 2 + 4 * t[1, 2] ** 2 * t[1, 3] ** 2 * t[2, 3] ** 2 == 4 * a * b

    aa = {i: (1, 0) for i in (1, 2, 3)}
    bb = {e: (1, 0) for e in ((1, 2), (1, 3), (2, 3))}
    for subset, factor in h.items():
        cut = subset & {1, 2, 3}
        if len(cut) == 1:
            i = next(iter(cut))
            aa[i] = mul(aa[i], factor)
        elif len(cut) == 2:
            e = tuple(sorted(cut))
            bb[e] = mul(bb[e], factor)
    w1 = mul((t[2, 3], 0), mul(aa[1], bar(bb[2, 3])))
    w2 = mul((-t[1, 3], 0), mul(aa[2], bar(bb[1, 3])))
    w3 = mul((t[1, 2], 0), mul(aa[3], bar(bb[1, 2])))
    assert add(add(w1, w2), w3) == (0, 0)
    assert (norm(w1), norm(w2), norm(w3)) == (a, b, c)
    assert det(w1, w2) == -t[1, 2] * t[1, 3] * t[2, 3]


if __name__ == "__main__":
    check_abstract_integer_lifts()
    check_boolean_triple()
    print("PASS: rational metric, Gaussian integer lift, primitive obstruction, and Boolean triple normalization")
