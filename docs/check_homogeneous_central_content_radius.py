"""Exact checks for homogeneous central content and least-radius formulas."""

from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import gcd, isqrt, lcm


Gaussian = tuple[int, int]


def mat_vec(q, v):
    return q[0] * v[0] + q[1] * v[1], q[1] * v[0] + q[2] * v[1]


def qeval(q, v):
    w = mat_vec(q, v)
    return v[0] * w[0] + v[1] * w[1]


def polar(q, u, v):
    w = mat_vec(q, v)
    return u[0] * w[0] + u[1] * w[1]


def det(u, v):
    return u[0] * v[1] - u[1] * v[0]


def primitive_vector(values):
    den = lcm(*(x.denominator for x in values))
    raw = [int(den * x) for x in values]
    content = reduce(gcd, (abs(x) for x in raw))
    return tuple(x // content for x in raw), Fraction(den, content)


def prime_factors(n):
    n = abs(n)
    out = set()
    p = 2
    while p * p <= n:
        if n % p == 0:
            out.add(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        out.add(n)
    return out


def valuation(n, p):
    n = abs(n)
    e = 0
    while n % p == 0:
        e += 1
        n //= p
    return e


def mul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def conj(z):
    return z[0], -z[1]


def norm(z):
    return z[0] * z[0] + z[1] * z[1]


def nearest(n, d):
    return (2 * n + d) // (2 * d) if n >= 0 else -((-2 * n + d) // (2 * d))


def divmod_gaussian(z, w):
    n = mul(z, conj(w))
    d = norm(w)
    quotient = nearest(n[0], d), nearest(n[1], d)
    product = mul(quotient, w)
    return quotient, (z[0] - product[0], z[1] - product[1])


def gcd_gaussian(z, w):
    while w != (0, 0):
        _, remainder = divmod_gaussian(z, w)
        z, w = w, remainder
    return z


def exact_div(z, w):
    quotient, remainder = divmod_gaussian(z, w)
    assert remainder == (0, 0)
    return quotient


def lcm_gaussian(z, w):
    return mul(exact_div(z, gcd_gaussian(z, w)), w)


def check_case(q, vectors):
    a, b, c = q
    r = isqrt(a * c - b * b)
    assert r * r == a * c - b * b and gcd(gcd(a, b), c) == 1
    assert all(gcd(abs(x), abs(y)) == 1 for x, y in vectors)
    assert all(det(u, v) != 0 for u, v in combinations(vectors, 2))

    qs = [qeval(q, v) for v in vectors]
    ds = [reduce(lambda x, j: x * det(vectors[i], vectors[j]),
                 (j for j in range(6) if j != i), 1)
          for i in range(6)]
    weights = [Fraction(qs[i] ** 2, ds[i]) for i in range(6)]
    primitive, kappa = primitive_vector(weights)

    primes = set()
    for value in qs + ds + [r]:
        primes |= prime_factors(value)
    predicted_kappa = Fraction(1)
    predicted_valuations = [[0] * 6 for _ in primes]
    for row, p in enumerate(sorted(primes)):
        alpha = [2 * valuation(qs[i], p) - valuation(ds[i], p) for i in range(6)]
        minimum = min(alpha)
        predicted_kappa *= Fraction(p) ** (-minimum)
        predicted_valuations[row] = [x - minimum for x in alpha]
        assert [valuation(x, p) for x in primitive] == predicted_valuations[row]

        signed_degrees = []
        for i in range(6):
            signed_degrees.append(sum(
                valuation(qs[i], p) + valuation(qs[j], p)
                - 2 * valuation(r, p)
                - 2 * valuation(det(vectors[i], vectors[j]), p)
                for j in range(6) if j != i
            ))
        degree_minimum = min(signed_degrees)
        assert [degree - degree_minimum for degree in signed_degrees] == [
            2 * value for value in predicted_valuations[row]
        ]
    assert predicted_kappa == kappa

    edge_ns = []
    chord_sq = {}
    epsilon_seen = set()
    for i, j in combinations(range(6), 2):
        delta = det(vectors[i], vectors[j])
        bij = polar(q, vectors[i], vectors[j])
        assert bij * bij + r * r * delta * delta == qs[i] * qs[j]
        gij = gcd(abs(bij), r * abs(delta))
        aa, tt = bij // gij, r * delta // gij
        epsilon = 2 if aa % 2 and tt % 2 else 1
        epsilon_seen.add(epsilon)
        nij = (aa * aa + tt * tt) // epsilon
        assert gij * gij == gcd(qs[i] * qs[j], r * r * delta * delta)
        assert nij == qs[i] * qs[j] // (gij * gij * epsilon)
        quotient = qs[i] * qs[j] // gcd(qs[i] * qs[j], r * r * delta * delta)
        while quotient % 2 == 0:
            quotient //= 2
        assert nij == quotient and nij % 2 == 1
        assert Fraction(4 * r * r * delta * delta, qs[i] * qs[j]).denominator == nij
        edge_ns.append(nij)
        chord_sq[i, j] = Fraction(4 * r * r * delta * delta, qs[i] * qs[j])

    # Squared form of the five-chord star identity (10).
    common_sq = Fraction((2 * r) ** 10, reduce(lambda x, y: x * y, qs, 1))
    for i in range(6):
        product_sq = Fraction(1)
        for j in range(6):
            if i != j:
                product_sq *= chord_sq[min(i, j), max(i, j)]
        assert weights[i] * weights[i] * product_sq == common_sq

    # Independent Gaussian reconstruction relative to the first phase.
    hs = [(a * x + b * y, r * y) for x, y in vectors]
    denominators = [(1, 0)]
    for h in hs[1:]:
        relative = mul(h, conj(hs[0]))
        common = gcd_gaussian(relative, conj(relative))
        reduced = exact_div(relative, common)
        denominators.append(conj(reduced))
    gaussian_lcm = reduce(lcm_gaussian, denominators)
    all_edge_n = lcm(*edge_ns)
    assert norm(gaussian_lcm) == all_edge_n
    for p in sorted(primes):
        if p == 2:
            assert valuation(all_edge_n, p) == 0
        else:
            radius_order = max(
                0,
                *(
                    valuation(qs[i], p) + valuation(qs[j], p)
                    - 2 * valuation(r, p)
                    - 2 * valuation(det(vectors[i], vectors[j]), p)
                    for i, j in combinations(range(6), 2)
                ),
            )
            assert valuation(all_edge_n, p) == radius_order

    # A nontrivial simultaneous SL_2(Z) change leaves every raw weight fixed.
    # S=[[2,1],[1,1]], S^(-1)=[[1,-1],[-1,2]].
    transformed_q = (4 * a + 4 * b + c, 2 * a + 3 * b + c, a + 2 * b + c)
    transformed_v = [(x - y, -x + 2 * y) for x, y in vectors]
    transformed_weights = []
    for i, v in enumerate(transformed_v):
        di = reduce(lambda x, j: x * det(v, transformed_v[j]),
                    (j for j in range(6) if j != i), 1)
        transformed_weights.append(Fraction(qeval(transformed_q, v) ** 2, di))
    assert transformed_weights == weights

    # Independent projective rescaling changes all raw weights by one scalar.
    scales = (2, -3, 5, 7, -11, 13)
    scaled = [(scales[i] * x, scales[i] * y) for i, (x, y) in enumerate(vectors)]
    scaled_weights = []
    for i, v in enumerate(scaled):
        di = reduce(lambda x, j: x * det(v, scaled[j]),
                    (j for j in range(6) if j != i), 1)
        scaled_weights.append(Fraction(qeval(q, v) ** 2, di))
    scalar = Fraction(1, reduce(lambda x, y: x * y, scales, 1))
    assert scaled_weights == [scalar * w for w in weights]
    assert primitive_vector(scaled_weights)[0] in (primitive, tuple(-x for x in primitive))

    return len(primes), len(edge_ns), all_edge_n, epsilon_seen


def check_evaluation_clearer_is_not_weight_clearer():
    """The p=13 balanced cut has ell_E larger than the minimal ell_W."""
    p = 13
    e = 1
    rho = 5
    modulus = p**e
    bad = next(
        shift for shift in range(p)
        if ((rho + shift * modulus) ** 2 + 1) % (p * modulus) == 0
    )
    start = next(
        shift for shift in range(p)
        if all((shift + offset) % p != bad for offset in range(3))
    )
    clustered = [rho + (start + offset) * modulus for offset in range(3)]
    nodes = [0, 1, *(value - 1 for value in clustered)]

    evaluation_rows = []
    weights = []
    for i, node in enumerate(nodes):
        denominator = reduce(
            lambda value, j: value * (node - nodes[j]),
            (j for j in range(5) if j != i),
            1,
        )
        evaluation_rows.append(
            [Fraction(node**degree, denominator) for degree in range(5)]
        )
        # These signs match homogeneous v_i=(node,1), with v_infinity=(1,0).
        weights.append(Fraction(-(node * node + 2 * node + 2) ** 2, denominator))
    weights.append(Fraction(1))

    ell_e = lcm(*(entry.denominator for row in evaluation_rows for entry in row))
    ell_w = lcm(*(weight.denominator for weight in weights))
    content_e = reduce(gcd, (abs(int(ell_e * weight)) for weight in weights))
    content_w = reduce(gcd, (abs(int(ell_w * weight)) for weight in weights))
    raw_minimum = min(
        valuation(weight.numerator, p) - valuation(weight.denominator, p)
        for weight in weights
    )

    assert raw_minimum == 0
    assert valuation(ell_w, p) == valuation(content_w, p) == 0
    assert valuation(ell_e, p) == valuation(content_e, p) == 2 * e
    assert Fraction(ell_w, content_w) == Fraction(ell_e, content_e)


def main():
    cases = [
        ((1, 0, 1), [(0, 1), (1, 1), (2, 1), (3, 2), (5, 3), (7, 4)]),
        ((1, 1, 2), [(1, 0), (0, 1), (1, 1), (2, 1), (3, 2), (5, 2)]),
        ((2, 1, 5), [(0, 1), (1, 0), (1, 1), (2, 1), (3, 1), (4, 3)]),
    ]
    totals = [check_case(q, vectors) for q, vectors in cases]
    check_evaluation_clearer_is_not_weight_clearer()
    assert set().union(*(item[3] for item in totals)) == {1, 2}
    print("PASS: exact homogeneous central-content and radius checks")
    print(f"cases={len(cases)}, prime rows={sum(x[0] for x in totals)}, "
          f"edge checks={sum(x[1] for x in totals)}")
    print("least squared radii:", *(x[2] for x in totals))
    print("PASS: polar-free edge gcd identity, including both parity cases")
    print("PASS: p=13 balanced cut distinguishes ell_E from ell_W")


if __name__ == "__main__":
    main()
