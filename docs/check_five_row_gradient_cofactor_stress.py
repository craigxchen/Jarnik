"""Exact checks for the five-row gradient cofactor obstruction.

This constructs stresses on literal source tuples, not gradients of a
single low-height invariant and not unbounded endpoint examples.
Only the Python standard library is required.
"""

from itertools import combinations
from math import gcd, prod


ROWS = range(5)
CUTS = range(1, 32)
TRIPLES = list(combinations(ROWS, 3))
MINIMUM = {1: 0, 2: 0, 3: 3, 4: 9, 5: 15}


def popcount(mask):
    return bin(mask).count("1")


def grad_order(mask, i):
    return max(0, 2 * popcount(mask & ~(1 << i)) - 4)


def minor_order(mask, triple):
    a = sum(bool(mask & (1 << i)) for i in triple)
    return sum(grad_order(mask, i) for i in triple) + a * (a - 1) // 2


def check_orders():
    for mask in CUTS:
        s = popcount(mask)
        assert min(minor_order(mask, t) for t in TRIPLES) == MINIMUM[s]
        for triple in TRIPLES:
            a = sum(bool(mask & (1 << i)) for i in triple)
            expected = int((s == 2 and a == 2) or (s == 3 and a == 1))
            assert minor_order(mask, triple) - MINIMUM[s] == expected
    assert sum(MINIMUM[popcount(mask)] for mask in CUTS) == 90
    for i in ROWS:
        assert sum(grad_order(mask, i) for mask in CUTS) == 24
    for triple in TRIPLES:
        assert sum(minor_order(mask, triple) for mask in CUTS) == 96
        assert sum(minor_order(mask, triple) - MINIMUM[popcount(mask)]
                   for mask in CUTS) == 6


def mul(a, b):
    return (a[0] * b[0] - a[1] * b[1],
            a[0] * b[1] + a[1] * b[0])


def conj(a):
    return (a[0], -a[1])


def norm(a):
    return a[0] ** 2 + a[1] ** 2


def power(a, exponent):
    ans = (1, 0)
    for _ in range(exponent):
        ans = mul(ans, a)
    return ans


def gaussian_gcd(a, b):
    while b != (0, 0):
        numerator = mul(a, conj(b))
        denominator = norm(b)
        # Exact nearest-integer quotients, valid for negative numerators.
        quotient = tuple((2 * n + denominator) // (2 * denominator)
                         for n in numerator)
        product = mul(quotient, b)
        a, b = b, (a[0] - product[0], a[1] - product[1])
    return a


def is_prime(n):
    if n < 2:
        return False
    return all(n % d for d in range(2, int(n ** 0.5) + 1))


def split_primes(count):
    found = {}
    for a in range(1, 60):
        for b in range(1, a):
            n = a * a + b * b
            if n % 4 == 1 and is_prime(n):
                found.setdefault(n, (a, b))
    return [found[n] for n in sorted(found)[:count]]


def bracket(a, b):
    return a[0] * b[1] - a[1] * b[0]


def det3(columns):
    a, b, c = columns
    return (a[0] * (b[1] * c[2] - b[2] * c[1])
            - b[0] * (a[1] * c[2] - a[2] * c[1])
            + c[0] * (a[1] * b[2] - a[2] * b[1]))


def matvec(columns, vector):
    return tuple(sum(columns[i][r] * vector[i] for i in ROWS)
                 for r in range(3))


def literal_source(correction_exponents):
    primes = split_primes(32)
    assert len(primes) == 32
    blocks = dict(zip(CUTS, primes[:31]))
    weights = {mask: norm(blocks[mask]) for mask in CUTS}
    corrections = [power(primes[31], e) for e in correction_exponents]
    points = []
    for i in ROWS:
        p = corrections[i]
        for mask in CUTS:
            if mask & (1 << i):
                p = mul(p, blocks[mask])
        assert norm(gaussian_gcd(p, conj(p))) == 1
        assert gcd(*p) == 1
        points.append(p)

    brackets, residues, b_factors = {}, {}, {}
    for i, j in combinations(ROWS, 2):
        delta = bracket(points[i], points[j])
        assert delta != 0
        brackets[i, j] = delta
        pair_gcd_norm = norm(gaussian_gcd(points[i], points[j]))
        assert delta % pair_gcd_norm == 0
        residues[i, j] = delta // pair_gcd_norm
        core = prod(weights[mask] for mask in CUTS
                    if mask & (1 << i) and mask & (1 << j))
        assert pair_gcd_norm % core == 0
        assert delta % core == 0
        b_factors[i, j] = abs(delta) // core
        assert pair_gcd_norm // core <= norm(corrections[i]) * norm(corrections[j])
    residue_bound = max(map(abs, residues.values()))
    correction_norm_bound = max(map(norm, corrections))
    assert max(b_factors.values()) <= correction_norm_bound ** 2 * residue_bound

    divisors = [prod(weights[mask] ** grad_order(mask, i) for mask in CUTS)
                for i in ROWS]
    columns = [(divisors[i] * p[0] ** 2,
                divisors[i] * p[0] * p[1],
                divisors[i] * p[1] ** 2)
               for i, p in enumerate(points)]
    common = prod(weights[mask] ** MINIMUM[popcount(mask)] for mask in CUTS)
    minors = {}
    for triple in TRIPLES:
        determinant = det3([columns[i] for i in triple])
        i, j, k = triple
        assert determinant == (divisors[i] * divisors[j] * divisors[k]
                               * brackets[i, j] * brackets[i, k] * brackets[j, k])
        assert determinant != 0 and determinant % common == 0
        minors[triple] = determinant // common
        residual_core = prod(weights[mask] for mask in CUTS
                             if ((popcount(mask) == 2
                                  and sum(bool(mask & (1 << h)) for h in triple) == 2)
                                 or (popcount(mask) == 3
                                     and sum(bool(mask & (1 << h)) for h in triple) == 1)))
        assert abs(minors[triple]) == (residual_core * b_factors[i, j]
                                      * b_factors[i, k] * b_factors[j, k])

    def circuit(support):
        vector = [0] * 5
        for j, i in enumerate(support):
            triple = tuple(h for h in support if h != i)
            vector[i] = (-1) ** j * minors[triple]
            assert vector[i] != 0
        assert matvec(columns, vector) == (0, 0, 0)
        return vector

    u, v = circuit((0, 1, 2, 3)), circuit((0, 1, 2, 4))
    candidates = [[a + t * b for a, b in zip(u, v)] for t in range(1, 5)]
    ell = next(vector for vector in candidates if all(vector))
    assert matvec(columns, ell) == (0, 0, 0)
    assert max(map(abs, ell)) <= 5 * max(map(abs, minors.values()))
    lambdas = [divisors[i] * ell[i] for i in ROWS]
    assert all(lambdas)
    assert all(lambdas[i] % divisors[i] == 0 for i in ROWS)
    assert sum(lambdas[i] * points[i][0] ** 2 for i in ROWS) == 0
    assert sum(lambdas[i] * points[i][0] * points[i][1] for i in ROWS) == 0
    assert sum(lambdas[i] * points[i][1] ** 2 for i in ROWS) == 0
    stresses = [lambdas[i] * norm(points[i]) for i in ROWS]
    assert sum(stresses) == 0
    # c_i*(P_i/bar P_i)=lambda_i*P_i^2, so no rational rounding occurs.
    squared = [power(p, 2) for p in points]
    assert sum(lambdas[i] * squared[i][0] for i in ROWS) == 0
    assert sum(lambdas[i] * squared[i][1] for i in ROWS) == 0


if __name__ == "__main__":
    check_orders()
    for exponents in ((0, 0, 0, 0, 0), (1, 1, 1, 1, 1), (0, 1, 2, 1, 3)):
        literal_source(exponents)
    print("PASS: all 310 cut/triple orders; common exponent90 and exact six-block quotients.")
    print("PASS: three literal 31-cut Gaussian sources, including shared correction primes.")
    print("PASS: actual gcd residues, five nonzero multipliers, and exact affine stresses.")
    print("No gradient of a common small invariant or unbounded endpoint family is constructed.")
