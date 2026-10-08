"""Exact cubic, integral basis, and raw-cofactor modular content checks."""

from functools import reduce
from itertools import combinations, permutations
from math import gcd, prod
from random import Random

from check_primitive_distance_smith_form import conj, divide, gp, mul, norm
from check_quartet_matching_gcd_cut_budget import ggcd


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def scale(k, a):
    return k * a[0], k * a[1]


def phi(v):
    a, b, c, d, e = v
    return b * c * e - a * d * sum(v)


def gphi(v):
    a, b, c, d, e = v
    return add(mul(mul(b, c), e), scale(-1, mul(mul(a, d), reduce(add, v))))


def det(a):
    n = len(a)
    return sum((-1) ** sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n))
               * prod(a[i][p[i]] for i in range(n)) for p in permutations(range(n)))


def cofactors(R):
    return [(-1) ** j * det([[row[k] for k in range(5) if k != j] for row in R])
            for j in range(5)]


def p_add(a, b, coefficient=1):
    out = dict(a)
    for exponent, value in b.items():
        out[exponent] = out.get(exponent, 0) + coefficient * value
        if not out[exponent]:
            del out[exponent]
    return out


def p_mul(a, b):
    out = {}
    for u, x in a.items():
        for v, y in b.items():
            exponent = tuple(s + t for s, t in zip(u, v))
            out[exponent] = out.get(exponent, 0) + x * y
    return {u: x for u, x in out.items() if x}


def difference(i, j):
    a, b = [0] * 5, [0] * 5
    a[i], b[j] = 1, 1
    return {tuple(a): -1, tuple(b): 1}


def check_polynomial():
    edges = [((0, 1), (2, 3)), ((0, 1), (3, 4)), ((0, 3), (1, 2)),
             ((1, 2), (3, 4)), ((1, 4), (2, 3))]
    w = [p_mul(difference(*a), difference(*b)) for a, b in edges]
    A, B, C, D, E = w
    total = reduce(p_add, w)
    assert not p_add(p_mul(p_mul(B, C), E), p_mul(p_mul(A, D), total), -1)
    pivots = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)]
    exponents = [tuple(int(k in p) for k in range(5)) for p in pivots]
    assert abs(det([[q.get(u, 0) for u in exponents] for q in w])) == 1


def matching_vector(rows):
    def bracket(i, j):
        a, b = rows[i], rows[j]
        return a[0] * b[1] - a[1] * b[0]
    edges = [((0, 1), (2, 3), 4), ((0, 1), (3, 4), 2),
             ((0, 3), (1, 2), 4), ((1, 2), (3, 4), 0),
             ((1, 4), (2, 3), 0)]
    return [scale(bracket(*a) * bracket(*b), rows[j]) for a, b, j in edges]


def divisible(z, divisor):
    numerator, denominator = mul(z, conj(divisor)), norm(divisor)
    return numerator[0] % denominator == numerator[1] % denominator == 0


def assert_lemma(R, v, H, kappa):
    assert gphi(v) == (0, 0)
    assert all(divisible(mul(kappa, reduce(add, (scale(a, x) for a, x in zip(row, v)))), H)
               for row in R)
    G = reduce(ggcd, v)
    M = divide(H, ggcd(H, mul(kappa, G)))
    c = cofactors(R)
    assert all(sum(a * b for a, b in zip(row, c)) == 0 for row in R)
    value = phi(c)
    assert divisible((value, 0), M)
    assert value % norm(M) == 0
    return value


def main():
    check_polynomial()
    rng = Random(20261002)
    v0 = matching_vector([(2 * j + 1, j * j + 1) for j in range(1, 6)])
    v0 = [divide(x, reduce(ggcd, v0)) for x in v0]
    assert gphi(v0) == (0, 0) and norm(reduce(ggcd, v0)) == 1
    checked = nonzero = 0
    for pi in [(2, 1), (3, 2), (4, 1)]:
        p = norm(pi)
        for e in range(1, 7):
            for a, b in [(0, 0), (1, 0), (0, 2), (2, 2)]:
                H, kappa = gp(pi, e), gp(pi, a)
                v = [mul(gp(pi, b), x) for x in v0]
                depth = max(0, e - a - b)
                modulus = p ** depth
                if depth:
                    M = gp(pi, depth)
                    root = -M[0] * pow(M[1], -1, modulus) % modulus
                    values = [(x + root * y) % modulus for x, y in v0]
                    pivot = next(j for j, x in enumerate(values) if gcd(x, p) == 1)
                else:
                    values, pivot = [0] * 5, 0
                R = []
                for _ in range(4):
                    row = [rng.randrange(-20, 21) for _ in range(5)]
                    if depth:
                        row[pivot] = -sum(row[j] * values[j] for j in range(5) if j != pivot)
                        row[pivot] = row[pivot] * pow(values[pivot], -1, modulus) % modulus
                    R.append(row)
                if checked % 4 == 0:
                    R[0] = [p * x for x in R[0]]
                nonzero += assert_lemma(R, v, H, kappa) != 0
                checked += 1
    # Multiple coprime oriented primes in the same modulus; the vector is an
    # actual matching evaluation at outside nodes 0,1,2,3,4.
    v_integer = [(x, 0) for x in (1, 1, 3, 1, 3)]
    for e in range(1, 5):
        H = mul(gp((2, 1), e), gp((3, 2), 5 - e))
        modulus = norm(H)
        R = []
        for _ in range(4):
            row = [rng.randrange(-20, 21) for _ in range(5)]
            row[0] = -sum(row[j] * v_integer[j][0] for j in range(1, 5)) % modulus
            R.append(row)
        nonzero += assert_lemma(R, v_integer, H, (1, 0)) != 0
        checked += 1
    R = [[5 * (int(j == i) if j < 4 else i + 1) for j in range(5)] for i in range(4)]
    assert_lemma(R, v_integer, (2, 1), (1, 0))
    c = cofactors(R)
    q = gcd(*c)
    assert q == 5 ** 4 and phi([x // q for x in c]) == 42
    print(f"verified universal cubic, integral basis, {checked} modular fixtures "
          f"({nonzero} nonzero cubic values), and primitive-cofactor content loss")


if __name__ == "__main__":
    main()
