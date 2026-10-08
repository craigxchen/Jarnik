"""Finite support and stabilizer checks for rational CM cell phases."""

from itertools import product
from math import gcd, prod


def support_q(q):
    return {(0, 0)} | {(a, b) for a in range(1, q) for b in range(1, q)}


def add_q(x, y, q):
    return ((x[0] + y[0]) % q, (x[1] + y[1]) % q)


def product_support(primes):
    factors = [sorted(support_q(q)) for q in primes]
    return set(product(*factors))


def add_product(x, y, primes):
    return tuple(add_q(a, b, q) for a, b, q in zip(x, y, primes))


def check_single_prime(q):
    support = support_q(q)
    assert len(support) == q * q - 2 * q + 2
    assert len(support) % q == 2 % q
    assert {(a, (-b) % q) for a, b in support} == support

    stabilizer = {
        h for h in product(range(q), repeat=2)
        if {add_q(x, h, q) for x in support} == support
    }
    assert stabilizer == {(0, 0)}, (q, stabilizer)


def check_product(primes):
    support = product_support(primes)
    assert len(support) == prod(q * q - 2 * q + 2 for q in primes)
    factor_stabilizers = []
    for q in primes:
        single = support_q(q)
        factor_stabilizers.append({
            h for h in product(range(q), repeat=2)
            if {add_q(x, h, q) for x in single} == single
        })
    stabilizer = {
        h for h in product(*[list(product(range(q), repeat=2)) for q in primes])
        if all(h_part in allowed
               for h_part, allowed in zip(h, factor_stabilizers))
    }
    assert stabilizer == {tuple((0, 0) for _ in primes)}, (primes, stabilizer)


def c_vector(bits):
    return tuple(x - bits[-1] for x in bits[:-1])


def check_boolean_vectors(m):
    vectors = {}
    for bits in product((0, 1), repeat=m):
        vector = c_vector(bits)
        vectors.setdefault(vector, []).append(bits)
    assert vectors[(0,) * (m - 1)] == [(0,) * m, (1,) * m]
    for vector, labels in vectors.items():
        if vector != (0,) * (m - 1):
            assert len(labels) == 1
            opposite = tuple(-x for x in vector)
            assert len(vectors[opposite]) == 1
            assert vectors[opposite][0] == tuple(1 - x for x in labels[0])
            assert gcd(*[abs(x) for x in vector if x]) == 1


def check_contact_multipliers():
    # For pairwise coprime m_e, b_e=M/m_e has gcd one, so the common
    # translation kernel is trivial.
    moduli = (5, 13, 17)
    M = 1
    for m in moduli:
        M *= m
    b_values = [M // m for m in moduli]
    common = 0
    for b in b_values:
        common = gcd(common, b)
    assert common == 1


def main():
    for q in (5, 13):
        check_single_prime(q)
    check_product((5, 13))
    check_product((5, 13))
    for m in (2, 3, 4, 5):
        check_boolean_vectors(m)
    check_contact_multipliers()
    print("q=5,13 support sizes, negation stability, and trivial translation stabilizers pass.")
    print("Cartesian products inherit trivial stabilizers; Boolean valuation vectors are primitive and uniquely paired.")
    print("For moduli 5,13,17, gcd_e(M/m_e)=1, so the common translation kernel is trivial.")


if __name__ == "__main__":
    main()
