"""Exact projection and even-squareclass checks on actual mapped tuples."""

from fractions import Fraction
from itertools import combinations
from math import gcd, isqrt, prod

from check_mobius_conductor_transfer import fixtures
from check_mobius_reciprocal_stretch_grid import clips, content, gram, realize
from check_gaussian_reflection_replacement import gnorm
from check_isometric_virtual_anchor import factor


def check_relation_system(E, Eprime, delta, T, solutions):
    """Check root collisions and every available even square-product relation."""
    H = 4 * delta**2 * E * Eprime
    B = 2 * delta * E
    Z = 2 * E * Eprime * T
    roots = {}
    for n, y in solutions:
        assert n > 0 and H % n == 0
        assert y*y == 2*T*E*n - n*n - B*B
        c = Eprime*n + E*(H//n)
        assert Eprime*y*y == n*(Z-c)
        assert 0 < c <= H*(E+Eprime) <= H*H
        if c in roots:
            old_n, _ = roots[c]
            assert old_n == n or old_n*n == B*B
        else:
            roots[c] = n, y
    labels = list(roots.values())
    for (n, _), (v, _) in combinations(labels, 2):
        assert n*v != B*B

    primes = list(factor(H))
    vectors = []
    for n, _ in labels:
        bits = 1 << len(primes)
        residual = n
        for j, p in enumerate(primes):
            exponent = 0
            while residual % p == 0:
                residual //= p
                exponent += 1
            bits |= (exponent % 2) << j
        assert residual == 1
        vectors.append(bits)

    relations = 0
    short_relation = False
    for size in range(2, len(labels)+1, 2):
        for indices in combinations(range(len(labels)), size):
            product_n = prod(labels[j][0] for j in indices)
            square_root = isqrt(product_n)
            xor = 0
            for j in indices:
                xor ^= vectors[j]
            assert (xor == 0) == (square_root*square_root == product_n)
            if xor:
                continue
            cs = [Eprime*labels[j][0] + E*(H//labels[j][0]) for j in indices]
            assert len(set(cs)) == size
            square = prod(Z-c for c in cs)
            rational_root = Fraction(Eprime**(size//2)
                                     * prod(labels[j][1] for j in indices), square_root)
            assert rational_root*rational_root == square
            assert rational_root.denominator == 1
            assert square >= 0 and isqrt(square)**2 == square
            k = size//2
            assert Z <= (8*k*H*H)**(2*k+1)
            relations += 1
            short_relation |= size <= len(primes)+2
    if len(labels) > len(primes)+1:
        assert short_relation
    return relations


def main():
    cases = identities = relations = 0
    row_sets = fixtures() + [[(1, 0), (-26, 7), (3, -1), (-17, -1)]]
    for a in range(1, 4):
        for d in range(1, 4):
            for b in range(-3, 4):
                if gcd(a, b, d) != 1 or (a == d and b == 0):
                    continue
                delta = a*d
                T, W = gram(a, b, d)
                for rows in row_sets:
                    source, N = realize(rows)
                    mapped = [(a*x+b*y, d*y) for x, y in rows]
                    _, Nprime = realize(mapped)
                    C, g = content(W, source[0])
                    Q = prod(p**e for p, e in clips(a, b, d, source[0], N).items())
                    assert g % Q == N % Q == Nprime % Q == 0
                    E, Eprime = N//Q, Nprime//Q
                    solutions = []
                    for h, mh, z in zip(rows, mapped, source):
                        re = C[0]*z[0] + C[1]*z[1]
                        im = C[0]*z[1] - C[1]*z[0]
                        assert re % Q == im % Q == 0
                        n = T*E + re//Q
                        y = im//Q
                        assert n == 2*E*Fraction(gnorm(mh), gnorm(h))
                        assert (n-T*E)**2 + y*y == (T*T-4*delta*delta)*E*E
                        solutions.append((n, y))
                        identities += 1
                    relations += check_relation_system(E, Eprime, delta, T, solutions)
                    cases += 1

    # A literal squareclass relation with distinct roots, used only as an
    # algebraic fixture (no endpoint map existence is asserted).
    # E=E'=1, delta=3, T=31, H=36: labels 1 and 4 give Y=5 and 14.
    assert check_relation_system(1, 1, 3, 31, [(1, 5), (4, 14)]) == 1
    # The old six-point Pell degeneracy must retain its repeated root.
    for s, t in [(1, 1), (5, 7), (29, 41), (169, 239)]:
        assert t*t-2*s*s == -1
        assert check_relation_system(1, 1, 1, s*s+2,
                                     [(1, t), (2, 2*s), (4, 2*t)]) == 0
    print(f"PASS: {cases} actual mapped tuples; {identities} integer square projections; "
          f"{relations} actual even-squareclass relations; distinct-root and Pell controls")


if __name__ == "__main__":
    main()
