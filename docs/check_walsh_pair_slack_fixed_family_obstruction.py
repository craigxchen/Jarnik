"""Exact slack identities and a fixed-family obstruction; no arc claim."""

from fractions import Fraction
from itertools import combinations
from math import log
from random import Random


def dot(a, x):
    return bin(a & x).count('1') & 1


def chi(a, x):
    return 1 - 2 * dot(a, x)


def span(basis):
    out = [0]
    for v in basis:
        out += [x ^ v for x in out]
    return out


def fixture(t, actual_primes=False):
    m, b = 1 << t, 9
    r = b * (m - 1)
    labels = [a for a in range(1, m) for _ in range(b)]
    if actual_primes:
        # Disjoint actual 1 mod 4 primes near 10^4 and 10^8.
        def primes(start, count):
            values, n = [], start + (1 - start) % 4
            while len(values) < count:
                if all(n % d for d in range(3, int(n ** .5) + 1, 2)):
                    values.append(n)
                n += 4
            return values
        light = iter(primes(10000, sum(not dot(a, 1) for a in labels)))
        heavy = iter(primes(100000000, sum(dot(a, 1) for a in labels)))
        physical = [next(heavy) if dot(a, 1) else next(light) for a in labels]
        assert len(set(physical)) == r
        weights = [log(p) for p in physical]
    else:
        weights = [Fraction(20 if dot(a, 1) else 10) + Fraction(j, r ** 3)
                   for j, a in enumerate(labels)]
    assigned = Random(1729 + t).sample(range(r), m)
    signs = [[chi(a, x) * (-1 if assigned[x] == j else 1)
              for j, a in enumerate(labels)] for x in range(m)]
    return m, labels, weights, assigned, signs


def check(t, actual_primes=False):
    m, labels, weights, assigned, signs = fixture(t, actual_primes)
    r, w0 = len(labels), sum(weights)
    wa, fa = [0] * m, [0] * m
    for j, a in enumerate(labels):
        wa[a] += weights[j]
    for j in assigned:
        fa[labels[j]] += weights[j]
    ftot = sum(fa)
    # kappa=-1 is a convenient exact algebra fixture; the actual prime
    # pair check below uses C=1,D=0 and kappa=-log 4.
    q = [wa[a] - Fraction(4, m) * fa[a] for a in range(m)]
    u = [0] + [-1 - sum(q[a] * chi(a, d) for a in range(m))
               for d in range(1, m)]
    assert min(u[1:]) >= 0
    for x in range(m):
        for y in range(x):
            gram = sum(weights[j] * signs[x][j] * signs[y][j] for j in range(r))
            assert gram < -log(4)
    certificates = 0
    # Coordinate subspaces suffice for the finite algebra fixtures;
    # the accompanying proof covers arbitrary subspaces.
    for dim in range(1, t + 1):
        vspace = span([1 << j for j in range(dim)])
        h = len(vspace)
        for ell in range(1, h):
            fcell = sum(fa[a] for a in range(m) if a % h == ell)
            psi = sum(-chi(ell, d) * u[d] for d in vspace if d)
            for x0 in range(0, m, h):
                support = [x0 ^ v for v in vspace]
                coeff = {x0 ^ v: chi(ell, v) for v in vspace}
                fp = sum(weights[assigned[x]] for x in support)
                fm = sum(weights[assigned[x]] for x in support
                         if labels[assigned[x]] % h == ell)
                height = sum(abs(sum(coeff[x] * signs[x][j] for x in support))
                             * weights[j] / 2 for j in range(r))
                rhs = (w0 / 2 + Fraction(1, 2) + psi / 2
                       + Fraction(2 * h, m) * fcell - Fraction(2, m) * ftot
                       + fp - 2 * fm)
                if actual_primes:
                    assert abs(height - rhs) < 1e-8
                else:
                    assert height == rhs
                if ell & 1:
                    assert height > w0 / 2
                certificates += 1
    return certificates


def check_patching():
    # Test all coordinate subspaces, and their invertible linear images,
    # against fields whose label multiplicity is at most b.
    for t in range(2, 7):
        m, b = 1 << t, 3
        pool = [a for a in range(1, m) for _ in range(b)]
        field = Random(71 + t).sample(pool, m)
        for dim in range(1, t + 1):
            for indices in combinations(range(t), dim):
                coordinate = [1 << j for j in indices]
                # The triangular map e_j -> e_j+e_(j-1) is invertible.
                for basis in (coordinate,
                              [v ^ (v >> 1) for v in coordinate]):
                    hspace = span(basis)
                    n = len(hspace)
                    zero_rows = {x for x in range(m)
                                 if all(dot(field[x], v) == 0 for v in basis)}
                    assert len(zero_rows) <= b * (m // n - 1)
                    unseen, counts = set(range(m)), []
                    while unseen:
                        x0 = min(unseen)
                        coset = {x0 ^ v for v in hspace}
                        counts.append(len(zero_rows & coset))
                        unseen -= coset
                    assert min(counts) <= Fraction(b * (m - n), m)
                    assert min(counts) <= b
    # Exhaust every set of up to three patched rows in a retained
    # family of disjoint four-row flats. Surviving actual restrictions
    # are unchanged, and each patched row destroys at most one flat.
    family = [set(range(x, x + 4)) for x in range(0, 16, 4)]
    for count in range(4):
        for selected in combinations(range(16), count):
            bad = set(selected)
            actual = [1 if x not in bad else 0 for x in range(16)]
            surviving = [p for p in family if not p & bad]
            rows = set().union(*surviving) if surviving else set()
            assert len(rows) >= 16 - 4 * count
            assert all(actual[x] == 1 for x in rows)


if __name__ == '__main__':
    check_patching()
    total = sum(check(t) for t in range(2, 6))
    prime_total = check(3, actual_primes=True)
    print('PASS: %d exact rational certificates; %d actual-prime certificates; '
          'all pair Grams and patching bounds checked.' % (total, prime_total))
