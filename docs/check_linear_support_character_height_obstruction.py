"""Exact saturation, Gaussian heights, and four-row endpoint exclusions."""
from fractions import Fraction as Q
from itertools import combinations, product
from math import prod

from check_strict_obtuse_prime_box_countermodels import (
    nearby_split_primes, gaussian_prime, gaussian_mul, gaussian_conj,
    gaussian_norm, gaussian_gcd,
)


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def model(M):
    H = [[(-1)**(bin(a & x).count('1') % 2) for x in range(M)]
         for a in range(1, M)]
    assignments = [x % (M-1) for x in range(M)]
    F = [h.copy() for h in (H[a] for a in assignments)]
    for x in range(M):
        F[x][x] *= -1
    mult = [5-assignments.count(a) for a in range(M-1)]
    return H, F, assignments, mult


def character(c, H, F):
    return [Q(dot(c, h), 2) for h in H+F]


def finite_models():
    counts = 0
    for M in (8, 16, 32, 64):
        H, F, assignments, mult = model(M)
        S = H+F
        canonical = {tuple(s if s[0] == 1 else [-x for x in s]) for s in S}
        assert len(canonical) == 2*M-1
        assert all(abs(sum(s)) < M for s in S)
        for a, h in enumerate(H):
            assert sum(h) == 0
            for b, k in enumerate(H):
                assert dot(h, k) == (M if a == b else 0)
        for i, j in combinations(range(M), 2):
            assert sum(w*h[i]*h[j] for h, w in zip(H, mult)) + sum(f[i]*f[j] for f in F) <= -1
            c = [int(x == i)-int(x == j) for x in range(M)]
            v = character(c, H, F)
            assert all(x.denominator == 1 for x in v)
            assert M//2 <= sum(abs(x) for x in v) <= M+3
            u = v[:M-1]
            assert sum(abs(x) for x in u) == M//2
            assert all(Q(2*sum(u[a]*H[a][x] for a in range(M-1)), M) == c[x]
                       for x in range(M))
            counts += 1
        supports = {tuple(s[i] != s[j] for s in S) for i, j in combinations(range(M), 2)}
        assert len(supports) == M*(M-1)//2
    return counts


def saturation_fixtures():
    M = 8
    H, F, assignments, mult = model(M)
    count = 0
    # Every tested integer old-coordinate vector reconstructs a rational row
    # coefficient. New coordinates are integral exactly when twice that row is.
    for u in product((-1, 0, 1), repeat=M-1):
        c = [Q(2*sum(u[a]*H[a][x] for a in range(M-1)), M) for x in range(M)]
        assert sum(c) == 0
        v = character(c, H, F)
        assert v[:M-1] == list(u)
        assert all(v[M-1+x] == u[assignments[x]]-c[x]*H[assignments[x]][x]
                   for x in range(M))
        assert all(x.denominator == 1 for x in v) == all(x.denominator == 1 for x in c)
        if any(c) and all(x.denominator == 1 for x in c):
            assert sum(abs(x) for x in v) >= M*max(abs(x) for x in c)/2
        count += 1
    return count


def gaussian_fixture():
    M = 8
    H, F, assignments, mult = model(M)
    _, primes = nearby_split_primes(M)
    factors = [gaussian_prime(p) for p in primes]
    blocks = [(1, 0) for _ in H+F]
    norms = [1 for _ in blocks]
    # The first M physical columns are flipped at rows 0,...,M-1.
    for j, (p, pi) in enumerate(zip(primes, factors)):
        block = M-1+j if j < M else j % (M-1)
        blocks[block] = gaussian_mul(blocks[block], pi)
        norms[block] *= p
    assert all(gaussian_norm(z) == n for z, n in zip(blocks, norms))
    assert all(gaussian_norm(gaussian_gcd(z, gaussian_conj(z))) == 1 for z in blocks)
    N = prod(primes)
    assert prod(norms) == N
    pmin, pmax = min(primes), max(primes)
    # Stronger exact finite replacement for the near-log comparison in (7).
    assert pmin**(15*M) > N**3
    count = 0
    for i, j in combinations(range(M), 2):
        c = [int(x == i)-int(x == j) for x in range(M)]
        v = [int(x) for x in character(c, H, F)]
        beta = (1, 0)
        for z, power in zip(blocks, v):
            for _ in range(abs(power)):
                beta = gaussian_mul(beta, z if power > 0 else gaussian_conj(z))
        expected = prod(n**abs(power) for n, power in zip(norms, v))
        assert gaussian_norm(beta) == expected
        assert gaussian_norm(gaussian_gcd(beta, gaussian_conj(beta))) == 1
        assert expected**10 > N**3
        count += 1
    return count



def four_row_exclusions():
    count = 0
    for relabelled, Ms in ((False, (32, 64)), (True, (16, 32, 64))):
        for M in Ms:
            H, _, _, _ = model(M)
            base, p, q = (1, 3, 12) if relabelled else (2, 6, 24)
            rows = [base, base ^ p, base ^ q, base ^ p ^ q]
            c = [0]*M
            for x, value in zip(rows, (1, -1, 1, -1)):
                c[x] = value
            assert sum(c) == 0 and sum(abs(x) for x in c) == 4
            assert bin(p & p).count('1') % 2 == 0
            assert bin(q & q).count('1') % 2 == 0
            assert bin(p & q).count('1') % 2 == 0
            u = [dot(c, h)//2 for h in H]
            assert all(abs(x) in (0, 2) for x in u)
            assert sum(x != 0 for x in u) == M//4
            assigned = {x: (x-1 if x else M-1) if relabelled else x for x in range(M)}
            assert len(set(assigned.values())) == M
            inverse = {j: x for x, j in assigned.items()}
            physical = []
            reductions = []
            for j in range(5*(M-1)):
                a = j % (M-1)
                old = u[a]
                new = old-c[inverse[j]]*H[a][inverse[j]] if j in inverse else old
                physical.append(new)
                if old != new:
                    reductions.append(abs(old)-abs(new))
            assert reductions == [1]*4
            assert sum(abs(x) for x in physical) == 5*M//2-4
            assert max(abs(x) for x in physical) == 2
            if M == 32 and not relabelled:
                _, primes = nearby_split_primes(M)
                N = prod(primes)
                norm_beta = prod(prime**abs(power) for prime, power in zip(primes, physical))
                # Exact C=1 fixed-grid obstruction, without taking logarithms.
                assert 4*norm_beta**2 < N
                # Prime-block and Gaussian norm identities use disjoint primes.
                beta = (1, 0)
                for prime, power in zip(primes, physical):
                    if not power:
                        continue
                    z = gaussian_prime(prime)
                    if power < 0:
                        z = gaussian_conj(z)
                    for _ in range(abs(power)):
                        beta = gaussian_mul(beta, z)
                assert gaussian_norm(beta) == norm_beta
                assert gaussian_norm(gaussian_gcd(beta, gaussian_conj(beta))) == 1
                assert beta[0] and beta[1] and abs(beta[0]) != abs(beta[1])
            count += 1
    return count



def higher_coset_certificates():
    count = 0
    for M in (32, 64, 256, 512):
        H, _, _, _ = model(M)
        t = M.bit_length()-1
        for relabelled in (False, True):
            x0 = 1 if relabelled else 2
            dimension = t//2 if relabelled else (t-1)//2
            basis = [(3 if relabelled else 6) << (2*j) for j in range(dimension)]
            for chosen_dimension in range(1, dimension+1):
                V = [0]
                for p in basis[:chosen_dimension]:
                    V += [v ^ p for v in V]
                rows = [x0 ^ v for v in V]
                h = len(rows)
                coeff = [(-1)**(bin(x0 & v).count('1') % 2) for v in V]
                assert sum(coeff) == 0 and 0 not in rows
                u = [sum(c*col[x] for x, c in zip(rows, coeff))//2 for col in H]
                assert sum(abs(x) for x in u) == M//2
                assert sum(x != 0 for x in u) == M//h
                assert max(abs(x) for x in u) == h//2
                flips = {}
                for x, c in zip(rows, coeff):
                    j = (x-1 if x else M-1) if relabelled else x
                    a = j % (M-1)
                    assert u[a] and c*H[a][x] == (1 if u[a] > 0 else -1)
                    flips[j] = c*H[a][x]
                for copies in (5, 8, 17, 64):
                    r = copies*(M-1)
                    physical = [u[j % (M-1)]-flips.get(j, 0) for j in range(r)]
                    assert sum(abs(x) for x in physical) == copies*M//2-h
                    assert max(abs(x) for x in physical) <= h//2
                    assert Q(sum(abs(x) for x in physical), r) < 1
                    count += 1
    # Integer forms of the threshold inequalities used in (19).
    for C in (1, 2, 10, 10**6):
        threshold = 16
        while C**8 > 5**threshold:
            threshold += 1
        h = 1
        while h < threshold:
            h *= 2
        assert h**8 <= 5**h and C**8 <= 5**h
        assert 3*h*h*C*C < 8*5**(h//2)  # e<3 bounds RHS of (18).
    for t in range(10, 20):
        M = 2**t
        h = 2**((t-1)//2)
        copies = h+1
        r = copies*(M-1)
        assert 4*h*h >= M and 16*r*r >= M**3
    return count


if __name__ == '__main__':
    models = finite_models()
    saturation = saturation_fixtures()
    gaussian = gaussian_fixture()
    exclusions = four_row_exclusions()
    cosets = higher_coset_certificates()
    print(f'PASS: {models} exact characters at M=8,16,32,64; '
          f'{saturation} rational saturation fixtures; {gaussian} literal Gaussian heights; '
          'distinct classes/pair supports and negative weighted Gram; '
          f'{exclusions} four-row exclusions and literal M32 Gaussian grid obstruction; '
          f'{cosets} higher-coset/copy certificates and integer growth thresholds.')
