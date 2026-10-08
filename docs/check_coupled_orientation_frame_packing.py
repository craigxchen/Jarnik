"""Exact certificates for coupled_orientation_frame_packing.md (stdlib)."""

from itertools import combinations
from math import gcd, isqrt
from random import Random

from check_oriented_contact_lattice_rigidity import (
    add, complex_det, conj, crt, exact, iota, lin, matmul, mul, norm,
    power, q, upper, lower, valuation, check_system,
)


def frame(U):
    return ((U[0][0], U[1][0]), (U[0][1], U[1][1]))


def gram(x):
    return norm(x[0]), mul(conj(x[0]), x[1])[0], norm(x[1])


def det3(rows):
    a, b, c = rows
    return sum(a[j] * (b[(j+1) % 3] * c[(j+2) % 3]
                       - b[(j+2) % 3] * c[(j+1) % 3]) for j in range(3))


def identity_check(x, y):
    assert q(x) == q(y) == 1
    u = complex_det(x, y)
    v = complex_det(x, tuple(map(conj, y)))
    assert norm(v) - norm(u) == 4
    A, B, C = gram(x)
    a, b, c = gram(y)
    assert A*C - B*B == a*c - b*b == 1
    assert C*a + A*c - 2*B*b == norm(u) + 2
    assert norm(u) <= (A+C)*(a+c)
    assert A*A+B*B+C*C <= (A+C)**2
    return u, v


def contact_pair_check(x, y, contacts_x, contacts_y):
    u, v = identity_check(x, y)
    agreeing, flipped, D = (1, 0), (1, 0), 1
    for (g, lam), (g2, lam2) in zip(contacts_x, contacts_y):
        assert lam == lam2 and norm(g) == norm(g2)
        assert gcd(*lam) == 1
        exact(lin(lam[0], x[0], lam[1], x[1]), g)
        exact(lin(lam[0], y[0], lam[1], y[1]), g2)
        if g == g2:
            agreeing = mul(agreeing, g)
        else:
            assert conj(g) == g2
            flipped = mul(flipped, g)
        D *= norm(g)
    aa, ff = norm(agreeing), norm(flipped)
    assert aa*ff == D and gcd(aa, ff) == 1
    a = norm(u)
    assert a % aa == 0 and (a+4) % ff == 0
    assert a % D == crt([0, -4], [aa, ff]) % D
    cofactor = mul(exact(u, agreeing), exact(v, flipped))
    assert a*(a+4) == D*norm(cofactor)
    assert (a+2)**2 - D*norm(cofactor) == 4
    if a:
        assert a*(a+4) >= D
    return D


def random_identities():
    rng = Random(60392)
    frames = []
    for _ in range(100):
        U = ((1, 0), (0, 1))
        for j in range(5):
            U = matmul(U, (upper if j % 2 else lower)(rng.randrange(-8, 9)))
        frames.append(frame(U))
    for x, y in combinations(frames, 2):
        identity_check(x, y)
    print('4,950 determinant, Gram, and Hadamard identities pass.')


def shared_five_contacts():
    primes = [(2, 1), (3, 2), (4, 1), (5, 2), (6, 1),
              (5, 4), (7, 2), (6, 5), (8, 3), (8, 5)]
    data, residues, moduli = [], [[] for _ in range(5)], []
    for j, (edge, pi) in enumerate(zip(combinations(range(5), 2), primes)):
        exponent = 2 if j == 0 else 1
        p, root = norm(pi), iota(pi, exponent+1)
        median = int(sum(a in edge for a in range(3)) >= 2)
        cluster = [a for a in range(5) if int(a in edge) != median]
        outside = [a for a in range(5) if a not in cluster]
        allowed = [a for a in range(1, p) if a not in (2*root % p, -2*root % p)]
        local = {a: p**exponent*(a+1) for a in cluster}
        local.update({a: allowed[t] for t, a in enumerate(outside)})
        for a in range(5):
            residues[a].append(local[a])
        moduli.append(p**(exponent+1))
        data.append((pi, exponent, cluster))
    moduli.append(2)
    vectors = [(crt(row+[0], moduli), 1) for row in residues]
    assert len(set(vectors)) == 5
    frames, contact_sets = [], []
    for mask in [0, 1, 341, 682, 1023]:
        values, contacts = [], []
        for j, (pi, exponent, cluster) in enumerate(data):
            chosen = conj(pi) if (mask >> j) & 1 else pi
            p, pe = norm(pi), norm(pi)**exponent
            root, r = iota(chosen, exponent+1), vectors[cluster[0]][0]
            t = next(t for t in range(p)
                     if all(((vectors[a][0]-r)//pe+t) % p for a in cluster))
            values.append(-root-r+pe*t)
            contacts.append((power(chosen, exponent), vectors[cluster[0]]))
        B = crt(values+[0], moduli)
        x = ((1, 0), (B, 1))
        for (pi, exponent, cluster), (g, _) in zip(data, contacts):
            chosen = pi if g == power(pi, exponent) else conj(pi)
            for a, (r, s) in enumerate(vectors):
                value = lin(r, x[0], s, x[1])
                assert valuation(value, chosen) == (exponent if a in cluster else 0)
                assert valuation(value, conj(chosen)) == 0
        # The same source rows are normalized once, independently of mask.
        r0 = vectors[0][0]
        U = matmul(upper(B), ((r0, -1), (1, 0)))
        normalized = [(1, r0-r) for r, _ in vectors]
        check_system(*frame(U), normalized)
        frames.append(x)
        contact_sets.append(contacts)
    D = 1
    for pi, exponent, _ in data:
        D *= norm(pi)**exponent
    for i, j in combinations(range(len(frames)), 2):
        assert contact_pair_check(frames[i], frames[j], contact_sets[i], contact_sets[j]) == D
    for indices in combinations(range(len(frames)), 3):
        rows = [gram(frames[j]) for j in indices]
        assert det3(rows) % D == 0
    print('Five fixed source vectors, all ten clean contact assignments, five orientations, and full depth-two CRT/Pell/three-Gram divisibility pass.')
    print('The constructed five-vector frames have no asserted short profile height.')


def sharp_two_frame_example():
    x, y = ((1, 8), (1, 9)), ((-5, 9), (1, -2))
    pi, rho = (42, -103), (64, 91)
    for p in [norm(pi), norm(rho)]:
        assert all(p % d for d in range(2, isqrt(p)+1))
    contacts = [(pi, (1, 11729)), (rho, (1, 5690))]
    other = [(pi, (1, 11729)), (conj(rho), (1, 5690))]
    D = contact_pair_check(x, y, contacts, other)
    H2 = max(sum(map(norm, x)), sum(map(norm, y)))
    assert D == 153140621 and H2 == 147 and D > H2**3
    assert gram(x) != gram(y)
    print('Sharp general two-frame example passes: D=153140621 > H^6=147^3.')


if __name__ == '__main__':
    random_identities()
    shared_five_contacts()
    sharp_two_frame_example()
