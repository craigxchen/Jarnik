"""Exact shared edge labels, cycle roots, and all-edge affine contents."""

from fractions import Fraction as F
from itertools import combinations
from math import gcd, isqrt, prod

from check_gaussian_reflection_replacement import ggcd, gdivexact, gnorm, gmul, gconj
from check_direct_monic_remainder_local_obstruction import (
    remainder, prime_after, choose_units,
)


def binary_rank(vectors):
    basis = {}
    for vector in vectors:
        while vector:
            pivot = vector.bit_length()-1
            if pivot not in basis:
                basis[pivot] = vector
                break
            vector ^= basis[pivot]
    return len(basis)


def halfrows(points):
    N, anchor = gnorm(points[0]), points[0]
    result = []
    for point in points:
        dot = anchor[0]*point[0]+anchor[1]*point[1]
        cross = anchor[0]*point[1]-anchor[1]*point[0]
        if dot == -N:
            result.append((0, 1))
        else:
            g = gcd(N+dot, cross)
            result.append(((N+dot)//g, cross//g))
    return result


def audit(points):
    N = gnorm(points[0])
    rows = halfrows(points)
    fs = [gnorm(h) for h in rows]
    data, Gplus = {}, 0
    common = points[0]
    for point in points[1:]:
        common = ggcd(common, point)
    primitive = [gdivexact(z, common) for z in points]
    assert gnorm(primitive[0]) % 2 == 1
    orientation = len({z[0] % 2 for z in primitive})
    for i, j in combinations(range(len(points)), 2):
        dot = sum(a*b for a, b in zip(points[i], points[j]))
        c = N-dot
        a, b = rows[i]
        x, y = rows[j]
        delta, A = a*y-x*b, a*x+b*y
        assert F(2*N*delta*delta, fs[i]*fs[j]) == c > 0
        assert F(2*N*A*A, fs[i]*fs[j]) == 2*N-c
        Gplus = gcd(Gplus, 2*N-c)
        data[i, j] = c, delta, A
    assert Gplus == gnorm(common)*(2 if orientation == 1 else 1)
    reference = next(iter(data.values()))[0]
    Geff = 2*N-reference
    for c, delta, A in data.values():
        Geff = gcd(Geff, c-reference)
    assert Geff == Gplus
    cycles = 0
    for a, b, c, d in combinations(range(min(len(points), 8)), 4):
        for walk in ((a, b, c, d), (a, b, d, c), (a, c, b, d)):
            edges = [tuple(sorted((walk[i], walk[(i+1) % 4]))) for i in range(4)]
            entries = [data[e] for e in edges]
            denominator = prod(fs[i] for i in walk)
            root_c = F((2*N)**2*prod(t[1] for t in entries), denominator)
            root_P = F((2*N)**2*prod(t[2] for t in entries), denominator)
            assert root_c.denominator == root_P.denominator == 1
            assert root_c*root_c == prod(t[0] for t in entries)
            assert root_P*root_P == prod(2*N-t[0] for t in entries)
            cycles += 1
    return len(data), cycles


def rank_fixtures():
    cases = kernels = 0
    for m in range(3, 8):
        hs, fs, t = [], [], 2
        for _ in range(m-1):
            h, f = (1, t), 1+t*t
            assert f % 2 and isqrt(f)**2 != f
            assert all(gcd(f, old) == 1 for old in fs)
            assert gnorm(ggcd(h, gconj(h))) == 1
            hs.append(h)
            fs.append(f)
            t = 2*prod(fs)
        anchor = (1, 0)
        for h in hs:
            anchor = gmul(anchor, gconj(h))
        points = [anchor] + [gdivexact(gmul(anchor, h), gconj(h)) for h in hs]
        assert len(set(points)) == m and all(gnorm(z) == prod(fs) for z in points)
        audit(points)
        # Basis consists of independent odd f_i classes and the prime 2.
        w = ((1 << m)-1) | (1 << m)
        vertices = [0] + [1 << j for j in range(m-1)]
        edges = list(combinations(range(m), 2))
        labels = [vertices[i] ^ vertices[j] ^ w for i, j in edges]
        assert binary_rank(labels) == m
        if m <= 5:
            for mask in range(1 << len(edges)):
                vector, degree, size = 0, [0]*m, 0
                for index, (i, j) in enumerate(edges):
                    if mask >> index & 1:
                        vector ^= labels[index]
                        degree[i] ^= 1
                        degree[j] ^= 1
                        size += 1
                assert (vector == 0) == (size % 2 == 0 and not any(degree))
                kernels += 1
        cases += 1
    return cases, kernels


def local_cycle_models():
    cases = 0
    for k in range(2, 9):
        p = prime_after(8*k)
        _, units = choose_units(p, k, 2)
        for e in (1, 2, 3):
            N, modulus = p**e, p**(e+2)
            ts = list(range(1, 2*k+1))
            ts[1], ts[-1] = -2*units[0] % modulus, -2*units[1] % modulus
            coordinates = [(N, 1)] + [(t, N*pow(t, -1, modulus) % modulus)
                                      for t in ts[1:]]
            assert all((a*b-N) % modulus == 0 for a, b in coordinates)
            roots = []
            for index in range(2*k):
                a, b = coordinates[index]
                c, d = coordinates[(index+1) % (2*k)]
                roots.append((N-(a*d+b*c)*pow(2, -1, modulus)) % modulus)
            assert sum(c % p != 0 for c in roots) == 2
            assert remainder(roots, 2*N) % p != 0
            cases += 1
    roots = [2, 40, 9, 5]
    assert prod(roots) == 60**2 and prod(130-c for c in roots) == 13200**2
    assert remainder(roots, 130) == -41433616
    assert remainder(roots, 130) % 5 == 4
    return cases



def graph_audit(points, edges):
    """Primitive binary tuples; no arc restriction for the exact formulas."""
    N, m = gnorm(points[0]), len(points)
    common = points[0]
    for z in points[1:]:
        common = ggcd(common, z)
    assert gnorm(common) == 1 and N % 2
    assert all(gcd(*z) == 1 for z in points)
    adjacency = [[] for _ in points]
    G = L = 0
    for i, j in edges:
        G = gcd(G, N + sum(a*b for a, b in zip(points[i], points[j])))
        L = gcd(L, points[i][0]+points[j][0], points[i][1]+points[j][1])
        adjacency[i].append(j)
        adjacency[j].append(i)
    colors, queue, bipartite = {0: 0}, [0], True
    for i in queue:
        for j in adjacency[i]:
            if j not in colors:
                colors[j] = 1-colors[i]
                queue.append(j)
            elif colors[j] == colors[i]:
                bipartite = False
    assert len(colors) == m and G > 0 and L > 0
    assert gcd(G, N) == 1
    uniform = len({z[0] % 2 for z in points}) == 1
    if uniform:
        assert L % 2 == 0 and 2*G == L*L
        K = L//2
    else:
        assert L % 2 == 1 and G == L*L
        K = L
    certificates = 0
    if not bipartite:
        assert K == 1 and G <= 2
    else:
        for color in (0, 1):
            part = [i for i in range(m) if colors[i] == color]
            short = True
            for i, j in combinations(part, 2):
                a, b = halfrows([points[i], points[j]])[1]
                f = a*a+b*b
                n = f if f % 2 else f//2
                assert (2*N) % f == 0 and N % n == 0 and b % K == 0
                c = N-sum(x*y for x, y in zip(points[i], points[j]))
                if c*c <= N:
                    assert b**4*N <= n*n
                else:
                    short = False
            if short and len(part) >= 2:
                assert K**(4*(len(part)-1)) <= N
                certificates += 1
    return int(K > 1), certificates


def sparse_graphs():
    cases = nontrivial = certificates = 0
    for N in (5, 13, 25, 65, 85):
        points = []
        for x in range(-isqrt(N), isqrt(N)+1):
            y = isqrt(N-x*x)
            if x*x+y*y == N and gcd(x, y) == 1:
                points.extend([(x, y)] + ([(x, -y)] if y else []))
        for subset in combinations(points, 3):
            common = ggcd(ggcd(subset[0], subset[1]), subset[2])
            if gnorm(common) != 1:
                continue
            for edges in (((0, 1), (1, 2)), ((0, 1), (0, 2)),
                          ((0, 1), (1, 2), (0, 2))):
                n, c = graph_audit(subset, edges)
                cases, nontrivial, certificates = cases+1, nontrivial+n, certificates+c
        for m in (4, 5, 6):
            subset = points[:m]
            common = subset[0]
            for z in subset[1:]:
                common = ggcd(common, z)
            if gnorm(common) != 1:
                continue
            graphs = [list(combinations(range(m), 2)),
                      [(i, (i+1) % m) for i in range(m)],
                      [(i, j) for i in range(m//2) for j in range(m//2, m)]]
            for edges in graphs:
                n, c = graph_audit(subset, edges)
                cases, nontrivial, certificates = cases+1, nontrivial+n, certificates+c
    # Nontrivial ramified contribution on a connected spanning binary graph.
    assert graph_audit([(1, 8), (-1, 8), (7, 4)], [(0, 1), (0, 2)])[0] == 1
    # Repeated-factor cancellation can destroy connected-graph normalization.
    cycle = [(-8, -1), (-7, -4), (-4, -7), (-1, -8), (4, -7), (7, -4)]
    graph_audit(cycle, [(i, (i+1) % 6) for i in range(6)])
    deficits = [65-sum(a*b for a, b in zip(cycle[i], cycle[(i+1) % 6]))
                for i in range(6)]
    assert deficits == [5, 9, 5, 13, 9, 117]
    retained = [c for c in set(deficits) if deficits.count(c) % 2]
    assert sorted(retained) == [13, 117]
    assert gcd(*(130-c for c in deficits)) == 1
    assert gcd(*(130-c for c in retained)) == 13
    assert nontrivial and certificates
    return cases, nontrivial, certificates


def main():
    edges = cycles = 0
    for N in (5, 9, 13, 18, 25, 50, 65, 100, 125):
        points = []
        for x in range(-isqrt(N), isqrt(N)+1):
            y = isqrt(N-x*x)
            if x*x+y*y == N:
                points.append((x, y))
                if y:
                    points.append((x, -y))
        e, c = audit(points)
        edges += e
        cycles += c
        # Proper subsets include single coordinate-parity orientations.
        for subset in combinations(points, 3):
            audit(subset)
    ranks, kernels = rank_fixtures()
    local = local_cycle_models()
    graphs, nontrivial, certificates = sparse_graphs()
    print(f'PASS: {edges} literal edge dictionaries; {cycles} cycle square roots; '
          f'{ranks} primitive sharp-rank fixtures; {kernels} exact formal kernels; '
          f'{local} shared-local-circle remainder models; actual quartic remainder; '
          'all full-circle and triple affine-gcd checks, including scaled norms; '
          f'{graphs} sparse graphs ({nontrivial} nontrivial K), '
          f'{certificates} exact short-pair Gram certificates; cancellation fixture.')


if __name__ == '__main__':
    main()
