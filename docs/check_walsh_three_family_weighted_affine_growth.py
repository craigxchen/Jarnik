"""Exact three-family affine identities, capacity, and literal source phases."""

from fractions import Fraction as Q
from itertools import combinations
from math import atan2, log, pi, prod

from check_walsh_arbitrary_binary_bilinear_isotropic_bound import images, parity
from check_walsh_ordinary_dot_weighted_phase_growth import split_primes
from check_walsh_symplectic_weighted_phase_growth import gaussian_product
from check_strict_obtuse_prime_box_countermodels import (
    gaussian_prime, gaussian_mul, gaussian_conj, gaussian_norm, gaussian_gcd,
)


def span(basis):
    out = {0}
    for v in basis:
        out |= {x ^ v for x in out}
    return out


def transpose_image(columns, x):
    return sum(parity(a & x) << j for j, a in enumerate(columns))


def find_direction(columns, d):
    im, M = images(columns), 1 << len(columns)
    def visit(basis, first):
        V = span(basis)
        if len(basis) == d:
            return V
        for u in range(first, M):
            if u in V or parity(u & im[u]):
                continue
            if any(parity(u & im[v]) or parity(v & im[u]) for v in V):
                continue
            enlarged = V | {v ^ u for v in V}
            if len({transpose_image(columns, v) for v in enlarged}) != len(enlarged):
                continue
            result = visit(basis+[u], u+1)
            if result is not None:
                return result
        return None
    result = visit([], 1)
    assert result is not None
    return result


def cosets(space, direction):
    unseen, result = set(space), []
    while unseen:
        x = min(unseen)
        P = {x ^ v for v in direction}
        assert P <= unseen
        unseen -= P
        result.append(P)
    return result


def difference(P, T):
    assert not (P & T)
    return {**{x: 1 for x in P}, **{x: -1 for x in T}}


def setup(columns, offset, edited):
    t, im = len(columns), images(columns)
    M = 1 << t
    k = t-(len(set(im)).bit_length()-1)
    d = max(0, (t-2)//2-k)
    V = find_direction(columns, d)
    h = len(V)
    assert h >= 2
    x0 = next(x for x in range(M)
              if all(parity(v & (im[x] ^ offset)) == 0 for v in V))
    residual = im[x0] ^ offset
    assert all(not parity(v & residual) for v in V)
    H = {x for x in range(M) if all(not parity(v & im[x]) for v in V)}
    n = M//(h*h)
    assert V <= H and len(H) == n*h and n > 1
    actual = {x: (im[x] ^ offset) or 1 for x in range(M)}
    for x, label in edited.items():
        assert 0 < label < M
        actual[x] = label
    bad_original = {x for x in range(M) if actual[x] != (im[x] ^ offset)}
    K = len(bad_original)
    predicted_translated = {y: im[y] ^ residual for y in range(M)}
    actual_translated = {x ^ x0: a for x, a in actual.items()}
    bad = {x ^ x0 for x in bad_original}
    assert bad == {y for y in range(M)
                   if actual_translated[y] != predicted_translated[y]}
    counts = {a: list(actual.values()).count(a) for a in range(1, M)}
    copies = max(5, max(counts.values()))
    physical = [(a, (-1)**parity(a & x0), None)
                for copy in range(copies) for a in range(1, M)]
    used = {a: 0 for a in range(1, M)}
    for x, a in actual_translated.items():
        j = used[a]*(M-1)+a-1
        used[a] += 1
        physical[j] = (a, physical[j][1], x)
    assert len({x for _, _, x in physical if x is not None}) == M
    # Image-capacity inequality uses actual agreement, not label injectivity.
    assert copies*len(set(im)) >= M-K
    if k:
        assert offset not in set(im)  # Deficient-rank fixture cannot solve Lx0=offset.
    groups = {}
    groups['g'] = []
    for P in cosets(range(M), V):
        if P <= H:
            continue
        xp = min(P)
        label = predicted_translated[xp]
        c = {x: (-1)**parity(label & (x ^ xp)) for x in P}
        assert sum(c.values()) == 0 and len(c) == h
        groups['g'].append(c)
    groups['p'] = [difference(P, T) for P, T in combinations(cosets(H, V), 2)]
    groups['f'] = [difference(P, T) for P, T in combinations(cosets(range(M), H), 2)]
    assert len(groups['g']) == n*(h-1)
    Vperp = {a for a in range(M) if all(not parity(a & v) for v in V)}
    Hperp = {a for a in range(M) if all(not parity(a & x) for x in H)}
    blocks = [set(range(1, M))-Vperp, Vperp-Hperp, Hperp-{0}]
    assert set.union(*blocks) == set(range(1, M))
    return dict(t=t, M=M, h=h, n=n, k=k, K=K, H=H, bad=bad,
                copies=copies, physical=physical, groups=groups, blocks=blocks)


def character(c, physical, flips=True):
    result = []
    for a, orientation, row in physical:
        raw = sum(value*orientation*((-1)**parity(a & x))
                  for x, value in c.items())
        assert raw % 2 == 0
        v = raw//2
        if flips and row in c:
            v -= c[row]*orientation*((-1)**parity(a & row))
        result.append(v)
    return result


def exact_fixture(data):
    M, h, n = data['M'], data['h'], data['n']
    H, bad, physical = data['H'], data['bad'], data['physical']
    weights = [1+(19*j*j+23*j)%211 for j in range(len(physical))]
    W0 = sum(weights)
    f = {x: w for w, (_, _, x) in zip(weights, physical) if x is not None}
    Fout = sum(w for x, w in f.items() if x not in H)
    FH = sum(w for x, w in f.items() if x in H)
    Fbad = sum(w for x, w in f.items() if x not in H and x in bad)
    Ftot = sum(f.values())
    Wa = {a: sum(w for w, (label, _, _) in zip(weights, physical) if label == a)
          for a in range(1, M)}
    Fa = {a: sum(w for w, (label, _, row) in zip(weights, physical)
                 if label == a and row is not None) for a in range(1, M)}
    coefficients = [Q(h, 2*(h-1)), Q(h*n, 2*(n-1)), Q(M, 2*(h-1))]
    lambdas = [Q(h-1, h), Q(n-1, h*n), Q(h-1, M)]
    assert sum(lambdas) == Q(M-1, M)
    assert [a*b for a, b in zip(coefficients, lambdas)] == [Q(1, 2)]*3
    averages, count = [], 0
    for index, name in enumerate(('g', 'p', 'f')):
        certs = data['groups'][name]
        base_sum = [0]*len(physical)
        heights = []
        for c in certs:
            base, changed = character(c, physical, False), character(c, physical)
            assert sum(c.values()) == 0
            assert sum(abs(v) for v in base) == data['copies']*M//2
            assert sum(abs(v) for v in changed) >= data['copies']*M//2-len(c) > 0
            for j, (a, z) in enumerate(zip(base, changed)):
                assert abs(abs(z)-abs(a)) <= 1
                base_sum[j] += abs(a)
            if name == 'g':
                for j, (_, _, row) in enumerate(physical):
                    if row in c and row not in bad:
                        assert abs(changed[j]) == abs(base[j])-1
            heights.append(sum(abs(v)*w for v, w in zip(changed, weights)))
            count += 1
        for total, (a, _, _) in zip(base_sum, physical):
            assert Q(total, len(certs)) == coefficients[index]*int(a in data['blocks'][index])
        actual = Q(sum(heights), len(heights))
        block_weight = sum(Wa[a] for a in data['blocks'][index])
        correction = [Q(-Fout+2*Fbad, n*(h-1)), Q(2*FH, n), Q(2*Ftot, h)][index]
        assert actual <= coefficients[index]*block_weight+correction
        averages.append(actual)
    mixed = sum(a*b for a, b in zip(lambdas, averages))
    upper = Q(W0, 2)+Q(-h*Fout+2*h*Fbad, M)+Q(2*(n-1)*FH, h*n*n)+Q(2*(h-1)*Ftot, M*h)
    assert mixed <= upper
    # Check the exact pair-Gram averaging and cone inversion for arbitrary weights.
    signs = [[orientation*((-1)**(parity(a & x)+int(row == x)))
              for a, orientation, row in physical] for x in range(M)]
    gram = {(x, y): sum(w*s*t for w, s, t in zip(weights, signs[x], signs[y]))
            for x in range(M) for y in range(x+1, M)}
    kap = max(gram.values())
    q = {a: Q(Wa[a])-Q(4*Fa[a], M) for a in Wa}
    assert min(q.values()) > 0
    hat = {d: sum(w*((-1)**parity(a & d)) for a, w in q.items())
           for d in range(1, M)}
    for d in range(1, M):
        actual = Q(sum(gram[tuple(sorted((x, x ^ d)))] for x in range(M)), M)
        assert actual == hat[d] <= kap
    u = {d: kap-value for d, value in hat.items()}
    assert sum(u.values()) == sum(q.values())+(M-1)*kap
    for a in q:
        assert q[a]+kap == Q(2, M)*sum(value for d, value in u.items() if parity(a & d))
        assert Wa[a] <= Q(2*W0+(M-2)*kap, M-4)
    Bmax = Q(2*W0+(M-2)*max(0, kap), M-4)
    assert FH <= (M//h)*Bmax and Fbad <= data['K']*Bmax
    # Full quantitative q algebra, using exact rational bounds before logarithms.
    s = Q(h*data['K'], M)
    assert Q(47, 6)+Q(16, 3)*s <= 8*(1+s)
    assert Q(2*M*(n-1), n)*Bmax <= Q(16, 3)*W0+Q(7, 3)*M*max(0, kap)
    assert 2*h*data['K']*Bmax <= Q(16, 3)*s*W0+Q(7, 3)*s*M*max(0, kap)
    return count


def angle(z):
    shift = max(abs(z[0]).bit_length(), abs(z[1]).bit_length())-500
    if shift > 0:
        z = (z[0] >> shift, z[1] >> shift)
    return atan2(z[1], z[0]) % (2*pi)


def literal_fixture(data):
    physical, M = data['physical'], data['M']
    primes = split_primes(len(physical))
    factors = [gaussian_prime(p) for p in primes]
    factors = [gaussian_conj(z) if j % 3 == 1 else z for j, z in enumerate(factors)]
    common, units = (2, 1), [(1, 0), (0, 1), (-1, 0), (0, -1)]
    points = []
    for x in range(M):
        raw = gaussian_product(z if orientation*((-1)**(parity(a & x)+int(row == x))) == 1
                               else gaussian_conj(z)
                               for z, (a, orientation, row) in zip(factors, physical))
        points.append(gaussian_mul(gaussian_mul(common, units[x % 4]), raw))
    W0, D = sum(map(log, primes)), log(5)
    assert all(gaussian_norm(z) == 5*prod(primes) for z in points)
    angles = sorted(map(angle, points))
    gaps = [angles[i+1]-angles[i] for i in range(M-1)]+[angles[0]+2*pi-angles[-1]]
    delta = 2*pi-max(gaps)
    logC = log(delta)+(W0+D)/4+0.01
    count = 0
    for name in ('g', 'p', 'f'):
        # Check the first and last members, including possible origin/exception supports.
        chosen = data['groups'][name]
        for c in [chosen[0], chosen[-1]]:
            v = character(c, physical)
            beta = gaussian_product(z if exponent > 0 else gaussian_conj(z)
                                    for z, exponent in zip(factors, v)
                                    for _ in range(abs(exponent)))
            assert gaussian_norm(beta) == prod(p**abs(e) for p, e in zip(primes, v)) > 1
            assert gaussian_norm(gaussian_gcd(beta, gaussian_conj(beta))) == 1
            assert beta[0] and beta[1] and abs(beta[0]) != abs(beta[1])
            positive = gaussian_product(points[x] for x, value in c.items() if value == 1)
            negative = gaussian_product(points[x] for x, value in c.items() if value == -1)
            unit = units[sum(x*value for x, value in c.items()) % 4]
            assert gaussian_mul(positive, gaussian_conj(beta)) == gaussian_mul(
                gaussian_mul(unit, beta), negative)
            height = sum(abs(e)*log(p) for p, e in zip(primes, v))
            assert height >= (W0+D)/2+log(8)-2*log(len(c))-2*logC-1e-9
            count += 1
    # These constants refer to the fixture's actual containing arc; no fixed-C cluster claim.
    kplus = max(0, 4*logC-log(4))
    LC = 2+((7/3)*kplus+2*max(0, logC))/log(16)
    q = 1+data['h']*data['K']/M
    f = {x: log(p) for p, (_, _, x) in zip(primes, physical) if x is not None}
    Fout = sum(w for x, w in f.items() if x not in data['H'])
    assert data['h']*Fout <= 8*q*W0+LC*q*M*log(M)+1e-8
    return count


def main():
    specs = [
        ((1, 2, 4, 9), 3, {1: 7, 5: 3}),
        ((1, 2, 4, 8, 17), 11, {0: 5, 7: 19}),
        ((1, 2, 4, 8, 16, 33), 9, {2: 17, 6: 41, 31: 5}),
        ((1, 2, 4, 8, 17, 0), 32, {0: 3, 9: 7, 31: 55}),
    ]
    total = literal = 0
    for columns, offset, edited in specs:
        data = setup(columns, offset, edited)
        total += exact_fixture(data)
        literal += literal_fixture(data)
    print(f'PASS: {total} exact three-family certificates on four affine fixtures; '
          'all baseline/mixed coefficients, arbitrary-weight corrections, pair-Gram '
          'inversion, class masses, exception constants, and copy capacities; '
          f'{literal} literal Gaussian signed products with actual source-phase bounds, '
          'mixed orientations, row units, common content, and deficient-rank zero fibres.')


if __name__ == '__main__':
    main()
