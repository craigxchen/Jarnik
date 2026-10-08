"""Small exact checks; no factorization of large Pell coordinates."""

from fractions import Fraction
from itertools import combinations
from math import gcd, isqrt, prod


def primes_up_to(D):
    return [p for p in range(2, D+1)
            if all(p % d for d in range(2, isqrt(p)+1))]


def parity(value, primes):
    answer = 1 << len(primes)
    for j, p in enumerate(primes):
        while value % p == 0:
            answer ^= 1 << j
            value //= p
    assert value == 1
    return answer


def audit(points):
    N = sum(v*v for v in points[0])
    identities = relations = dependencies = 0
    for anchor in points:
        g = gcd(*anchor)
        v = tuple(a//g for a in anchor)
        Q, Z = N//(g*g), 2*N//g
        assert Q == sum(a*a for a in v)
        selected = {}
        for point in points:
            if point == anchor:
                continue
            s = N//g-v[0]*point[0]-v[1]*point[1]
            t = v[0]*point[1]-v[1]*point[0]
            assert s > 0 and t*t == s*(Z-s)
            assert 2*g*s == sum((a-b)**2 for a, b in zip(point, anchor))
            selected[s] = t
            identities += 1
        D = max(selected, default=0)
        if D == 0:
            assert len(points) == 1
            continue
        primes = primes_up_to(D)
        r = len(primes)
        values = list(selected)
        assert len(points) <= 2*len(values)+1
        assert len(points) <= 2*r+3 or Z <= (8*(r+2)*D)**(r+3)
        basis = {}
        for index, s in enumerate(values):
            vector, support = parity(s, primes), 1 << index
            while vector:
                pivot = vector.bit_length()-1
                if pivot not in basis:
                    basis[pivot] = vector, support
                    break
                prior, labels = basis[pivot]
                vector, support = vector ^ prior, support ^ labels
            if vector == 0:
                chosen = [values[j] for j in range(len(values)) if support >> j & 1]
                assert 0 < len(chosen) <= r+2 and len(chosen) % 2 == 0
                product = prod(chosen)
                assert isqrt(product)**2 == product
                square = Fraction(prod(selected[s] for s in chosen), isqrt(product))
                P = prod(Z-s for s in chosen)
                assert square*square == P and square.denominator == 1
                dependencies += 1
                break
        for size in (2, 4):
            for chosen in combinations(values, size):
                product = prod(chosen)
                if isqrt(product)**2 != product:
                    continue
                square = Fraction(prod(selected[s] for s in chosen), isqrt(product))
                assert square*square == prod(Z-s for s in chosen)
                assert square.denominator == 1
                relations += 1
    return identities, relations, dependencies


def main():
    counts = [0, 0, 0]
    for N in (5, 25, 50, 100, 125, 200):
        points = []
        for x in range(-isqrt(N), isqrt(N)+1):
            y = isqrt(N-x*x)
            if x*x+y*y == N:
                points.append((x, y))
                if y:
                    points.append((x, -y))
        counts = [a+b for a, b in zip(counts, audit(points))]
    # Five actual integer-radius points, with only two positive deficits.
    s, t = 3, 2
    for index in range(1, 13):
        g = t*t+1
        points = [(g, 0), (g-1, s), (g-1, -s), (g-2, 2*t), (g-2, -2*t)]
        assert all(x*x+y*y == g*g for x, y in points)
        # At the first anchor Q=1, and the normalized deficits are 1,2.
        assert [g-x for x, y in points[1:]] == [1, 1, 2, 2]
        assert all(y*y == (g-x)*(2*g-(g-x)) for x, y in points)
        if 2*g > 48**4:
            assert len(points) == 2*len(primes_up_to(2))+3
        s, t = 3*s+4*t, 2*s+3*t
    for D in range(16, 129):
        r = len(primes_up_to(D))
        assert 8*(r+2)*D <= D**4
    assert counts[2] > 0
    print(f'PASS: {counts[0]} arbitrary-circle normalized identities; '
          f'{counts[1]} even square-product reconstructions; {counts[2]} '
          'short parity dependencies; 12 literal Pell configurations; '
          '113 explicit threshold-base inequalities.')


if __name__ == '__main__':
    main()
