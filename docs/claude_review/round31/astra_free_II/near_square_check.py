"""Exact near-square separation tests. Requires only the standard library."""
from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import gcd, isqrt, prod


def mul(z, w):
    return z[0]*w[0] - z[1]*w[1], z[0]*w[1] + z[1]*w[0]


def conj(z):
    return z[0], -z[1]


def norm(z):
    return z[0]**2 + z[1]**2


def quotient(z, w):
    n = norm(w)
    x, y = mul(z, conj(w))
    assert x % n == y % n == 0
    return x//n, y//n


def gaussian_gcd(z, w):
    while w != (0, 0):
        n = norm(w)
        x, y = mul(z, conj(w))
        # Round each coordinate of the quotient to a nearest integer.
        q = ((2*x+n)//(2*n), (2*y+n)//(2*n))
        a, b = mul(q, w)
        z, w = w, (z[0]-a, z[1]-b)
    return z


def primary(z):
    # z - 1 is divisible by 2+2i, an associate of (1+i)^3.
    x, y = z
    return (x-1+y) % 4 == (y-x+1) % 4 == 0


def check_cluster():
    n = 10**30+1
    x, y, z, u = -3, -2, 6, -6
    blocks = []
    for a, b, c in [(x, y, z), (x, z, y), (y, z, x)]:
        assert (u*(1+c)) % c == 0
        blocks.append({(0, 1), (u, a-1), (u, b-1), (u*(1+c)//c, c)})
    roots = sorted(set.union(*blocks))
    points = []
    for chosen in [set(), *blocks]:
        factors = [(n-4*re, (-4 if (re, im) in chosen else 4)*im) for re, im in roots]
        points.append(reduce(mul, factors, (1, 0)))
    common = reduce(gaussian_gcd, points)
    common = next(w for unit in [(1, 0), (0, 1), (-1, 0), (0, -1)]
                  if primary(w := mul(common, unit)))
    points = [quotient(v, common) for v in points]
    N = norm(points[0])
    assert len(set(points)) == 4 and all(primary(v) and norm(v) == N for v in points)
    assert norm(reduce(gaussian_gcd, points)) == 1

    # tan((theta_i-theta_0)/2). Choose theta_i-theta_0 in (-pi,pi).
    tangents = []
    for v in points:
        re, im = mul(v, conj(points[0]))
        assert N+re > 0
        tangents.append(Fraction(im, N+re))
    spread = max(tangents)-min(tangents)
    # atan is 1-Lipschitz. Arc length <= 2 sqrt(N) spread < N^(1/4)/2.
    assert (4*spread)**4*N < 1

    divisors, offsets = {}, []
    for i, j in combinations(range(len(points)), 2):
        v, w = points[i], points[j]
        re, im = mul(v, conj(w))
        assert N+re > 0
        common = gcd(N+re, im)
        s, t = (N+re)//common, im//common
        d = s*s+t*t
        assert gcd(s, t) == 1 and d % 2 and N % d == 0
        assert d == N//norm(gaussian_gcd(v, w))
        assert t and t*t < 2*s+1 and d-isqrt(d)**2 == t*t
        divisors[i, j] = d
        offsets.append(t*t)
    assert len(set(divisors.values())) == 6
    for i, j, k in combinations(range(len(points)), 3):
        triple = divisors[i, j] * divisors[i, k] * divisors[j, k]
        assert isqrt(triple)**2 == triple
    print('Actual cluster: norm digits', len(str(N)), 'distinct separation divisors', len(divisors))
    print('Square offsets in lexicographic edge order:', offsets)
    print('All four triangular products are squares.')
    for power in range(31):
        if (4*spread)**4*N >= Fraction(1, 10**power):
            break
    print(f'Certified arc length < N^(1/4)/2; (4*spread)^4*N < 10^(-{power-1})')


def check_imprimitive():
    y = 200
    primes = [p for p in range(y//2, y+1) if p % 4 == 1 and
              all(p % q for q in range(2, isqrt(p)+1))]
    L = prod(primes)
    N = L*L*((32*L)**2+1)
    seen = set()
    for p, q in combinations(primes, 2):
        k = p*q
        d = k*k*((32*L)**2+1)
        s = 32*k*L
        assert N % d == 0 and d not in seen and isqrt(d) == s
        assert d-s*s == k*k and gcd(s, k) == k and k > 1
        seen.add(d)
    print('Imprimitive stress test: primes', primes, 'distinct divisors', len(seen))
    print('Norm digits', len(str(N)), 'largest square offset', max(d-isqrt(d)**2 for d in seen))


if __name__ == '__main__':
    check_cluster()
    check_imprimitive()
