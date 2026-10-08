"""Exact rational audit; no large-cardinality or map-existence claim."""

from fractions import Fraction
from math import gcd, isqrt

from check_gaussian_reflection_replacement import gconj, gmul, gnorm
from check_mobius_conductor_transfer import least_norm
from check_mobius_source_aligned_compression import gaussian_valuation


def factor(n):
    out = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def gaussian_prime(p):
    assert p % 4 == 1
    for y in range(1, isqrt(p) + 1):
        x = isqrt(p - y*y)
        if x*x + y*y == p:
            return x, y
    raise AssertionError(p)


def physical(reference, halfrow):
    numerator = gmul(reference, gmul(halfrow, halfrow))
    return tuple(Fraction(v, gnorm(halfrow)) for v in numerator)


def matmul(A, B):
    return tuple(tuple(sum(A[i][k] * B[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))


def normalize_audit():
    count = 0
    for a in range(1, 9):
        for d in range(1, 9):
            if a == d:
                continue
            for b in range(1, 151):
                delta = a*d
                D = b*b + (a-d)**2
                s = isqrt(delta*D)
                if s*s != delta*D:
                    continue
                norms = []
                for sign in (-1, 1):
                    x, y = -a*b + sign*s, a*(a-d)
                    n = x*x+y*y
                    norms.append(n)
                    u, v = a*x+b*y, d*y
                    assert u*u+v*v == delta*n
                    result = matmul(((u,v),(-v,u)),
                                    matmul(((a,b),(0,d)), ((x,-y),(y,x))))
                    assert result[0][0] == result[1][1] == delta*n
                    assert result[1][0] == 0
                    assert abs(result[0][1]) == s*n
                    count += 1
                assert norms[0]*norms[1] == a*a*(a-d)**2*D*(D+4*delta)
    return count


def pell_audit():
    b, c = 1, 1
    cases = local = 0
    for _ in range(6):
        assert b*b - 2*c*c == -1
        N = c*c*(c*c+4)
        z0, z1 = (-(c*c+1), -b), (-(c*c+1), b)
        w0, w1 = (c*c-2, -2*b), (c*c-2, 2*b)
        assert all(gnorm(z) == N and gcd(*z) == 1 for z in (z0,z1,w0,w1))
        source = [(1,0),z1]
        target = [(1,0),w1]
        assert least_norm(source) == least_norm(target) == N
        assert physical(z0,z1) == z1
        assert physical(w0,w1) == w1
        Hs = [(-b+2*c,-1),(-b-2*c,-1)]
        for sign,H in zip((1,-1),Hs):
            raw = (H[0]+b*H[1],2*H[1])
            assert gnorm(raw) == 2*gnorm(H)
            vp = physical(z0,H)
            wp = physical(w0,raw)
            assert vp == (-c*c,sign*2*c)
            assert wp == (c*c,sign*2*c)
            assert gcd(*(int(x) for x in vp)) == c
            assert gcd(*(int(x) for x in wp)) == c
            source.append(H)
            target.append(raw)
            assert least_norm(source) == least_norm(target) == N
        U,V = (3,-b),(-1,b)
        assert gmul(U,gconj(V)) == tuple(2*x for x in z0)
        primes = factor(c)
        for p,e in factor(c*c+4).items():
            assert p not in primes
            primes[p] = e
        Nclip = 1
        for p in primes:
            pi = gaussian_prime(p)
            for orient in (pi,gconj(pi)):
                e0 = gaussian_valuation(z0,orient)
                e1 = gaussian_valuation(z1,orient)
                J = sorted((-e0, e1))
                I = sorted((gaussian_valuation(V,orient)-gaussian_valuation(U,orient),
                            gaussian_valuation(gconj(U),orient)
                            -gaussian_valuation(gconj(V),orient)))
                assert J == I
                e = J[1]-J[0]
                assert e > 0
                for H in Hs:
                    r = (gaussian_valuation(H,orient)
                         -gaussian_valuation(gconj(H),orient))
                    assert J[0] <= r <= J[1]
                    depth = min(r-J[0],J[1]-r)
                    assert depth == factor(c).get(p,0)
                    if c % p == 0:
                        assert 2*r == sum(J)
                local += 1
            Nclip *= p**e
        assert Nclip == N == c*(c*c+4)*c
        if c >= 5:
            # Certificates used with arctan(t)<=t for the C<=8 arc bounds.
            assert 625*N <= (6*c)**4
            assert 2*b <= 3*c
            assert 10*(c*c-2) >= 9*c*c
        cases += 1
        b,c = 3*b+4*c,2*b+3*c
    return cases,local


def virtual_width_audit():
    cases = 0
    for L in range(-5,1):
        for U in range(0,6):
            for r in range(-8,9):
                outside = max(L-r,0,r-U)
                depth = max(0,min(r-L,U-r))
                assert max(U,r)-min(L,r) == U-L+outside
                assert min(r-min(L,r),max(U,r)-r) == depth
                assert outside*depth == 0
                assert depth <= min(-L,U)+abs(r)
                for t in range(1,6):
                    clipped = max(0,min(U-r,t)-max(L-r,-t))
                    assert clipped <= t+depth
                cases += 1
    return cases


if __name__ == '__main__':
    count = normalize_audit()
    pell,local = pell_audit()
    widths = virtual_width_audit()
    print(f'PASS: {count} exact isometric normalizations; {pell} Pell tuples; '
          f'{local} oriented prime profiles; {widths} virtual width profiles.')
