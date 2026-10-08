"""Exact physical secant gauges, clipping, and Gaussian denominator costs."""

from fractions import Fraction as F
from itertools import combinations
from math import gcd, lcm

from check_gaussian_reflection_replacement import gnorm, gmul, gconj, ggcd, gdivexact
from check_mobius_conductor_transfer import least_norm
from check_mobius_source_aligned_compression import gaussian_valuation
from check_isometric_virtual_anchor import factor, gaussian_prime, physical
from check_mobius_reciprocal_stretch_grid import clips


def sub(a,b):
    return a[0]-b[0],a[1]-b[1]


def det(a,b):
    return a[0]*b[1]-a[1]*b[0]


def fraction(n,d):
    numerator = gmul(n,gconj(d))
    return tuple(F(x,gnorm(d)) for x in numerator)


def halfrows(points):
    N = gnorm(points[0])
    out = []
    for z in points:
        q = gmul(z,gconj(points[0]))
        h = (q[0]+N,q[1])
        if h == (0,0):
            h = (0,1)
        content = gcd(*h)
        out.append(tuple(x//content for x in h))
    assert out[0] == (1,0)
    assert least_norm(out) == N
    return out


def center(points,j,k,l):
    u = sub(points[j],points[0])
    v = sub(points[l],points[k])
    D = det(u,v)
    if D == 0:
        return None
    t = det(sub(points[k],points[0]),v)
    raw = tuple(D*x+t*y for x,y in zip(points[0],u))
    r0 = gcd(*raw,D)
    sign = 1 if D > 0 else -1
    A = tuple(sign*x//r0 for x in raw)
    return A,abs(D)//r0,D,t,r0,u


def main():
    datasets = [
        [(5,0),(3,4),(-3,4),(-5,0),(-3,-4),(3,-4)],
        [(1,8),(4,7),(7,4),(8,1),(8,-1),(7,-4)],
        [(1,8),(4,7),(-7,4),(-8,-1),(8,-1),(7,-4)],
    ]
    cases = partners = prime_checks = seed_loss_checks = 0
    for points in datasets:
        N = gnorm(points[0])
        rows = halfrows(points)
        common = points[0]
        for z in points[1:]:
            common = ggcd(common,z)
        assert gnorm(common) == 1
        for j in range(1,len(points)):
            others = [k for k in range(1,len(points)) if k != j]
            for k,l in combinations(others,2):
                data = center(points,j,k,l)
                if data is None:
                    continue
                A,s,D,t,r0,u = data
                if A == (0,0):
                    continue
                K = gnorm(A)-N*s*s
                assert K != 0
                h0 = sub(tuple(s*x for x in points[0]),A)
                D0 = gnorm(h0)
                assert D0*r0*r0 == t*t*gnorm(u)
                assert K*r0*r0 == t*(t-D)*gnorm(u)
                beta = gmul(A,gconj(points[0]))[1]
                content = gcd(D0,2*s*beta,K)
                a,b,d = D0//content,-2*s*beta//content,-K//content
                assert abs(a*d)*content**2 == abs(K)*D0
                images = [(a*x+b*y,d*y) for x,y in rows]
                actual = []
                reduced = []
                hs = []
                L = (1,0)
                for z,image in zip(points,images):
                    h = sub(tuple(s*x for x in z),A)
                    hs.append(h)
                    numerator = sub(gmul(A,gconj(z)),(N*s,0))
                    denominator = gconj(h)
                    w = fraction(numerator,denominator)
                    assert sum(x*x for x in w) == N
                    assert physical(points[j],image) == w
                    actual.append(w)
                    g = gcd(*h)
                    v = tuple(x//g for x in h)
                    S = gnorm(v)
                    criterion = K % (g*S) == 0
                    if criterion:
                        q = K//(g*S)
                        criterion = all((x+q*y) % s == 0 for x,y in zip(A,v))
                    assert criterion == all(x.denominator == 1 for x in w)
                    T = s*g*S
                    n = tuple(x*g*S+K*y for x,y in zip(A,v))
                    rational_denominator = T//gcd(T,*n)
                    assert rational_denominator == lcm(*(x.denominator for x in w))
                    common = ggcd(numerator,denominator)
                    nr,br = gdivexact(numerator,common),gdivexact(denominator,common)
                    assert gnorm(ggcd(nr,br)) == 1
                    reduced.append((nr,br))
                    L = gdivexact(gmul(L,br),ggcd(L,br))
                    partners += 1
                assert actual[0] == points[j] and actual[j] == points[0]
                assert actual[k] == points[l] and actual[l] == points[k]
                integral_targets = [gdivexact(gmul(L,n),b) for n,b in reduced]
                G = integral_targets[0]
                for w in integral_targets[1:]:
                    G = ggcd(G,w)
                assert gnorm(ggcd(G,L)) == 1
                assert gnorm(gdivexact(points[j],G))*gnorm(G) == N
                Nprime = N*gnorm(L)//gnorm(G)
                assert Nprime == least_norm(images)
                assert Nprime % gnorm(L) == 0
                union_common = G
                for z in points:
                    union_common = ggcd(union_common,gmul(L,z))
                assert gnorm(union_common) == 1
                Q = 1
                clipped_exponents = {}
                for p,e in factor(N).items():
                    pi = gaussian_prime(p)
                    widths = []
                    for prime in (pi,gconj(pi)):
                        r = gaussian_valuation(A,prime)
                        rp = gaussian_valuation(A,gconj(prime))
                        power_s = 0
                        ss = s
                        while ss % p == 0:
                            ss //= p
                            power_s += 1
                        lo,hi = sorted((r-power_s,e+power_s-rp))
                        widths.append(max(0,min(e,hi)-max(0,lo)))
                    assert widths[0] == widths[1]
                    clipped_exponents[p] = widths[0]
                    Q *= p**widths[0]
                direct = 1
                for p,e in clips(a,b,d,points[0],N).items():
                    direct *= p**e
                assert direct == Q
                Ared = gdivexact(A,ggcd(A,(s,0)))
                assert gcd(N,gnorm(Ared)) % (N//Q) == 0
                for p,e in factor(N).items():
                    pi = gaussian_prime(p)
                    allocations = [gaussian_valuation(z,pi) for z in points]
                    if all(v in (0,e) for v in allocations):
                        seed_ones = sum(allocations[index] == e for index in (0,j,k,l))
                        if seed_ones in (1,3):
                            assert clipped_exponents[p] == 0
                            seed_loss_checks += 1
                # Exhaust small denominator norms, plus a fixed range of larger-case probes.
                primes = {5,13,17,29,37,41,53,61,73,89,97,101,109,113}
                if gnorm(L) <= 10**9:
                    primes.update(factor(gnorm(L)))
                for p in primes:
                    if p % 4 != 1 or N*s*K % p == 0:
                        continue
                    pi = gaussian_prime(p)
                    expected = max(gaussian_valuation(h,pi) for h in hs)
                    expected += max(gaussian_valuation(h,gconj(pi)) for h in hs)
                    n = Nprime
                    actual_exponent = 0
                    while n % p == 0:
                        n //= p
                        actual_exponent += 1
                    assert actual_exponent == expected
                    prime_checks += 1
                cases += 1
    print(f'PASS: {cases} actual four-seed secant maps; {partners} partner and '
          f'integrality identities; {prime_checks} oriented new-prime checks; '
          f'{seed_loss_checks} three-against-one seed erasures; '
          'exact union/target norms and source clipping in every case.')


if __name__ == '__main__':
    main()
