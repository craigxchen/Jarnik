"""Exact checks of the actual-anchor reciprocal stretch obstruction."""

from collections import Counter
from fractions import Fraction
from math import gcd, lcm

from check_gaussian_reflection_replacement import gconj, gmul, gnorm, ggcd, gdivexact
from check_mobius_conductor_transfer import fixtures, least_norm
from check_mobius_source_aligned_compression import gaussian_valuation
from check_isometric_virtual_anchor import factor, gaussian_prime


def realize(rows):
    """Realize all phases, then remove the full common Gaussian gcd."""
    ratios = [tuple(Fraction(v,gnorm(h)) for v in gmul(h,h)) for h in rows]
    denominator = lcm(*(v.denominator for q in ratios for v in q))
    raw = [tuple(int(denominator*v) for v in q) for q in ratios]
    common = raw[0]
    for z in raw[1:]:
        common = ggcd(common,z)
    points = [gdivexact(z,common) for z in raw]
    N = gnorm(points[0])
    assert all(gnorm(z) == N for z in points)
    assert N == least_norm(rows)
    return points,N


def gram(a,b,d):
    return a*a+b*b+d*d, (a*a-b*b-d*d,2*a*b)


def content(W,z):
    C = gmul(W,z)
    return C,gcd(*C)


def clips(a,b,d,z0,N):
    U,V = (a+d,-b),(a-d,b)
    _,W = gram(a,b,d)
    result = {}
    for p in factor(gnorm(W)):
        if p % 4 != 1 or N % p:
            continue
        pi = gaussian_prime(p)
        oriented = []
        for prime in (pi,gconj(pi)):
            e0 = gaussian_valuation(z0,prime)
            e = 0
            n = N
            while n % p == 0:
                n //= p
                e += 1
            lo,hi = sorted((gaussian_valuation(V,prime)-gaussian_valuation(U,prime),
                            gaussian_valuation(gconj(U),prime)
                            -gaussian_valuation(gconj(V),prime)))
            width = max(0,min(hi,e-e0)-max(lo,-e0))
            oriented.append(width)
        assert oriented[0] == oriented[1]
        if oriented[0]:
            result[p] = oriented[0]
    return result


def main():
    count = identities = profiles = 0
    rows_sets = fixtures()
    # Retain the exact mixed-parity virtual Pell samples as another source.
    rows_sets.append([(1,0),(-26,7),(3,-1),(-17,-1)])
    for a in range(1,5):
        for d in range(1,5):
            for b in range(-5,6):
                if gcd(a,b,d) != 1 or (a == d and b == 0):
                    continue
                delta = a*d
                T,W = gram(a,b,d)
                Tprime,Wprime = gram(d,-b,a)
                assert gnorm(W) == gnorm(Wprime) > 0
                for rows in rows_sets:
                    source,N = realize(rows)
                    images = [(a*x+b*y,d*y) for x,y in rows]
                    target,Nprime = realize(images)
                    C,g = content(W,source[0])
                    Cprime,gprime = content(Wprime,target[0])
                    Delta = Fraction(g,2*delta*N)
                    Deltaprime = Fraction(gprime,2*delta*Nprime)
                    forward = clips(a,b,d,source[0],N)
                    inverse = clips(d,-b,a,target[0],Nprime)
                    assert forward == inverse
                    clipped = 1
                    for p,e in forward.items():
                        clipped *= p**e
                    assert g % clipped == gprime % clipped == 0
                    assert N % clipped == Nprime % clipped == 0
                    E,Eprime = N//clipped,Nprime//clipped
                    kernel = 4*delta**2*E*Eprime
                    values = []
                    divisors = []
                    for H,J,z,w in zip(rows,images,source,target):
                        x = Fraction(gnorm(J),delta*gnorm(H))
                        formula = Fraction(T,2*delta)+Fraction(C[0]*z[0]+C[1]*z[1],2*delta*N)
                        inverse_formula = (Fraction(Tprime,2*delta)
                                           +Fraction(Cprime[0]*w[0]+Cprime[1]*w[1],
                                                     2*delta*Nprime))
                        assert x == formula
                        assert x*inverse_formula == 1
                        assert ((x-Fraction(T,2*delta))/Delta).denominator == 1
                        assert ((1/x-Fraction(Tprime,2*delta))/Deltaprime).denominator == 1
                        n,nprime = 2*delta*E*x,2*delta*Eprime/x
                        assert n.denominator == nprime.denominator == 1
                        assert n > 0 and nprime > 0 and n*nprime == kernel
                        assert kernel % int(n) == 0
                        divisors.append(int(n))
                        values.append(x)
                        identities += 1
                    assert max(Counter(values).values()) <= 2
                    assert max(Counter(divisors).values()) <= 2
                    m = len(rows)
                    if kernel <= 10**8:
                        tau = 1
                        for e in factor(kernel).values():
                            tau *= e+1
                        assert m <= 2*tau
                    if m > 4:
                        assert (m-4)**2*Delta*Deltaprime <= 16
                        assert (m-4)**2*clipped**2 <= 64*delta**2*N*Nprime
                    profiles += len(forward)
                    count += 1
    print(f'PASS: {count} actual-frame matrix/tuple cases; {identities} physical '
          f'projection, reciprocal, and positive divisor-product identities; '
          f'{profiles} retained prime profiles.')


if __name__ == '__main__':
    main()
