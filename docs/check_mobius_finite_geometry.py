"""Exact finite geometry audits for lattice-circle projective maps."""
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, product
from math import gcd, isqrt, lcm

from check_gaussian_reflection_replacement import gnorm
from check_mobius_conformal_normalization import image, invariants
from check_mobius_conductor_transfer import fixtures, least_norm


def det(x,y):
    return x[0]*y[1]-x[1]*y[0]


def sine2(x,y):
    return Fraction(det(x,y)**2,gnorm(x)*gnorm(y))


def exponent(q):
    return Fraction(q*q//4,q*(q-1))


def diameter_checks(rows,N):
    count = 0
    for q in range(2,len(rows)+1):
        for subset in combinations(rows,q):
            d2 = max(sine2(x,y) for x,y in combinations(subset,2))
            assert d2**(q*(q-1))*(2*N)**(2*(q*q//4)) >= 1
            count += 1
    return count


def main():
    l1_cases = 0
    for m in range(3,10):
        F = (m-1)**2//4
        for xs in combinations_with_replacement(range(-3,4),m):
            assert m*sum(abs(x+y) for x,y in combinations(xs,2)) >= 2*F*sum(map(abs,xs))
            l1_cases += 1
    sources = fixtures()
    sources.append([(1,0),(1,1),(1,2),(3,1),(3,2),(5,1),(5,2),(2,5)])
    source_norms = [least_norm(zs) for zs in sources]
    diameters = sum(diameter_checks(zs,N) for zs,N in zip(sources,source_norms))
    matrices = []
    for M in product(range(-2,3),repeat=4):
        if M[0]*M[3]-M[1]*M[2] and gcd(*M) == 1:
            matrices.append(M)
    matrices += [(1,10**k,0,1) for k in range(1,7)]
    cases = identities = 0
    for M in matrices:
        a,b,c,d = M
        delta = abs(a*d-b*c)
        T,K = invariants(M)
        # A rational upper enclosure of the exact condition number.
        root_ceiling = isqrt(K)
        if root_ceiling*root_ceiling < K:
            root_ceiling += 1
        condition_upper = Fraction(T+root_ceiling,2*delta)
        for rows,N in zip(sources,source_norms):
            targets = image(rows,M)
            Nprime = least_norm(targets)
            raw = [(a*x+b*y,c*x+d*y) for x,y in rows]
            stretches2 = [Fraction(gnorm(z),delta*gnorm(h)) for z,h in zip(raw,rows)]
            m = len(rows)
            for i,j in combinations(range(m),2):
                assert sine2(targets[i],targets[j])*stretches2[i]*stretches2[j] == sine2(rows[i],rows[j])
                identities += 1
            for n in range(2,m):
                r = m-n+1
                hn,hr = exponent(n),exponent(r)
                power = lcm(hn.denominator,hr.denominator)
                assert (condition_upper/4)**power <= (2*N)**int(hn*power)*(2*Nprime)**int(hr*power)
            cases += 1
    print(f'PASS: {l1_cases} exact signed L1 profiles and {diameters} inherited-radius subset diameter bounds.')
    print(f'PASS: {cases} matrix/configuration cases and {identities} exact pair stretch identities.')
    print('PASS: all shared-pivot condition-number upper bounds, including mixed units and shears through 10^6.')


if __name__ == '__main__':
    main()
