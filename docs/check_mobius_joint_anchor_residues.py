"""Exact finite audits of anchor intervals and joint endpoint residues.

These checks verify algebraic identities, not the unproved uniform bound.
"""
from functools import lru_cache
from fractions import Fraction
from itertools import combinations, product
from math import gcd, isqrt, prod

from check_gaussian_reflection_replacement import gmul, gconj, gnorm, ggcd, gdivexact
from check_mobius_conductor_transfer import fixtures, least_norm, residues
from check_mobius_conformal_normalization import invariants, image
from check_mobius_actual_anchor import reduced_coefficients, minimal_matrix


@lru_cache(None)
def factors(n):
    out = []
    p = 2
    while p*p <= n:
        if n % p == 0:
            out.append(p)
            while n % p == 0:
                n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        out.append(n)
    return tuple(out)


@lru_cache(None)
def gaussian_prime(p):
    assert p % 4 == 1
    for x in range(1, isqrt(p)+1):
        y = isqrt(p-x*x)
        if x*x+y*y == p:
            return x, y
    raise AssertionError(p)


def gv(z, pi):
    n = gnorm(pi)
    out = 0
    assert z != (0, 0)
    while True:
        a, b = gmul(z, gconj(pi))
        if a % n or b % n:
            return out
        z = (a//n, b//n)
        out += 1


def gpow(z, k):
    out = (1, 0)
    for _ in range(k):
        out = gmul(out, z)
    return out


def determinant(M):
    a,b,c,d = M
    return a*d-b*c


def target_pair_bound(rows):
    norms = [gnorm(z) for z in rows]
    pair_norms = [gnorm(x)*gnorm(y)//gnorm(ggcd(x,y))**2
                  for x,y in combinations(rows, 2)]
    N = least_norm(rows)
    m, r = len(rows), sum(n % 2 == 0 for n in norms)
    assert prod(pair_norms)**4 <= 2**(4*r*(m-r))*N**(m*m)
    return N


def audit(rows, M):
    u,v = reduced_coefficients(M)
    M, eps = minimal_matrix(u,v)
    T,K = invariants(M)
    U,V = gnorm(u), gnorm(v)
    k0 = gcd(U,V)
    assert determinant(M) % k0 == 0
    core = abs(determinant(M))//k0
    N = least_norm(rows)
    targets = image(rows,M)
    Nprime = target_pair_bound(targets)
    old_b, new_b = residues(rows), residues(targets)
    m = len(rows)
    k = [1]*m
    clipped_product = 1
    clipped_norm = 1
    external = 1
    ps = set(factors(U*V))
    for h in rows:
        ps.update(factors(gnorm(h)))
    for p in sorted(ps):
        if p % 4 != 1:
            continue
        pi = gaussian_prime(p)
        barpi = gconj(pi)
        d = gv(v,pi)-gv(u,pi)
        e = gv(v,barpi)-gv(u,barpi)
        left,right = sorted((d,-e))
        ts = [gv(h,pi)-gv(h,barpi) for h in rows]
        lo,hi = min(ts),max(ts)
        external *= p**max(0,left-hi,lo-right)
        clips = [max(left,min(right,t)) for t in ts]
        output_ts = [gv(h,pi)-gv(h,barpi) for h in targets]
        assert all(abs(x-y) >= abs(c-d)
                   for (x,y),(c,d) in zip(combinations(output_ts,2),combinations(clips,2)))
        clipped_norm *= p**(max(clips)-min(clips))
        for j,(t,clipped) in enumerate(zip(ts,clips)):
            k[j] *= p**abs(t-clipped)
        clipped_product *= p**sum(abs(x-y) for x,y in combinations(clips,2))
    assert k[0] == k0
    assert Nprime % clipped_norm == 0
    assert gcd(external,Nprime) == 1
    assert all(b % external == 0 for b in new_b)
    assert all(kj % external == 0 and external*N % kj == 0 for kj in k)
    assert N**(m-1) % prod(kj//external for kj in k) == 0
    pair_product = prod(gnorm(x)*gnorm(y)//gnorm(ggcd(x,y))**2
                        for x,y in combinations(rows,2))
    source_r = sum(gnorm(z) % 2 == 0 for z in rows)
    assert pair_product % 2**(source_r*(m-source_r)) == 0
    odd_pair_product = pair_product//2**(source_r*(m-source_r))
    overlap = prod(k)**(m-1)*clipped_product
    assert prod(new_b)**2*pair_product % overlap == 0
    if m >= 4:
        cut_k = m//4
        cut_d = m-2*cut_k-1
        cut_c = cut_k*((m-1)*cut_d+1)
        lhs = N**cut_c*odd_pair_product**(4*cut_k)
        rhs = clipped_norm**cut_c*prod(new_b)**(4*cut_d)*N**(cut_k*m*m)
        assert rhs % lhs == 0
    for h,kj in zip(rows,k):
        uh,vh = gmul(u,h),gmul(v,gconj(h))
        common = ggcd(uh,vh)
        assert gcd(gnorm(gdivexact(uh,common)),gnorm(gdivexact(vh,common))) == kj
    odd_core = core
    while odd_core % 2 == 0:
        odd_core //= 2
    F = (m-1)**2//4
    assert prod(old_b)*prod(new_b) % odd_core**F == 0
    r = sum(gnorm(z) % 2 == 0 for z in targets)
    parity_charge = r*(r-1)//2+source_r*(source_r-1)//2
    assert (2**parity_charge*prod(old_b)*prod(new_b)) % core**F == 0
    J = gcd(K,abs(determinant(M))*N)
    for lam in ((1,1),(2,1),(1,2)):
        x,y = lam
        a,b,c,d = M
        padded = (x*a-y*c,x*b-y*d,y*a+x*c,y*b+x*d)
        s = gnorm(lam)
        tp,kp = invariants(padded)
        jp = gcd(kp,abs(determinant(padded))*N)
        assert jp == s*gcd(s*K,abs(determinant(M))*N)
        assert jp % (s*J) == 0
        assert tp == s*T
    return int(external > 1), int(odd_core > 1)


def main():
    for m in range(4,101):
        k = m//4
        d = m-2*k-1
        A = Fraction(k,d)
        c = Fraction(k*((m-1)*d+1),4*d)
        for n in range(1,m):
            assert Fraction(n*(n-1),2)+A*(n-Fraction(m,2))**2 >= c
        retained = (c-A*Fraction(m,4))/(c+Fraction(m,8))
        assert retained == Fraction(2*k*(m-2*k-2),m+2*k*(m-2*k-2))
        assert retained >= 1-Fraction(4,m)
        erased = (Fraction(m,8)+A*Fraction(m,4))/c
        F = (m-1)**2//4
        singularity = (1-erased)/2-Fraction(m,4*F)
        if m >= 8:
            assert singularity >= Fraction(1,66)
            assert c >= Fraction(m*(m-4),16)
            full_parity_H = Fraction(m*m,16*c)*(1+2*A)+Fraction(m*(m-1),F)
            assert full_parity_H < 8
        assert Fraction(m*m,8)*(1+2*A)/(c+Fraction(m,8)) <= 4
        if m == 8:
            assert singularity == Fraction(1,66)
    cases = extcases = corecases = 0
    source_sets = fixtures()
    source_sets.append([(1,0),(1,1),(1,2),(3,1),(3,2),(5,1),(5,2),(2,5)])
    for entries in product(range(-2,3),repeat=4):
        if not determinant(entries) or gcd(*entries) != 1 or invariants(entries)[1] == 0:
            continue
        for rows in source_sets:
            e,c = audit(rows,entries)
            cases += 1
            extcases += e
            corecases += c
    for ell in range(1,5):
        pi = (4,1)
        u,v = gmul((2,0),gpow(gconj(pi),ell)),gpow(pi,ell)
        M,_ = minimal_matrix(u,v)
        for rows in source_sets:
            e,c = audit(rows,M)
            assert e == 1
            cases += 1
            extcases += e
            corecases += c
    for p in (2,3,5,7,17):
        for ell in range(1,4):
            q = p**ell
            M = (1,q,1,0)
            for rows in source_sets + [[(1,0),(q,2),(2*q,1)]]:
                # For p=2 replace the nonprimitive targeted second row.
                if any(gcd(*z) != 1 for z in rows):
                    continue
                e,c = audit(rows,M)
                cases += 1
                extcases += e
                corecases += c
    print(f'PASS: {cases} joint matrix/configuration cases; {extcases} with external conductor, {corecases} with nontrivial odd determinant core.')
    print('PASS: exact anchor distances, external target residues, internal overlap divisibility, and full core packing including two.')
    print('PASS: every clipped pair valuation and the clipped source conductor survive in the target.')
    print('PASS: exact removed-cut balance inequality, integer conductor charges, and optimized retention exponents for m=4,...,100.')
    print('PASS: resulting condition-number exponents are positive for every checked m>=8, with z_8=1/66.')
    print('PASS: arbitrary source and target parity, their full core charge, and the universal retention and condition-number constants.')
    print('PASS: monotonicity of source-capped content under conformal padding.')


if __name__ == '__main__':
    main()
