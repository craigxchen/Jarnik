"""Exact algebraic checks of the four-row angular squareclass proof.

The analytic inequality and nonzero-word assertion are proved in the note;
these checks cover the allocation inequality and every Gaussian root identity.
"""
from itertools import combinations, product
from math import prod
from check_gaussian_reflection_replacement import gmul, gconj, gnorm


def power(z,n):
    out = (1,0)
    for _ in range(n):
        out = gmul(out,z)
    return out


def halfsums(es):
    a,b,c,d = es
    assert (a+b+c+d) % 2 == 0
    return ((a+b-c-d)//2,(a+c-b-d)//2,(a+d-b-c)//2)


def main():
    valuation_cases = 0
    for es in product(range(10),repeat=4):
        if sum(es) % 2:
            continue
        ss = halfsums(es)
        assert 2*sum(map(abs,ss)) <= 3*(max(es)-min(es))
        valuation_cases += 1
    blocks = [(2,1),(3,2)]
    rows = list(product(range(3),repeat=2))
    N = prod(gnorm(z)**2 for z in blocks)
    partitions = [(0,1,2,3),(0,2,1,3),(0,3,1,2)]
    root_cases = 0
    for chosen in combinations(rows,4):
        if any(sum(row[j] for row in chosen) % 2 for j in range(2)):
            continue
        source = []
        for row in chosen:
            z = (1,0)
            for h,e in zip(blocks,row):
                z = gmul(z,gmul(power(h,e),power(gconj(h),2-e)))
            assert gnorm(z) == N
            source.append(z)
        norms = []
        for j,(a,b,c,d) in enumerate(partitions):
            u = (1,0)
            for k,h in enumerate(blocks):
                s = halfsums([row[k] for row in chosen])[j]
                u = gmul(u,power(h if s >= 0 else gconj(h),abs(s)))
            norms.append(gnorm(u))
            for unit in [(1,0),(0,1),(-1,0),(0,-1)]:
                w = gmul(unit,u)
                # (w/conj(w))^2 = z_a z_b/(z_c z_d), with no denominators.
                assert gmul(gmul(w,w),gmul(source[c],source[d])) == gmul(gmul(gconj(w),gconj(w)),gmul(source[a],source[b]))
                root_cases += 1
        assert prod(norms)**2 <= N**3
    print(f'PASS: {valuation_cases} integer four-row half-sum inequalities with entries zero through nine.')
    print(f'PASS: {root_cases} exact Gaussian root identities, including all four unit choices.')


if __name__ == '__main__':
    main()
