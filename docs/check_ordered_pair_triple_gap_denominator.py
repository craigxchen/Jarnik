"""Exact gap/overlap and actual subgroup-gcd checks; no endpoint inference."""

from functools import reduce
from itertools import product
from math import gcd, prod
from random import Random

from check_ordered_affine_circuit_conductor_index import (
    band_index, gconj, ggcd, gmul, gnorm, gpower, gdiv_exact, minors,
)


def gap_overlap(t):
    lo, hi = sorted((t[0], t[4]))
    li, hii = min(t[1:4]), max(t[1:4])
    return max(0, lo-hii, li-hi), max(0, min(hi,hii)-max(lo,li))


def check_tuple(points, expected=None):
    common = reduce(ggcd, points)
    points = [gdiv_exact(z, common) for z in points]
    n = gnorm(points[0])
    assert all(gnorm(z) == n for z in points)
    ao = gnorm(ggcd(points[0], points[4]))
    ai = gnorm(reduce(ggcd, points[1:4]))
    no, ni = n//ao, n//ai
    k, ell = n//gcd(n,no*ni), no*ni//gcd(n,no*ni)
    assert k == ao*ai//gcd(n,ao*ai)
    assert ell == n//gcd(n,ao*ai)
    assert gcd(k,ell) == 1 and n%(k*ell) == 0
    assert no*ni*k == n*ell
    assert ai%k == 0 and ao%k == 0
    if expected is not None:
        assert (k,ell) == expected
    if len(set(points)) == 5:
        index,_,_ = band_index(points)
        assert gcd(index,n)%k == 0
        d = abs(minors(points)[(1,2,3)])
        assert d%(2*ai) == 0 and 2*k <= d
        assert n%2 == 1
        assert all((x+y)%2 for x,y in points)
    return k,ell


def main():
    allocations = 0
    for e in range(1,7):
        for t in product(range(e+1), repeat=5):
            if min(t) != 0 or max(t) != e:
                continue
            h,l = gap_overlap(t)
            ro = abs(t[0]-t[4])
            ri = max(t[1:4])-min(t[1:4])
            assert ro+ri-e == l-h
            assert h == max(0,e-ro-ri)
            count = sum({i for i in range(5) if t[i]>=a}
                        in ({0,4},{1,2,3}) for a in range(1,e+1))
            assert h == count
            allocations += 1
    rng = Random(921517)
    pis = ((2,1),(3,2),(4,1))
    units = ((1,0),(0,1),(-1,0),(0,-1))
    actual = 0
    gap_cases = overlap_cases = 0
    for _ in range(600):
        points = [rng.choice(units) for _ in range(5)]
        k = ell = 1
        for pi in pis:
            e = rng.randrange(1,5)
            t = [rng.randrange(e+1) for _ in range(5)]
            h,l = gap_overlap(t)
            k *= gnorm(pi)**h
            ell *= gnorm(pi)**l
            for j in range(5):
                points[j] = gmul(points[j], gmul(gpower(pi,t[j]),
                                               gpower(gconj(pi),e-t[j])))
        # Common-layer removal preserves both interval gap and overlap.
        got = check_tuple(points,(k,ell))
        gap_cases += got[0]>1
        overlap_cases += got[1]>1
        actual += 1
    for points in (
        [(4,-33),(9,-32),(12,-31),(23,-24),(24,-23)],
        [(-110,-35),(-109,-38),(-98,-61),(-94,-67),(-86,-77)],
    ):
        check_tuple(points)
        actual += 1
    assert gap_cases and overlap_cases
    print(f"Checked {allocations} allocations and {actual} actual tuples; "
          f"{gap_cases} gap and {overlap_cases} overlap cases.")


if __name__ == '__main__':
    main()
