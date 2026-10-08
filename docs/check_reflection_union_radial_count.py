"""Exact reflection-union and radial-spacing checks on actual Pell arcs."""
from fractions import Fraction
from itertools import combinations

from check_gaussian_near_real_divisor_reduction import tuple_from_columns
from check_gaussian_reflection_replacement import gmul, gconj, gnorm, ggcd, gcd_many, gdivexact


def main():
    count = 0
    for index in (1,11,21):
        u,v = 1,0
        for _ in range(index):
            u,v = 9*u+20*v,4*u+9*v
        blocks = [(-1,2),(u+v,2*v),(2*v,u-v),(2*u+8*v,3*u+v)]
        allocations = [(1,1,1,1),(0,1,0,0),(0,0,1,0),(0,0,0,1)]
        rows = tuple_from_columns(blocks,[1]*4,allocations)
        N = gnorm(rows[0])
        assert all(gnorm(z) == N for z in rows)
        assert gnorm(gcd_many(rows)) == 1
        first,last = rows[1],rows[2]
        x,y = gmul(last,gconj(first))
        tangent = Fraction(y,N+x)
        assert x > 0 and y > 0 and tangent**4*N < 1
        for i,j in combinations(range(4),2):
            numerator = gmul(rows[i],rows[j])
            common = ggcd(numerator,(N,0))
            A,D = gdivexact(numerator,common),gdivexact((N,0),common)
            Q = gnorm(D)
            union = set(gmul(D,z) for z in rows) | set(gmul(A,gconj(z)) for z in rows)
            assert len(union) == 6 and gnorm(gcd_many(union)) == 1
            assert all(gnorm(z) == N*Q for z in union)
            reflection = gdivexact(A,gconj(D))
            assert gnorm(reflection) == 1 and all(gmul(reflection,gconj(z)) in union for z in union)
            if reflection == (1,0):
                levels = sorted(set(x for x,y in union))
                spacing = 1
            elif reflection == (-1,0):
                levels = sorted(set(y for x,y in union))
                spacing = 1
            else:
                sign = reflection[1]
                levels = sorted(set(x+sign*y for x,y in union))
                assert all(n % 2 for n in levels)
                spacing = 2
            assert len(levels) == 3 and (levels[0] > 0 or levels[-1] < 0)
            assert min(b-a for a,b in zip(levels,levels[1:])) >= spacing
            assert levels[-1]-levels[0] >= spacing*(len(levels)-1)
            # C_upper=2*tan(Delta/2)*N^(1/4) bounds the actual C.
            assert Q*16*tangent**4*N >= 4*(4-2)**2
            if (i,j) == (1,2):
                assert Q*16*tangent**4*N >= 64*(4-2)**2
            count += 1
    print(f'PASS: {count} actual Pell reflection unions, each primitive with six points and a Gaussian-unit reflection.')
    print('PASS: exact coordinate/diagonal radial spacing, full parity accounting, and both radius-cost consequences.')


if __name__ == '__main__':
    main()
