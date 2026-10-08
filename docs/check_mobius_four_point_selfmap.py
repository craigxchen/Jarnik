"""Exact rational four-point selfmaps on actual shrinking Pell clusters."""
from fractions import Fraction
from math import gcd, lcm

from check_gaussian_near_real_divisor_reduction import pair_factor
from check_gaussian_reflection_replacement import gnorm
from check_mobius_actual_anchor import reduced_coefficients
from check_mobius_conformal_normalization import image, invariants
from check_mobius_conductor_transfer import least_norm


def main():
    for index in (1,11,21):
        u,v = 1,0
        for _ in range(index):
            u,v = 9*u+20*v,4*u+9*v
        blocks = [(-1,2),(u+v,2*v),(2*v,u-v),(2*u+8*v,3*u+v)]
        allocations = [(1,1,1,1),(0,1,0,0),(0,0,1,0),(0,0,0,1)]
        rows = [(1,0)]+[pair_factor(blocks,row,allocations[0]) for row in allocations[1:]]
        slopes = sorted(Fraction(y,x) for x,y in rows)
        t0,t1,t2,t3 = slopes
        c = t0+t2-t1-t3
        a = t0*t2-t1*t3
        b = c*t0*t2-a*(t0+t2)
        assert c < 0
        assert [(a*t+b)/(c*t-a) for t in slopes] == [t2,t3,t0,t1]
        h = a/c
        s2 = (h-t0)*(t2-h)
        width = t3-t0
        assert t1 < h < t2 and 0 < s2 <= width**2/4
        entries = (-a,c,b,a)  # Acts on columns (1,t).
        denominator = lcm(*(q.denominator for q in entries))
        M = tuple(int(q*denominator) for q in entries)
        content = gcd(*M)
        M = tuple(q//content for q in M)
        T,K = invariants(M)
        delta = M[0]*M[3]-M[1]*M[2]
        assert delta > 0 and K > 0 and T*width**2 >= 4*delta
        targets = image(rows,M)
        assert sorted(Fraction(y,x) for x,y in targets) == slopes
        N = least_norm(rows)
        assert least_norm(targets) == N
        cu,cv = reduced_coefficients(M)
        U,V = gnorm(cu),gnorm(cv)
        core_den = U//gcd(U,V)
        assert Fraction(V,U) == Fraction(T-2*delta,T+2*delta)
        assert core_den >= 1/width**2+Fraction(1,2)
        print(f'PASS: Pell index {index}, radius norm bits {N.bit_length()}, invariant core denominator bits {core_den.bit_length()}.')
    print('PASS: exact permutations and least-radius preservation with unbounded core on the actual four-point family.')


if __name__ == '__main__':
    main()
