"""Finite rational-quotient diagnostic; not a general rational-point obstruction.

Default bound 1000 controls the maximum absolute entry of the primitive
integer representative of the critical quintet (b,-c,d,-e,-w). It does not
bound the remaining coefficients. No search beyond the stated bound is
performed. All tests use integers and Fraction arithmetic.
"""

import argparse
from fractions import Fraction as Q
from itertools import permutations
from math import gcd,isqrt


def rational_sqrt(x):
    if x<0:return None
    a,b=isqrt(x.numerator),isqrt(x.denominator)
    return Q(a,b) if a*a==x.numerator and b*b==x.denominator else None


def primitive_quintets(bound):
    """All primitive critical quintets, up to permutation/global sign.

The critical five-root theorem forces signs ++--+ in increasing absolute
order. Put s=r1+r2. The sum/cube equations give
  r3+r4=s+r5, r3*r4=s*r5+r1*r2*s/(s+r5).
Thus the integer divisibility and discriminant tests are exhaustive.
"""
    for r1 in range(1,bound):
        for r2 in range(r1+1,bound):
            ss=r1+r2;pps=r1*r2*ss
            for r5 in range(r2+3,bound+1):
                pair_sum=ss+r5
                if pps%pair_sum:continue
                pair_product=ss*r5+pps//pair_sum
                discriminant=pair_sum**2-4*pair_product
                if discriminant<=0:continue
                root=isqrt(discriminant)
                if root*root!=discriminant or (pair_sum-root)%2:continue
                r3,r4=(pair_sum-root)//2,(pair_sum+root)//2
                if not r2<r3<r4<r5:continue
                if gcd(gcd(gcd(gcd(r1,r2),r3),r4),r5)!=1:continue
                roots=(r1,r2,-r3,-r4,r5)
                assert sum(roots)==sum(a**3 for a in roots)==0
                yield roots


def search(bound):
    quintets=quadratic_passes=h2_passes=0
    for roots in primitive_quintets(bound):
        quintets+=1
        # Overall sign reversal scales every ten-factor coefficient by -1,
        # so it creates no new rational existence case.
        for b,minus_c,d,minus_e,minus_w in permutations(roots):
            c,e=-minus_c,-minus_e
            w=b-c+d-e
            assert w==-minus_w
            assert b**3-c**3+d**3-e**3-w**3==0
            h,k=c+d-e,c-d+e
            if not h:
                # h=u-v=0 is excluded by distinct absolute values.
                continue
            discriminant=Q(k*k)+Q(4*d*e*(d-e),h)
            root=rational_sqrt(discriminant)
            if root is None:continue
            quadratic_passes+=1
            for y in {(-k+root)/2,(-k-root)/2}:
                t=b+y
                if not t:continue
                u,v,q=t-b+c,t-b-d+e,2*t-b-d
                assert u**3-v**3-c**3-d**3+e**3==0
                if u**3+v**3+w**3-q**3-d**3:continue
                h2_passes+=1
                # Report quotient candidates even if their final quadratic
                # cover fails. At the recorded bounds there are none.
                print("RATIONAL QUOTIENT CANDIDATE",(t,b,c,d,e))
    return quintets,quadratic_passes,h2_passes


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bound',type=int,default=1000)
    args=parser.parse_args()
    assert args.bound>=1
    counts=search(args.bound)
    expected={300:(593,232,0),1000:(3645,444,0)}
    if args.bound in expected:assert counts==expected[args.bound]
    print(f"Bound {args.bound}: {counts[0]} primitive quintets, "
          f"{counts[1]} rational quadratic passes, {counts[2]} H2 passes")
    print("This finite diagnostic makes no assertion beyond the stated quintet-height bound.")


if __name__=='__main__':
    main()
