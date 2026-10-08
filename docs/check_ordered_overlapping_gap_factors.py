"""Exact checks for coprime oriented gaps on overlapping five-row windows."""

from functools import reduce
from itertools import product

from check_ordered_affine_circuit_conductor_index import (
    gconj,gdiv_exact,ggcd,gmul,gnorm,gpower,gvaluation,
)
from check_ordered_gap_opposite_twist_descent import oriented_gap


def main():
    checked=0
    for shift in (1,2,3):
        size=5+shift
        for e in range(1,4):
            for t in product(range(e+1),repeat=size):
                if min(t)!=0 or max(t)!=e:
                    continue
                h,s=oriented_gap(t[:5])
                k,v=oriented_gap(t[shift:shift+5])
                assert not (h and k and s==v)
                # Each subset may have its own nonzero common valuation.
                ta=t[:5];tb=t[shift:shift+5]
                assert oriented_gap([q-min(ta) for q in ta])==(h,s)
                assert oriented_gap([q-min(tb) for q in tb])==(k,v)
                for j in set(range(1,4)) & set(range(shift+1,shift+4)):
                    # pi-exponent of conjugate(alpha_U alpha_V).
                    pival=(h if s==-1 else 0)+(k if v==-1 else 0)
                    barval=(h if s==1 else 0)+(k if v==1 else 0)
                    assert t[j]>=pival and e-t[j]>=barval
                if shift==1:
                    # The pair gcd norm has exponent e-|t_2-t_3|.
                    assert h+k<=e-abs(t[2]-t[3])
                checked+=1

    z=((-107,-54),(-98,-69),(-91,-78),(-78,-91),(-69,-98),(-54,-107))
    assert all(gnorm(w)==14365 for w in z)
    assert gnorm(reduce(ggcd,z))==1
    for i in range(6):
        for j in range(i+1,6):
            assert z[i][0]*z[j][1]-z[i][1]*z[j][0]>0
            assert z[i][0]*z[j][0]+z[i][1]*z[j][1]>0
    alphas=[]
    for group in (z[:5],z[1:]):
        common=reduce(ggcd,group)
        normalized=[gdiv_exact(w,common) for w in group]
        alpha=(1,0)
        for pi in ((2,1),(3,2),(4,1)):
            t=[gvaluation(w,pi) for w in normalized]
            h,s=oriented_gap(t)
            alpha=gmul(alpha,gpower(pi if s==1 else gconj(pi),h))
        alphas.append(alpha)
    assert alphas==[(3,2),(3,-2)]
    assert gnorm(ggcd(*alphas))==1
    assert gmul(*alphas)==(13,0)
    pair_gcd=ggcd(z[2],z[3])
    assert gnorm(pair_gcd)==169
    for w in z[2:4]:
        assert gnorm(gdiv_exact(w,gconj(gmul(*alphas))))==85
    chord=(z[3][0]-z[2][0],z[3][1]-z[2][1])
    assert gnorm(chord)==338
    end_chord=(z[-1][0]-z[0][0],z[-1][1]-z[0][1])
    # At C=1/2 the chord squared must be <= sqrt(N)/4.
    assert 16*gnorm(end_chord)**2>14365
    print(f"PASS: {checked} allocation vectors over shifts 1,2,3")
    print("PASS: separate primitive windows, opposite norm-13 twists, and exact pair-bound equality")


if __name__=='__main__':
    main()
