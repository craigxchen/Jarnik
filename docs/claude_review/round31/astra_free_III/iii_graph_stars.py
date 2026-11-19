"""Exact 19-mark partial-atom integer lift of the Paley cofactor model."""
from fractions import Fraction
from itertools import combinations
from math import gcd,isqrt,prod
from iii_graph_rounding import prime,skew_core


def sqrt_minus_one(p):
    u=2
    while pow(u,(p-1)//2,p)==1:
        u+=1
    r=pow(u,(p-1)//4,p)
    r2=(r-p*((r*r+1)//p)*pow(2*r,-1,p))%(p*p)
    assert (1+r2*r2)%(p*p)==0
    return r2


def main():
    M=20
    rows=[[int(x==1) for x in row] for row in skew_core(M-1)[1:]]
    p=[]
    n=2**40+1
    while len(p)<2*(M-1):
        if prime(n):
            p.append(n)
        n+=4
    # Eleven globally selected first-family blocks. Every pair shares at least one;
    # all 19 vertices have <=7 such obligations.
    S=[2*k for k in (0,1,2,3,4,5,6,7,8,9,15)]
    roots={j:sqrt_minus_one(p[j]) for j in S}
    q=2**34
    assert all((q*q+1)%v for v in p)
    G=p[0]**(M//2-1)
    H=isqrt(p[0]**(M-1))
    b=[4*(2*M*M+i)**2 for i in range(len(rows))]
    stars=[]
    largest_error=Fraction(0)
    for i,row in enumerate(rows):
        selected=[j for j in S if row[j]]
        assert len(selected)<=7
        mod=prod(p[j] for j in selected)
        s_center=Fraction(q*G*(b[i]-b[0]),H)
        s=2*round(s_center/2)
        while s==0 or gcd(s,mod)!=1:
            s+=2
        residue,step=1,abs(s)
        for j in selected:
            pj=p[j]
            r=roots[j]%pj
            dest=(-r*s*pow(1+r*q,-1,pj))%pj
            t=((dest-residue)*pow(step,-1,pj))%pj
            residue+=step*t
            step*=pj
        t0=(H//q-residue)//step
        for k in range(2*len(S)+2):
            x=residue+(t0+k)*step
            y=q*x+s
            if all((x-roots[j]*y)%p[j] and
                   bool((x+roots[j]*y)%p[j])==(not row[j]) and
                   (x+roots[j]*y)%(p[j]*p[j]) for j in S):
                break
        else:
            raise AssertionError("No lift found")
        assert gcd(x,y)==1 and x%2 and y%2==0
        B=prod(v for a,v in zip(row,p) if a)
        assert abs(x*x+y*y-B)*1000000<B
        # squared distance to the common rational line y=qx
        assert (y-q*x)**2 < 4*(q*q+1)
        stars.append((x,y))
    counts=[]
    for i,j in combinations(range(len(rows)),2):
        shared=[l for l in range(len(p)) if rows[i][l] and rows[j][l]]
        g=prod(p[l] for l in shared)
        forced=prod(p[l] for l in shared if l in S)
        assert len(shared)==9 and forced>1
        x,y=stars[i]
        u,v=stars[j]
        det=x*v-y*u
        assert det%forced==0
        largest_error=max(largest_error,abs(Fraction(det,g)-(b[j]-b[i])))
        counts.append(sum(1 for l in shared if l in S))
    assert largest_error<Fraction(1,64)
    assert any((x*v-y*u)%prod(p[l] for l in range(len(p)) if rows[i][l] and rows[j][l])
               for i,(x,y) in enumerate(stars) for j,(u,v) in enumerate(stars) if i<j)
    print('19 primitive primary-type stars, 38 prime atoms, 11 globally exact.')
    print('171 weighted scalar errors < 1/',largest_error.denominator//largest_error.numerator,
          '; exact shared atoms per pair:',min(counts), 'to',max(counts),'of 9.')
    print('All 19 norm relative errors < 1e-6 and every perpendicular distance < 2.')


if __name__=='__main__':
    main()
