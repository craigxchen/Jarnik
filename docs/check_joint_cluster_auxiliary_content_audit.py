"""Exact content checks for chord products and complementary stacked minors."""

from fractions import Fraction as Q
from functools import reduce
from itertools import combinations,combinations_with_replacement,permutations,product
from math import gcd,prod
from random import Random

from check_ordered_affine_circuit_conductor_index import (
    gconj,gdiv_exact,ggcd,gmul,gnorm,gpower,gvaluation,
)
from check_general_bipartition_opposite_twist_descent import cuts,oriented_gap

PERMS=tuple(permutations(range(4)))
SIGNS=tuple((-1)**sum(p[i]>p[j] for i in range(4) for j in range(i+1,4)) for p in PERMS)


def det4(rows):
    return sum(s*prod(rows[i][p[i]] for i in range(4)) for p,s in zip(PERMS,SIGNS))


def line(a,b):
    dx,dy=b[0]-a[0],b[1]-a[1]
    return (-dy,dx,dy*a[0]-dx*a[1])


def evaluate(coef,z):return coef[0]*z[0]+coef[1]*z[1]+coef[2]


def divide(z,a):
    return tuple(Q(x,gnorm(a)) for x in gmul(z,gconj(a)))


def valuation(n,p):
    assert n
    n=abs(n);v=0
    while n%p==0:n//=p;v+=1
    return v


def forced(e,t,u,v):
    order=sorted(range(4),key=lambda i:t[i])
    t,u,v=([a[i] for i in order] for a in (t,u,v))
    extra=min(u[1]-u[0],v[1]-v[0])+min(u[3]-u[2],v[3]-v[2])
    return e-(t[3]-t[0])+extra


def main():
    assignments=0
    for e in range(1,7):
        for t in combinations_with_replacement(range(e+1),4):
            for bits in product((0,1),repeat=e):
                a=sum(bits);b=e-a
                u=[sum(bits[:q]) for q in t];v=[q-r for q,r in zip(t,u)]
                rows=[(u[i],a-u[i],v[i],b-v[i]) for i in range(4)]
                actual=min(sum(rows[i][p[i]] for i in range(4)) for p in PERMS)
                assert actual==forced(e,t,u,v)
                assignments+=1
    rng=Random(923917)
    nonzero=0
    for _ in range(300):
        # Ambient six rows contain every level, so each prime threshold
        # has width one. Select its unordered labelled cut consistently
        # across all primes, including coincidences between prime cuts.
        selected={min(mask,63^mask):rng.randrange(2) for mask in range(1,32)}
        beta=[(1,0)]*6;w=[rng.choice(((1,0),(0,1),(-1,0),(0,-1))) for _ in range(6)]
        prime_data=[];n=qq=1
        for pi in ((2,1),(3,2),(4,1)):
            e=5;t=list(range(6));rng.shuffle(t)
            masks=[sum(1<<i for i in range(6) if t[i]>=h) for h in range(1,e+1)]
            bits=[selected[min(mask,63^mask)] for mask in masks]
            a=sum(bits);b=e-a
            u=[sum(bits[:q]) for q in t];v=[q-r for q,r in zip(t,u)]
            prime_data.append((gnorm(pi),e,t,u,v))
            n*=gnorm(pi)**e;qq*=gnorm(pi)**a
            for i in range(6):
                beta[i]=gmul(beta[i],gmul(gpower(pi,u[i]),gpower(gconj(pi),a-u[i])))
                w[i]=gmul(w[i],gmul(gpower(pi,v[i]),gpower(gconj(pi),b-v[i])))
        z=[gmul(a,b) for a,b in zip(beta,w)]
        assert all(gnorm(a)==qq for a in beta)
        assert all(gnorm(a)==n//qq for a in w)
        assert all(gnorm(a)==n for a in z)
        assert gnorm(reduce(ggcd,beta))==gnorm(reduce(ggcd,w))==1
        rows=rng.sample(range(6),4)
        dd=det4([beta[i]+w[i] for i in rows])
        divisor=prod(p**forced(e,[t[i] for i in rows],[u[i] for i in rows],[v[i] for i in rows])
                     for p,e,t,u,v in prime_data)
        assert dd%divisor==0
        nonzero+=dd!=0

    # Maximal-gap chord-product fixture, including line coefficient content.
    z=((4,-33),(9,-32),(12,-31),(24,-23),(31,-12))
    side={0,2};alpha=(4,1);k=17;n=1105
    w=tuple(gdiv_exact(a,alpha if i in side else gconj(alpha)) for i,a in enumerate(z))
    assert w==((-1,-8),(4,-7),(1,-8),(7,-4),(8,-1))
    assert len(set(w))==5 and gnorm(reduce(ggcd,w))==1
    assert all(gnorm(a)==65 for a in w)
    ks=[]
    for s,t in cuts(5):
        value=1
        for pi in ((2,1),(3,2),(4,1)):
            ts=[gvaluation(a,pi) for a in z]
            value*=gnorm(pi)**oriented_gap(ts,s,t)[0]
        ks.append(value)
    assert max(ks)==17
    ls,lt=line(z[0],z[2]),line(z[1],z[3])
    cs,ct=reduce(gcd,map(abs,ls)),reduce(gcd,map(abs,lt))
    assert (cs,ct)==(2,3)
    assert tuple(a//cs for a in ls)==(-1,4,136)
    assert tuple(a//ct for a in lt)==(-3,5,187)
    val=evaluate(ls,z[4])*evaluate(lt,z[4])//(cs*ct)
    assert val==1938 and valuation(val,17)==1
    lx,ly=line(w[0],w[2]),line(w[1],w[3])
    # Three affine-independent arguments verify each affine pullback identity.
    for zz in ((0,0),(1,0),(0,1)):
        assert evaluate(ls,zz)==k*evaluate(lx,divide(zz,alpha))
        assert evaluate(lt,zz)==k*evaluate(ly,divide(zz,gconj(alpha)))
    assert all(z[i][0]*z[j][1]-z[i][1]*z[j][0]>0 for i in range(5) for j in range(i+1,5))
    chord=tuple(z[-1][j]-z[0][j] for j in range(2))
    assert 16*gnorm(chord)**2>n  # incompatible with C=1/2

    # Exact half-norm sharpness, with complete selected cut blocks.
    pi=(4,1);p=17;levels=(0,1,2,3,6)
    beta=[];w=[];z=[]
    for t in levels:
        u=min(t,3);v=t-u
        beta.append(gmul(gpower(pi,u),gpower(gconj(pi),3-u)))
        w.append(gmul(gpower(pi,v),gpower(gconj(pi),3-v)))
        z.append(gmul(beta[-1],w[-1]))
    assert all(gnorm(a)==p**3 for a in beta+w)
    assert all(gnorm(a)==p**6 for a in z)
    assert gnorm(reduce(ggcd,beta))==gnorm(reduce(ggcd,w))==1
    dd=det4([beta[i]+w[i] for i in (0,1,2,4)])
    assert dd==2659072 and valuation(dd,p)==1
    assert forced(6,[0,1,2,6],[0,1,2,3],[0,0,0,3])==1
    assert gnorm(reduce(ggcd,[z[i] for i in (0,1,2,4)]))==1
    assert gpower(pi,3)==(52,47)
    chord=tuple(z[-1][j]-z[0][j] for j in range(2))
    assert 16*gnorm(chord)**2>p**6
    assert nonzero
    print(f"PASS: {assignments} exact assignment minima; 300 multi-prime tuples ({nonzero} nonzero minors)")
    print("PASS: primitive maximal-gap line-product fixture and exact half-norm content sharpness")


if __name__=='__main__':
    main()
