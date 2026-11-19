"""All matching-prime products for one quartet of the same 14-prime parent."""
from itertools import combinations
from math import gcd, isqrt, prod
from ii_paley import double_core
from iii_quartet import split_primes, modular_minimum, quartet_types


def norm(v):
    return sum(x*x for x in v)


def reduced_basis(a,b,c):
    # First and last coordinates parameterize the congruence; middle is exact.
    r = -c * pow(a, -1, b) % b
    v,w = (b,-a,0),(r,(-a*r-c)//b,1)
    while True:
        if norm(v)>norm(w):
            v,w=w,v
        n=norm(v)
        inner=sum(x*y for x,y in zip(v,w))
        if 2*abs(inner)<=n:
            return v,w
        k=(2*inner+n)//(2*n)
        w=tuple(y-k*x for x,y in zip(v,w))


def small_relation(a,b,c,u):
    v,w = reduced_basis(a,b,c)
    vv, ww = norm(v),norm(w)
    vw=sum(x*y for x,y in zip(v,w))
    if vv>3*u*u:
        return None
    det=vv*ww-vw*vw
    # m*v+n*w in the cube implies both inverse-Gram coefficient bounds.
    L=isqrt(3*u*u*ww//det)
    K=isqrt(3*u*u*vv//det)
    for n in range(-K,K+1):
        for m in range(-L,L+1):
            if gcd(m,n)!=1:
                continue
            xyz=tuple(m*x+n*y for x,y in zip(v,w))
            if all(0<abs(x)<=u for x in xyz):
                assert sum(x*y for x,y in zip(xyz,(a,b,c)))==0
                return xyz
    return None


def count(k):
    ps=split_primes(100000,14)
    parent=prod(ps)
    # All eight separating cuts of any pair can bear the eight largest primes.
    denom=16*isqrt(parent)
    u=(prod(ps[-8:])+denom-1)//denom
    groups=[(sum(1<<i for i in g),prod(ps[i] for i in g)) for g in combinations(range(14),k)]
    total=passing=individual=primary=0
    bound=isqrt(u)
    bounded_products={i*j for i in range(1,bound+1) for j in range(1,bound+1)}
    halved_products={i*j for i in range(1,bound//2+1) for j in range(1,bound//2+1)}
    first=None
    for i,(mask_a,a) in enumerate(groups):
        for j in range(i+1,len(groups)):
            mask_b,b=groups[j]
            if mask_b&mask_a:
                continue
            used=mask_a|mask_b
            for mask_c,c in groups[j+1:]:
                if used&mask_c:
                    continue
                total+=1
                if modular_minimum(a,b,c)>2*u*u:
                    continue
                xyz=small_relation(a,b,c,u)
                if xyz is not None:
                    passing+=1
                    assert b>2*u*u
                    individual+=all(abs(v) in bounded_products for v in xyz)
                    primary+=all(abs(v) in halved_products for v in xyz)
                    if first is None:
                        first=(a,b,c,xyz)
    print('k, number, relaxed integer passes, global U, first witness:', k,total,passing,u,first,flush=True)
    print('triples with three products factorable below the single-edge bound; same for even edges:',individual,primary,flush=True)
    return total,passing,individual,primary


def five_mark_control():
    cols=[c[1:] for family in double_core(7) for c in family]
    for five in combinations(range(7),5):
        assert sum(len(quartet_types(cols,q)[0])==3 for q in combinations(five,4))==2


if __name__=='__main__':
    five_mark_control()
    assert count(2)==(45045,156,16,12)
    assert count(3)==(560560,0,0,0)
