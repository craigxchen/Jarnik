# independent verification of an 8-point cluster from the VI.3 search: V=F_2[t]_{<4}, B, T, sides
import sys, itertools
from f4lib import *
from clusters import span
from test_vi3 import coords, parity
def build(B, T, choice_bits):
    basis=[1,2,4,8]; V=span(basis); co=coords(basis)
    forms=[l for l in range(1,16) if parity(l & co[B])]
    Y=(T,B); pts=[]
    for i,l in enumerate(forms):
        c=(choice_bits>>i)&1
        S=[v for v in V if parity(l&co[v])==c]
        pts.append(prod4([add4(Y,(v,0)) for v in S]))
    return pts
def reduce_and_check(pts):
    G=pts[0]
    for z in pts[1:]: G=gcd4(G,z)
    w=[]
    for z in pts:
        q,r=divmod4(z,G); assert d4(r)<0
        w.append(monic4(q))
    D=d4(w[0]); N=mul4(w[0],conj4(w[0]))
    assert is_real(N) and all(mul4(x,conj4(x))==N for x in w) and all(d4(x)==D for x in w)
    distinct=len(set(w))==len(w)
    md=max(d4(add4(a,b)) for a,b in itertools.combinations(w,2))
    return D, md, distinct, d4(G), N
if __name__=='__main__':
    B=int(sys.argv[1]); T=int(sys.argv[2],0)
    best=None
    for bits in range(256):
        D,md,dist,g,N=reduce_and_check(build(B,T,bits))
        if dist and md<=D//2 and (best is None or D<best[0]): best=(D,md,bits,g,N)
    D,md,bits,g,N=best
    Nr=N[0]; Np=0
    for k in range(1,deg(Nr)+1,2):
        if (Nr>>k)&1: Np|=1<<(k-1)
    print(f"B={bin(B)} T={bin(T)}: best sides={bin(bits)} point degree D={D} parent deg {2*D} max pair diff {md} (<= D/2={D/2}) gcd deg {g}; N squarefree: {gcd2(Nr,Np)==1}; deg gcd(N,N')={deg(gcd2(Nr,Np))}")
