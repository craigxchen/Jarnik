import itertools
from f4lib import *
from clusters import span
from test_vi3 import coords, parity
basis=[1,2,4]; V=span(basis); co=coords(basis)
def show(B,T):
    forms=[l for l in range(1,8) if parity(l&co[B])]
    Y=(T,B); res=[]
    for bits in range(16):
        pts=[]
        for i,l in enumerate(forms):
            S=[v for v in V if parity(l&co[v])==((bits>>i)&1)]
            pts.append((S,prod4([add4(Y,(v,0)) for v in S])))
        G=pts[0][1]
        for _,z in pts[1:]: G=gcd4(G,z)
        w=[monic4(divmod4(z,G)[0]) for _,z in pts]
        D=d4(w[0]); N=mul4(w[0],conj4(w[0]))
        assert all(mul4(x,conj4(x))==N for x in w) and is_real(N)
        if len(set(w))<4: continue
        md=max(d4(add4(a,b)) for a,b in itertools.combinations(w,2))
        if md<=D//2: res.append((D,bits,md,d4(G),[S for S,_ in pts],w,N))
    res.sort(key=lambda r:r[0]); return res[0] if res else None
r=show(1,0)
D,bits,md,g,Ss,w,N=r
print("B=1,T=0: D'=",D,"maxdiff",md,"gcd deg",g,"sides",Ss)
def fmt(x):
    A,Bb=x; n=d4(x); cs=[]
    for k in range(n,-1,-1):
        c=((A>>k)&1)|(((Bb>>k)&1)<<1); cs.append(".1ab"[c])  # a=theta, b=theta^2
    return "".join(cs)
for x in w: print("  point (coeffs high->low, a=theta,b=theta^2):",fmt(x))
print("  N =",bin(N[0]))
