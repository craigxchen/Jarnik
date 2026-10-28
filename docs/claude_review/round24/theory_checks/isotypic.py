"""Exact max vanishing order at e of iota-even / iota-odd sections of m(-K) on LM_5 (iota: x -> x^{-1}).
S_m = {s in Z^5: sum s = 0, sum|s| <= 2m}. Even part basis: x^s + x^{-s} (orbits {s,-s}, plus s=0);
odd part: x^s - x^{-s}. With x = exp(y), y_5 = 0: even f has Taylor terms only in even degrees j,
 coefficient vector of degree j: sum_O c_O * s_O^alpha (|alpha| = j), up to multinomial factors.
t_plus(m) = least even d such that the map c -> (sum_O c_O s_O^alpha)_{|alpha| even <= d} is injective.
(Injective at d means: no nonzero even f has order > d; order of even f is even, so t_plus <= d;
 non-injective at d-2 means some even f has order >= d-1, hence >= d.)
Rank certification: rank computed mod p = 2^31-1 on random linear combinations of the monomial rows
(evaluation at random points y_r); full column rank mod p => full rank over Q => injective (certified).
Non-injectivity (lower bound) is certified separately by explicit sections (sections.py), not here.
"""
import sys, random, itertools
import numpy as np
P = 2147483647
def S(m):
    out=[]
    rng=range(-2*m,2*m+1)
    for s in itertools.product(rng,repeat=4):
        s5=-sum(s); t=s+(s5,)
        if sum(abs(v) for v in t)<=2*m: out.append(t)
    return out
def orbits(Sm):
    seen=set(); orb=[]
    for s in Sm:
        if s in seen: continue
        ms=tuple(-v for v in s); seen.add(s); seen.add(ms); orb.append(s)
    return orb
def rank_mod(A):
    A=A.copy()%P; r=0; rows,cols=A.shape
    for c in range(cols):
        piv=None
        for i in range(r,rows):
            if A[i,c]: piv=i;break
        if piv is None: continue
        A[[r,piv]]=A[[piv,r]]
        inv=pow(int(A[r,c]),P-2,P)
        A[r]=(A[r]*inv)%P
        col=A[:,c].copy(); col[r]=0
        nz=np.nonzero(col)[0]
        if len(nz): A[nz]=(A[nz]-np.outer(col[nz],A[r]))%P
        r+=1
        if r==rows: break
    return r
def tpar(m,parity,seed=1):
    Sm=S(m); orb=orbits(Sm)
    if parity==1: orb=[s for s in orb if any(s)]  # odd part excludes s=0
    n=len(orb)
    random.seed(seed)
    d=parity  # degrees of this parity
    while True:
        degs=[j for j in range(parity,d+1,2)]
        # number of monomials of degree j in 4 vars
        rows=[]
        for j in degs:
            nmj=(j+3)*(j+2)*(j+1)//6
            for r in range(min(nmj,n)):
                y=[random.randrange(P) for _ in range(4)]
                rows.append([pow((sum(s[i]*y[i] for i in range(4)))%P, j, P) for s in orb])
        A=np.array(rows,dtype=np.int64)
        rk=rank_mod(A.T) if A.shape[0]>A.shape[1] else rank_mod(A)
        if rk==n: return d,n
        d+=2
if __name__=="__main__":
    M=int(sys.argv[1]) if len(sys.argv)>1 else 4
    for m in range(1,M+1):
        tp,npl=tpar(m,0); tm,nmi=tpar(m,1)
        print(f"m={m} |S_m|={len(S(m))} dim_even={npl} dim_odd={nmi} t_plus={tp} t_minus={tm} max={max(tp,tm)} floor(10m/3)={10*m//3}",flush=True)
