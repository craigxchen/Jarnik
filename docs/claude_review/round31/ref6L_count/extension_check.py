"""(d) check |q_i-q_j|>=exp(-d_ij/2), sum d_ij = 6W-2W_singleton, and the deletion profiles,
on random primitive Gaussian 5-tuples (exact arithmetic)."""
import random, math, itertools
from glib import *
from sympy import factorint, primerange
random.seed(7)
P1=[p for p in primerange(5,400) if p%4==1]
def cuts_of(z,H):
    """return dict frozenset(high side)->log weight, plus W"""
    fac=factorint(H); cuts={}
    for p,e in fac.items():
        pi=split_prime(p); a=[val(x,pi) for x in z]
        assert min(a)==0 and max(a)==e, (p,a,e)
        for t in range(1,e+1):
            hi=frozenset(k for k in range(len(z)) if a[k]>=t)
            lo=frozenset(range(len(z)))-hi
            key=min(hi,lo,key=lambda s:(len(s),sorted(s)))
            if len(hi)==len(lo): key=hi if 0 in hi else lo
            cuts[key]=cuts.get(key,0)+math.log(p)
    return cuts
bad=0;N=0
for trial in range(300):
    ps=random.sample(P1,random.randint(2,5))
    pis=[split_prime(p) for p in ps]; es=[random.randint(1,4) for _ in ps]
    z=[]
    for k in range(5):
        x=(1,0)
        for pi,e in zip(pis,es):
            a=random.randint(0,e)
            for _ in range(a): x=mul(x,pi)
            for _ in range(e-a): x=mul(x,conj(pi))
        z.append(mul(x,random.choice([(1,0),(0,1),(-1,0),(0,-1)])))
    if len(set(z))<5: continue
    g=z[0]
    for x in z[1:]: g=gcd(g,x)
    z=[exdiv(x,g) for x in z]; H=nrm(z[0]); W=math.log(H)
    cuts=cuts_of(z,H)
    Ws=sum(w for c,w in cuts.items() if len(c)==1)
    D=0
    for i,j in itertools.combinations(range(5),2):
        dij=sum(w for c,w in cuts.items() if (i in c)!=(j in c))
        D+=dij
        q=abs(complex(*z[i])-complex(*z[j]))/math.sqrt(H)
        if q< math.exp(-dij/2)*(1-1e-12): bad+=1
        # reduced word: z_j/z_i = xi C/conj C with |C|^2=e^{d_ij}
        C2=H//nrm(gcd(z[i],z[j]))
        assert abs(math.log(C2) - dij)<1e-9*W+1e-12
    assert abs(D-(6*W-2*Ws))<1e-8*W+1e-9,(D,6*W-2*Ws)
    N+=1
print("tuples",N,"violations of |q_i-q_j|>=e^{-d/2}:",bad,"; sum identity held on all")
