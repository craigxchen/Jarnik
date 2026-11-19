"""Random exact tests of Part V.1: factorisation, (V.1), N/C^2<=H<=NC^2, explicit circle,
common-multiplier formula, c^4 | prod t_ij, c=2 <-> c=1."""
import random
from itertools import combinations
from math import gcd, prod
from sympy import primerange, factorint
from gauss import *

random.seed(12345)
SPLIT=[p for p in primerange(5,3000) if p%4==1]
PI={}
for p in SPLIT:
    for a in range(1,p):
        b2=p-a*a
        if b2<=0: break
        b=int(round(b2**0.5))
        if b*b==b2: PI[p]=(a,b); break
UNITS=[(1,0),(0,1),(-1,0),(0,-1)]
def gpow(a,k):
    r=(1,0)
    for _ in range(k): r=mul(r,a)
    return r
def vals(z,p):
    """(v_pi, v_pibar) of Gaussian z at split p."""
    pi=PI[p]; out=[]
    for g in (pi,conj(pi)):
        k=0; w=z
        while divides(g,w): w=exdiv(w,g); k+=1
        out.append(k)
    return out
edges=list(combinations(range(5),2))
def trial(common=None, allow2=True, balanced_flip=False):
    primes=random.sample(SPLIT,40)
    blockprimes=[primes[2*k:2*k+random.choice([1,2])] for k in range(10)]
    blockprimes=[bp[:1] if random.random()<0.5 else bp for bp in blockprimes]
    used=set(p for bp in blockprimes for p in bp)
    extra=[p for p in primes if p not in used][:6]
    # marked side per block
    sides=[]
    for (i,j) in edges:
        S={i,j}
        if balanced_flip and random.random()<0.5: S=set(range(5))-S
        sides.append(S)
    G=[(1,0)]*5
    for S,bp in zip(sides,blockprimes):
        for p in bp:
            for i in S: G[i]=mul(G[i],PI[p])
    # cofactor structure U_i
    if common is None:
        U=[]
        for i in range(5):
            u=random.choice(UNITS)
            for p in random.sample(extra+list(used),random.randint(0,3)):
                marked = any(p in bp and i in S for S,bp in zip(sides,blockprimes))
                g = PI[p] if (marked or random.random()<0.5) else conj(PI[p])
                u=mul(u,gpow(g,random.randint(1,2)))
            if allow2 and random.random()<0.3: u=mul(u,(1,1))
            U.append(u)
    else:
        U=common(sides,blockprimes,extra,used)
    V=[mul(G[i],U[i]) for i in range(5)]
    assert all(gcd(*v)==1 for v in V), "primitivity"
    n=[prod(p for S,bp in zip(sides,blockprimes) for p in bp if i in S) for i in range(5)]
    N=prod(used)
    c=[norm(V[i])//n[i] for i in range(5)]
    assert all(norm(V[i])==c[i]*n[i] for i in range(5))
    # isotropy and divisibilities
    for S,bp in zip(sides,blockprimes):
        J=prod(bp)
        for i in S: assert norm(V[i])%J==0
    t={}
    for (i,j) in edges:
        g=gcd(n[i],n[j]); d=det(V[i],V[j]); assert d%g==0 and dot(V[i],V[j])%g==0 and d!=0
        t[(i,j)]=d//g
    # (V.1)
    allp=set(used)
    for ci in c:
        for p in factorint(ci): allp.add(p)
    Hf=1
    for p in allp:
        if p==2: continue
        s=[a-b for a,b in (vals(v,p) for v in V)]
        Hf*=p**(max(s)-min(s))
    zs,H=genuine(V)
    assert H==Hf, (H,Hf)
    C=1
    for ci in c: C=C*ci//gcd(C,ci)
    assert N*H>=N*N//(C*C) if False else (H*C*C>=N and H<=N*C*C)
    # explicit circle
    B=(1,0)
    for p in used: B=mul(B,conj(PI[p]))
    pts=[]
    for i in range(5):
        Gi=G[i]; num=mul(B,mul(Gi,mul(U[i],U[i])))   # B*G_i*U_i^2 ; divide by conj(G_i) and multiply C/c_i
        num=exdiv(num,conj(Gi)); num=(num[0]*(C//c[i]),num[1]*(C//c[i]))
        pts.append(num)
    assert all(norm(z)==N*C*C for z in pts)
    # doubled ratios realised: pts_i * conj(V_i) / V_i equal for all i
    r=[mul(z,mul(conj(V[i]),conj(V[i]))) for i,z in enumerate(pts)]  # z * conj(V)^2 = z conj(V)/V *|V|^2
    base=[(x*1, y*1) for x,y in r]
    for i in range(5):
        # z_i conj(V_i)/V_i = (r_i)/|V_i|^2 must be equal
        assert base[i][0]*norm(V[0])==base[0][0]*norm(V[i]) and base[i][1]*norm(V[0])==base[0][1]*norm(V[i])
    return dict(c=c,t=t,H=H,N=N,C=C,V=V,U=U,sides=sides,blockprimes=blockprimes,used=used)

for k in range(300):
    trial(balanced_flip=(k%2==1))
print("V.1 random: factorisation, (V.1), bounds, explicit circle OK (300 trials)")

# common multiplier
def common_factory(cprimes, with2):
    def f(sides,blockprimes,extra,used):
        U=[]
        for i in range(5):
            u=random.choice(UNITS)
            for p,e in cprimes:
                marked = any(p in bp and i in S for S,bp in zip(sides,blockprimes))
                g = PI[p] if (marked or random.random()<0.5) else conj(PI[p])
                u=mul(u,gpow(g,e))
            if with2: u=mul(u,(1,1))
            U.append(u)
        return U
    return f
cnt=0
for k in range(400):
    # common c from primes: sometimes parent primes, sometimes extra
    random.seed(1000+k)
    holder={}
    def f(sides,blockprimes,extra,used):
        pool=list(used)[:2]+extra[:2]
        cp=[(p,random.randint(1,2)) for p in random.sample(pool,random.randint(1,3))]
        holder['cp']=cp
        return common_factory(cp, random.random()<0.5)(sides,blockprimes,extra,used)
    d=trial(common=f, balanced_flip=(k%2==1))
    c=d['c']; assert len(set(c))==1; c=c[0]
    T=prod(abs(x) for x in d['t'].values())
    assert T%(c**4)==0
    # H = N prod_{p in D} p^{2 v_p(c)}
    Hpred=d['N']
    for p,e in holder['cp']:
        ors=set()
        for u in d['U']:
            a,b=vals(u,p); ors.add(a>0)
        if len(ors)==2: Hpred*=p**(2*e)
    assert Hpred==d['H'],(Hpred,d['H'])
    cnt+=1
print("common multiplier: H formula and c^4 | prod|t| OK (%d trials)"%cnt)

# sharpness of c^4: search minimal exponent ratio
# c=2 <-> c=1
for k in range(200):
    random.seed(5000+k)
    d=trial(common=common_factory([],False), balanced_flip=(k%2==1))
    assert d['c']==[1]*5
    W=[mul((1,1),v) for v in d['V']]
    assert all(gcd(*w)==1 and norm(w)==2*norm(v) for w,v in zip(W,d['V']))
    for (i,j) in edges:
        assert det(W[i],W[j])==2*det(d['V'][i],d['V'][j])
    zs1,H1=genuine(d['V']); zs2,H2=genuine(W)
    assert H1==H2==d['N']
print("c=2 <-> c=1 OK; exact stars give H=N")
