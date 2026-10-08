import itertools, math, random, sys
from fractions import Fraction
from math import gcd
from rlib import *
# (2) actual clusters
random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 3)
def lcm(a,b): return a*b//gcd(a,b)
def build(Nfac):
    N=1
    for p,e in Nfac: N*=p**e
    # generate all points of norm N from Gaussian factorizations
    pts=[(1,0)]
    for p,e in Nfac:
        pi=gauss_prime_above(p); cp=gconj(pi)
        new=[]
        for a in range(e+1):
            f=(1,0)
            for _ in range(a): f=gmul(f,pi)
            for _ in range(e-a): f=gmul(f,cp)
            for z in pts: new.append(gmul(z,f))
        pts=new
    allp=[]
    for z in pts:
        u=z
        for _ in range(4):
            allp.append(u); u=gmul(u,(0,1))
    return N, sorted(set(allp), key=angle)
results=[]
for trial in range(40):
    fac=[(5,random.randint(6,10))]+[(p,random.randint(1,2)) for p in random.sample([13,17,29,37,41],random.randint(1,3))]
    N,pts=build(fac)
    # restrict to an open half-plane window
    st=random.randrange(len(pts)); win=[]
    a0=angle(pts[st])
    for t in range(len(pts)):
        z=pts[(st+t)%len(pts)]
        if (angle(z)-a0)%(2*math.pi) < math.pi*0.95: win.append(z)
    random.shuffle(win)
    sub=[]
    for z in win:
        okz=True
        for w in sub:
            c=cot_half(w,z,N)
            if c.denominator%5==0: okz=False; break
        if okz: sub.append(z)
    P,g=primitive_cluster(sub)
    Np=gnorm(P[0]); M=len(P)
    if M<17: continue
    pairs=list(itertools.combinations(range(M),2))
    L0=1
    for (i,j) in pairs: L0=lcm(L0,cot_half(P[i],P[j],Np).denominator)
    assert L0%5!=0
    fN=factor(Np); W=math.log(Np)
    C=math.sqrt(max(gnorm((P[i][0]-P[j][0],P[i][1]-P[j][1])) for (i,j) in pairs))/Np**0.25
    lhsA=0; sumD=0; info=[]
    for p,e in fN.items():
        pi=gauss_prime_above(p); a=[vpi(z,pi) for z in P]
        v=0;t=L0
        while t%p==0: t//=p; v+=1
        c=(p-1)*p**v
        cnt=[a.count(h) for h in range(e+1)]
        assert max(cnt)<=c, (p,cnt,c)
        D=sum((sum(1 for x in a if x<tau)-M/2)**2 for tau in range(1,e+1))
        sumD+=D*math.log(p)
        term=math.log(p)*max(M-4*c,0)**3/(12*c); lhsA+=term
        if p==5: info=(e,cnt,D,max(M-4*c,0)**3/(12*c))
    rhsA=M*W/4+len(pairs)*math.log(C*C/2)
    assert lhsA<=sumD+1e-9<=rhsA+1e-6, (lhsA,sumD,rhsA)
    results.append((M,round(lhsA,2),round(sumD,2),round(rhsA,2),info))
for r in results[:12]: print(r)
print('Theorem A verified with active cubic term on',sum(1 for r in results if r[1]>0),'actual clusters')
