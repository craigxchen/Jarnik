# Independent specialization: own local search, own CRT, check every condition of D_B(N), Part I eq.(1)
import sys, random
from math import gcd, lcm, prod, isqrt, log
from sympy import primerange, factorint, resultant, symbols, Poly
from mpmath import mp, mpf, log as mlog, sqrt as msqrt
from sel import SEL1, SEL2, C1, C2, Bpoly
mp.dps=60
def val(n,p):
    if n==0: return 10**9
    e=0
    while n%p==0: n//=p; e+=1
    return e
def B_at(roots,v,n):
    a,b=1,0
    for j,s in enumerate(v):
        if s==0: continue
        x,y=roots[j]; r,i=n-x,-s*y
        a,b=a*r-b*i,a*i+b*r
    return a,b
def data(roots,sel,n):
    f=[(n-x)**2+y*y for x,y in roots]
    N0=prod(f); ds=[];ss=[];ts=[]
    for v in sel:
        a,b=B_at(roots,v,n)
        g=gcd(a,b); s,t=abs(a)//g,abs(b)//g
        ds.append(s*s+t*t); ss.append(s); ts.append(t)
        assert a*a+b*b==prod(f[j] for j in range(8) if v[j])
    return f,N0,ds,ss,ts
def gain(roots,sel,n,p):
    f,N0,ds,_,_=data(roots,sel,n)
    return val(N0,p)-max(val(d,p) for d in ds), val(N0,p), all(d%2 for d in ds) if p==2 else True
def run(name,roots,sel,seed):
    bs=[abs(B_at(roots,v,0)[1]) for v in sel]
    # check b constant in n
    for n in (1,7,-13,1000):
        assert [abs(B_at(roots,v,n)[1]) for v in sel]==bs
    target=256*max(bs)**4
    t=symbols('t')
    fp=[Poly((t-x)**2+y*y,t) for x,y in roots]
    res=[int(resultant(fp[i],fp[j])) for i in range(8) for j in range(i+1,8)]
    cand=set()
    for b in bs: cand|=set(factorint(b))
    for r in res: cand|=set(factorint(abs(r)))
    cand=sorted(cand)
    print(name,'b=',bs,'target 256 max b^4 =',target,'candidate primes',cand)
    local={}
    for p in cand:
        best=(0,None,None)
        lim = p**3 if p<50 else p
        if p==2: lim=2**14
        if p in (3,5): lim=p**8
        for n in range(lim):
            gn,v0,odd=gain(roots,sel,n,p)
            if not odd: continue
            if gn>best[0]: best=(gn,n,v0)
        if best[0]>0: local[p]=best
    print('local maxima found (p: gain, n_p, v_p(N0)):',local)
    L=prod(p**g for p,(g,n,v) in local.items())
    print('L=',L,'L>target?',L>target)
    # CRT
    r,M=0,1
    for p,(g,n,v) in local.items():
        q=p**(v+1)
        r=r+M*((n-r)*pow(M,-1,q)%q); M*=q
    rnd=random.Random(seed)
    for K in (rnd.randrange(10**6,10**7), rnd.randrange(10**30,10**31)):
        n=r+M*K
        f,N0,ds,ss,ts=data(roots,sel,n)
        N=lcm(*ds); G=N0//N
        assert N0%N==0 and G%L==0
        W=mlog(mpf(N))
        ok=[]
        for d,s,t in zip(ds,ss,ts):
            c=[N%d==0, d%2==1, isqrt(d)==s, d-s*s==t*t, t>0, gcd(s,t)==1,
               256*d*d>=256*N and d*d>=256*t**4*N]
            Bmin=float(mlog(mpf(d)/msqrt(mpf(N)))/mlog(W))
            ok.append((all(c),c,Bmin, float(mpf(d)/(16*t*t*msqrt(mpf(N))))))
        print(' n has',len(str(n)),'digits; N has',len(str(N)),'digits, odd',N%2==1,'distinct d',len(set(ds))==len(ds))
        print('  offsets t',ts,'G/L=',G//L)
        for o in ok: print('  all conditions',o[0],'minimal B for upper cutoff %.3f'%o[2],'ratio d/(16t^2sqrtN) %.3g'%o[3])
        # split-supported: every prime of N divides some primitive s^2+t^2, odd -> 1 mod 4; spot check small primes
        smallp=[p for p in primerange(3,200000) if N%p==0]
        print('  small prime factors of N (<2e5):',smallp,'all 1 mod 4:',all(p%4==1 for p in smallp))
    # G bound
    Gbound=lcm(*bs)**2*prod(abs(x) for x in res)
    print(' upper bound for G: (lcm b)^2*prod Res has',len(str(Gbound)),'digits')
    Gs=[]
    for K in range(1,40):
        n=r+M*K; f,N0,ds,_,_=data(roots,sel,n); Gs.append(N0//lcm(*ds))
    print(' G along progression: min/L',min(Gs)//L,'max/L',max(Gs)//L,'distinct values',len(set(Gs)))
    return local
run('C1',C1,SEL1,1)
run('C2',C2,SEL2,2)
