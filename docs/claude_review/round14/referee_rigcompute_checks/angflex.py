# verify an ANG-class FLEX witness at q0: level-1 ANG survivor (from onemod_angle) + random lift w/o level-2 collision;
# exact per-unit + exact per-v angle check with all-level demand; slab slack at phi = arg pi (mp, 40 digits).
import sys, json, math, random
import mpmath as mp
from reflib import *
from onemod_angle import run
name,K,q0,Fsrc=sys.argv[1],int(sys.argv[2]),int(sys.argv[3]),sys.argv[4]
res,ex=run(name,K,q0,Fsrc)
thr=json.load(open(Fsrc)); P=thr["P"]; r=len(P); F=sorted(set(thr["forced"])|set(thr.get("aborted",[]))); pis=[gprime(p) for p in P]
m=m_of(q0); n=m//4; L,iota,_=genuine_L(P,q0); k=(iota//n)%4
mp.mp.dps=40; phi=[mp.atan2(b,a) for a,b in pis]
random.seed(5)
for ell in ex:
    coll=[v for v in monos(r,K) if sum(a*b for a,b in zip(v,ell))%n==0]
    for _ in range(500):
        lift=[x+m*random.randrange(q0) for x in ell]
        if all(sum(a*b for a,b in zip(v,lift))%(n*q0)!=0 for v in coll): break
    else: print("no lift",ell); continue
    viol=0; angbad=0; slabneg=0
    for v in monos(r,K,half=True):
        Dv=Dvals(Zv(pis,v)); N=norm_v(P,v); d=sum(a*b for a,b in zip(v,lift)); dem=[]
        x=mp.fsum(t*p for t,p in zip(v,phi))
        for s in range(4):
            g=1
            for q in F:
                if q!=q0:
                    t=val(Dv[s],q)
                    if t: g*=q**t
            if d%n==0 and ((d%m)//n)*(1 if k==1 else 3)%4==s: g*=q0**(1+val(d//n,q0))
            dem.append(g); viol+= g*g>capv(P,v,s)
            slabneg+= abs(mp.sin(x-s*mp.pi/4))-(1 if s%2==0 else 1/mp.sqrt(2))*g/mp.sqrt(N) < -mp.mpf(10)**-30
        angbad+= not angle_ok_exact(dem,N)
    acc=(1<<2**(r-1))-1
    for q in F:
        if q>q0: break
        if q==q0: acc&=class_mask([x%n for x in ell],L,n,r)
        else:
            Lq,_,mq_=genuine_L(P,q); acc&=class_mask([x%(mq_//4) for x in Lq],Lq,mq_//4,r)
    print(name,K,"q0",q0,"ell",ell,"lift",lift,"| per-unit violations",viol,"| angle-incompatible v",angbad,"| slab violations at phi=arg pi",slabneg,"| non-coherent",acc==0,flush=True)
