import json, math, random
from reflib import *
thr=json.load(open('../rigcompute_checks/thr_first8_K4.json'))
P=thr['P']; r=len(P); K=4; F=sorted(set(thr['forced'])|set(thr['aborted'])); pis=[gprime(p) for p in P]
q0=179; m=m_of(q0); n=m//4
ell=[2, 18, 48, 62, 103, 173, 145, 168]
nu=nu_par(P,q0); assert all(x%2==p for x,p in zip(ell,nu))
L,iota,_=genuine_L(P,q0)
print("my iota",iota,"n",n)
V=monos(r,K)
coll=[v for v in V if sum(a*b for a,b in zip(v,ell))%n==0]
print("collisions",len(coll))
random.seed(1)
for trial in range(1000):
    lift=[x+m*random.randrange(q0) for x in ell]
    if all(sum(a*b for a,b in zip(v,lift))%(n*q0)!=0 for v in coll): break
print("lift",lift)
# both iota conventions (author's file not available): test both
for k in (1,3):
    nviol=0
    for v in V:
        z=Zv(pis,v); Dv=Dvals(z); d=sum(a*b for a,b in zip(v,lift))
        for s in range(4):
            dem=1
            for qq in F:
                if qq!=q0:
                    t=val(Dv[s],qq)
                    if t: dem*=qq**t
            if d%n==0:
                sd=((d%m)//n)*(1 if k==1 else 3)%4
                if sd==s: dem*=q0**(1+val(d//n,q0))
            if dem*dem>capv(P,v,s): nviol+=1
    print("iota=",k,"n: violations",nviol)
acc=(1<<2**(r-1))-1
for q in F:
    if q>q0: break
    if q==q0: acc&=class_mask([x%n for x in ell],L,n,r)
    else:
        Lq,_,mm=genuine_L(P,q); acc&=class_mask([x%(mm//4) for x in Lq],Lq,mm//4,r)
print("mask AND on F∩[3,179]",acc)
import mpmath as mp
mp.mp.dps=40
phi=[mp.atan2(b,a) for a,b in pis]; nneg=0; worst=mp.mpf(9); k=3
for v in monos(r,K,half=True):
    Dv=Dvals(Zv(pis,v)); N=norm_v(P,v); d=sum(a*b for a,b in zip(v,lift)); x=mp.fsum(t*p for t,p in zip(v,phi))
    for s in range(4):
        g=1
        for q in F:
            if q!=q0:
                t=val(Dv[s],q)
                if t: g*=q**t
        if d%n==0 and ((d%m)//n)*3%4==s: g*=q0**(1+val(d//n,q0))
        sl=abs(mp.sin(x-s*mp.pi/4))-(1 if s%2==0 else 1/mp.sqrt(2))*g/mp.sqrt(N)
        worst=min(worst,sl); nneg+= sl< -mp.mpf(10)**-30
print("first8 179 slab at phi=arg pi: violated",nneg,"worst",mp.nstr(worst,5))
