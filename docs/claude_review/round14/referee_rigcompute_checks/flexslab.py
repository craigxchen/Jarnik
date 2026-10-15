# slab slack at phi = arg(pi_j) for one-modulus FLEX witnesses (genuine on F \ {q0}, witness lift at q0), norm <= K
import sys, json, math
import mpmath as mp
from reflib import *
mp.mp.dps=40
SETS=None
def run(name,K,WI=None):
    D=json.load(open(f"../rigcompute_checks/def_wit_{name}_K{K}.json")); thr=json.load(open(f"../rigcompute_checks/thr_{name}_K{K}.json"))
    P=thr["P"]; r=len(P); F=sorted(set(thr["forced"])|set(thr["aborted"])); pis=[gprime(p) for p in P]
    e=min([e for e in D if e["verdict"]=="FLEX"],key=lambda e:e["q0"]); q0=e["q0"]; m=e["m"]; n=m//4; k=(e["iota"]//n)%4
    phi=[mp.atan2(b,a) for a,b in pis]
    out=[]
    for wi,w in enumerate(e["witnesses"]):
        lift=w["lift"]; worst=mp.mpf(9); nneg=0
        for v in monos(r,K,half=True):
            Dv=Dvals(Zv(pis,v)); N=norm_v(P,v); d=sum(a*b for a,b in zip(v,lift)); x=mp.fsum(t*p for t,p in zip(v,phi))
            for s in range(4):
                g=1
                for q in F:
                    if q!=q0:
                        t=val(Dv[s],q)
                        if t: g*=q**t
                if d%n==0 and ((d%m)//n)*(1 if k==1 else 3)%4==s: g*=q0**(1+val(d//n,q0))
                sl=abs(mp.sin(x-s*mp.pi/4))-(1 if s%2==0 else 1/mp.sqrt(2))*g/mp.sqrt(N)
                worst=min(worst,sl); nneg+= sl< -mp.mpf(10)**-30
        out.append((wi,mp.nstr(worst,5),nneg))
    print(name,K,"q0",q0,"[witness, worst slab slack at phi=arg pi, #violated (v,s)]:",out,flush=True)
for a in sys.argv[1:]:
    nm,K=a.split(":"); run(nm,int(K))
