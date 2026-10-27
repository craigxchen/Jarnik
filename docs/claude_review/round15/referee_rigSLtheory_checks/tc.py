# Referee's independent TC (triple condition) certificate search for Theorem P.
import sys, itertools, time, rlib
from math import gcd
def split_primes(r):
    out=[];p=5
    while len(out)<r:
        if rlib.isprime(p) and p%4==1: out.append(p)
        p+=1
    return out
class Lat:
    def __init__(s,r): s.r=r; s.piv={}
    def add(s,v):
        v=list(v)
        for c in range(s.r):
            if v[c]==0: continue
            if c not in s.piv:
                if v[c]<0: v=[-x for x in v]
                s.piv[c]=v; return True
            b=s.piv[c]
            if v[c]%b[c]==0:
                f=v[c]//b[c]; v=[x-f*y for x,y in zip(v,b)]; continue
            # extended gcd combine
            a,bb=b[c],v[c]
            g,x,y=egcd(a,bb)
            new=[x*p+y*q for p,q in zip(b,v)]
            v=[(a//g)*q-(bb//g)*p for p,q in zip(b,v)]  # has 0 at c
            if new[c]<0: new=[-t for t in new]
            s.piv[c]=new
            # v now zero at c; continue
        return False
    def index(s):
        if len(s.piv)<s.r: return 0
        d=1
        for c in range(s.r): d*=abs(s.piv[c][c])
        return d
def egcd(a,b):
    if b==0: return (abs(a), 1 if a>0 else -1, 0)
    g,x,y=egcd(b,a%b)
    return (g,y,x-(a//b)*y)
def smooth(n,Y,plist):
    n=abs(n)
    for q in plist:
        while n%q==0: n//=q
    return n==1
def run(r,mode,Y=None,tlimit=1800):
    P=split_primes(r); pis=[rlib.gpi(p) for p in P]; Nrad=1
    for p in P: Nrad*=p
    plist=[q for q in range(3,(Y or 3)+1) if rlib.isprime(q) and Nrad%q!=0]
    cache={}
    def inP0(u):
        if u in cache: return cache[u]
        z=rlib.Zv(pis,u); re=z[0]*z[0]-z[1]*z[1]
        ok = gcd(re,Nrad)==1 and (mode=='all' or smooth(re,Y,plist))
        cache[u]=ok; return ok
    if mode=='all': base=None; maxsum=3
    else: maxsum=4
    # generate base vectors: V_2 (Y mode) or V_3 small r (all mode); for all-mode large r use V_2 base with sums in V_3? (author: triples inside V_3)
    if mode=='all':
        base=rlib.VKfast(r,3) if r<=16 else rlib.VKfast(r,2)
        maxsum=3
    else:
        base=rlib.VKfast(r,2)
    good=[v for v in base if inP0(v)]
    L=Lat(r); t0=time.time(); ntrip=0
    for i,v in enumerate(good):
        for w in good[i:]:
            u=tuple(a+b for a,b in zip(v,w))
            if not any(u) or sum(map(abs,u))>maxsum: continue
            if inP0(u):
                ntrip+=1
                L.add(v);L.add(w);L.add(u)
        if L.index()==1: break
        if time.time()-t0>tlimit: break
        if i%500==0: print('  progress',r,mode,Y,i,len(good),'idx',L.index(),'rank',len(L.piv),round(time.time()-t0),file=sys.stderr,flush=True)
    return dict(r=r,mode=mode,Y=Y,nbase_good=len(good),ntriples_used=ntrip,index=L.index(),rank=len(L.piv),time=round(time.time()-t0,1))
if __name__=='__main__':
    r=int(sys.argv[1]); mode=sys.argv[2]; Y=int(sys.argv[3]) if len(sys.argv)>3 else None
    print(run(r,mode,Y),flush=True)
