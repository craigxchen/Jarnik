# Exhaustive check of X.1 on single triples.
import sys, math
from math import isqrt
def factor(n):
    f=[];p=2
    while p*p<=n:
        if n%p==0:
            f.append(p); n//=p
            if n%p==0: return None
        p+=1
    if n>1: f.append(n)
    return f
def H(a,b,c): return (a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c)
QMAX=int(sys.argv[1]); CAPF=float(sys.argv[2])
stats=dict(q=0,triples=0,nondeg=0,divfail=0,mod4fail=0,boundfail=0)
best=[]  # min ratio max/sqrt(qg) among nondegenerate admissible
for q in range(2,QMAX+1):
    ps=factor(q)
    if ps is None: continue
    stats['q']+=1
    N=int(CAPF*math.sqrt(q))
    for a in range(1,N+1):
     for b in range(a,N+1):
      for c in range(b,N+1):
        g=1; ok=True
        for p in ps:
            r=(a%p,b%p,c%p); z=sum(1 for x in r if x==0)
            if z==3: g*=p
            elif z==1:
                nz=[x for x in r if x]
                if not (nz[0]==nz[1] or (nz[0]+nz[1])%p==0): ok=False;break
            else: ok=False;break
        if not ok: continue
        stats['triples']+=1
        h=H(a,b,c)
        if h % (q*q*g*g): stats['divfail']+=1
        Q=q//g; u,v,w=a//g,b//g,c//g
        assert a%g==0 and b%g==0 and c%g==0
        hu=H(u,v,w); assert hu%(Q*Q)==0
        T=hu//(Q*Q)
        if Q%2==1 and T%4 not in (0,3): stats['mod4fail']+=1
        if Q%2==0 and T%4!=0: stats['mod4fail']+=1
        if h!=0:
            stats['nondeg']+=1
            if c*c < q*g: stats['boundfail']+=1
            best.append((c*c/(q*g),q,g,a,b,c,h//(q*q*g*g)))
print(stats)
best.sort(); print(best[:12])
