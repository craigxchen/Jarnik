# Judge [J]: independent exact check of the five-mark Pell family of VII.5.
# Own arithmetic: polynomials over Q(i) as lists of (Fraction re, Fraction im), no sympy.
from fractions import Fraction as Fr
import sys, cmath
def cm(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def ca(a,b): return (a[0]+b[0], a[1]+b[1])
def cn(a): return (-a[0],-a[1])
def cinv(a):
    d=a[0]*a[0]+a[1]*a[1]; return (a[0]/d, -a[1]/d)
Z=(Fr(0),Fr(0))
def trim(p):
    p=list(p)
    while p and p[-1]==Z: p.pop()
    return p
def padd(p,q):
    n=max(len(p),len(q)); r=[]
    for k in range(n):
        a=p[k] if k<len(p) else Z; b=q[k] if k<len(q) else Z
        r.append(ca(a,b))
    return trim(r)
def pneg(p): return [cn(c) for c in p]
def pmul(p,q):
    if not p or not q: return []
    r=[Z]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        if a==Z: continue
        for j,b in enumerate(q):
            r[i+j]=ca(r[i+j],cm(a,b))
    return trim(r)
def pbar(p): return [(c[0],-c[1]) for c in p]
def deg(p): return len(p)-1
def pmod(p,q):
    p=list(p); iq=cinv(q[-1])
    while len(p)>=len(q):
        c=cm(p[-1],iq); s=len(p)-len(q)
        for k,b in enumerate(q): p[s+k]=ca(p[s+k],cn(cm(c,b)))
        p=trim(p)
    return p
def pgcd(p,q):
    while q: p,q=q,pmod(p,q)
    ic=cinv(p[-1]); return [cm(c,ic) for c in p]
def C(x,y=0): return (Fr(x),Fr(y))
t=[Z,C(1)]
def const(c): return [c]
def blocks(n,seeds):
    out=[]
    for s in seeds:
        X=[C(1)]; Y=[const(C(s))][0]
        for _ in range(n):
            # A = [[t^2+1, t],[t,1]]
            X,Y = padd(pmul([C(1),Z,C(1)],X), pmul(t,Y)), padd(pmul(t,X),Y)
        out.append(padd(X, pmul([C(0,1)],Y)))
    return out
def marks(G,edges):
    beta=[C(0,2),C(1)]
    def prod(L):
        r=[C(1)]
        for f in L: r=pmul(r,f)
        return r
    z0=pmul(beta,prod([pbar(g) for g in G]))
    zs=[z0]
    for e in edges:
        zs.append(pmul(pbar(beta),prod([G[j] if j in e else pbar(G[j]) for j in range(len(G))])))
    return zs
def ev(p,x): 
    r=0
    for c in reversed(p): r=r*x+complex(float(c[0]),float(c[1]))
    return r
edges=[(0,1),(0,2),(1,2),(0,3)]
for seeds in [(Fr(1,2),Fr(-2),Fr(7,3),Fr(4)), (Fr(-1,3),Fr(5,2),Fr(6),Fr(-7))]:
  for n in (1,2,3):
    G=blocks(n,seeds); zs=marks(G,edges)
    N=pmul(zs[0],pbar(zs[0]))
    assert all(pmul(z,pbar(z))==N for z in zs), "norms"
    g=zs[0]
    for z in zs[1:]: g=pgcd(g,z)
    assert deg(g)==0, "gcd"
    D={deg(z) for z in zs}; assert len(D)==1; D=D.pop()
    lead={zs[k][-1] for k in range(5)}
    diffs=[deg(padd(zs[a],pneg(zs[b]))) for a in range(5) for b in range(a+1,5)]
    assert min(diffs)>=0, "distinct"
    # contact with marks scaled to common leading coefficient (all monic here?)
    h=D-max(diffs)
    # squarefree check of N
    dN=[ (c[0]*k, c[1]*k) for k,c in enumerate(N)][1:]
    sqf = deg(pgcd(N,dN))==0
    # Jacobian of seeds -> anchored ratios at t=1 (float, numerical rank)
    def ratios(sd):
        G=blocks(n,[Fr(x).limit_denominator(10**12) if not isinstance(x,Fr) else x for x in sd]); z=marks(G,edges)
        return [ev(z[k],1.0)/ev(z[0],1.0) for k in range(1,5)]
    import numpy as np
    base=ratios(seeds); eps=Fr(1,10**6); J=[]
    for k in range(4):
        sd=list(seeds); sd[k]=sd[k]+eps
        r=ratios(sd); J.append([ (cmath.phase(r[j]/base[j]))/1e-6 for j in range(4)])
    sv=np.linalg.svd(np.array(J),compute_uv=False)
    print(f"seeds={tuple(str(s) for s in seeds)} n={n} D={D} h={h} leads={sorted(lead)} pairdeg={sorted(diffs)} N_squarefree={sqf} |z(1)|-phase-Jacobian sv={['%.3g'%x for x in sv]}")
