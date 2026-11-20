# Exhaustive integral-cofactor configurations (one unit class) on small q_P, via the criterion
# diam(S) <= C q^{1/4} gcd_{pairs} T_ij, checked with exact integers: (diam^2)^2 <= C^4 q G^4.
import math,sys,itertools
def isprime(n): return n>1 and all(n%d for d in range(2,math.isqrt(n)+1))
def two_sq(p):
    for a in range(1,math.isqrt(p)+1):
        b=math.isqrt(p-a*a)
        if a*a+b*b==p: return (a,b)
def gm(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
P=int(sys.argv[1]); C=int(sys.argv[2])
primes=[p for p in range(P,2*P+1) if isprime(p) and p%4==1]
q=math.prod(primes); s=len(primes)
pts=[(1,0)]
for p in primes:
    a,b=two_sq(p); pts=[gm(z,f) for z in pts for f in ((a,b),(a,-b))]
n=len(pts); C4q=C**4*q
def d2(i,j): return (pts[i][0]-pts[j][0])**2+(pts[i][1]-pts[j][1])**2
T={}
adj=[set() for _ in range(n)]
for i in range(n):
    for j in range(i+1,n):
        dx=pts[j][0]-pts[i][0]; dy=pts[j][1]-pts[i][1]
        t=math.gcd(dx,dy)//2
        T[i,j]=T[j,i]=t
        if d2(i,j)**2<=C4q*t**4: adj[i].add(j); adj[j].add(i)
# genuine clusters: G=1 criterion
gen=max(1,max((2 for i in range(n) for j in adj[i] if d2(i,j)**2<=C4q),default=1))
best=2 if any(adj) else 1; bestS=None; cnt={}
def grow(S,g,dmax2,cand):
    global best,bestS
    k=len(S); cnt[k]=cnt.get(k,0)+1
    if k>best: best=k; bestS=(list(S),g)
    for c in cand:
        if c<=S[-1]: continue
        g2=g; dm=dmax2
        for x in S:
            g2=math.gcd(g2,T[x,c]); dm=max(dm,d2(x,c))
        if dm*dm<=C4q*g2**4:
            grow(S+[c],g2,dm,cand & adj[c])
for i in range(n):
    for j in adj[i]:
        if j>i: grow([i,j],T[i,j],d2(i,j),adj[i]&adj[j])
print(f"P={P} s={s} C={C} marks/class={n} pair-configs={sum(len(a) for a in adj)//2} of {n*(n-1)//2}; counts by size {dict(sorted(cnt.items()))}; max integral-cofactor M={best}; best={bestS}")
