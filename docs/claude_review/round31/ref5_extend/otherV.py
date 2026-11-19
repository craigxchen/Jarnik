# All-hyperplane families for random V (1 in V) with Y=t^L+theta: least L giving a c=0 cluster, D, and the theorem's bounds.
import random, itertools, math, sys
from f2lib import *
from irrlib import iota
def span(b):
    S=[0]
    for x in b: S+=[y^x for y in S]
    return S
random.seed(int(sys.argv[1])); res={}
for it in range(int(sys.argv[2])):
    n=random.randint(3,5); e=random.randint(n,9)
    b=[1]+[random.randrange(2,1<<e) for _ in range(n-1)]; V=span(b)
    if len(set(V))<len(V): continue
    co={}
    for x in range(1<<n):
        v=0
        for k in range(n):
            if x>>k&1: v^=b[k]
        co[v]=x
    forms=[c for c in range(1<<n) if c&1]
    H={c:[v for v in V if bin(co[v]&c).count('1')%2==0] for c in forms}
    maxP=0; sumi=0
    for c1,c2 in itertools.combinations(forms,2):
        J=set(H[c1])&set(H[c2]); P=1
        for v in J: P=mul(P,1^v)
        maxP=max(maxP,deg(P)); sumi+=iota(P)
    s=len(forms); L=max(2*maxP, max(deg(v) for v in V)+1)
    D=(s-1)*L   # H\{0} has s-1 elements, each factor degree L
    # verify cluster at this L directly
    Y=(1<<L,1); z={c:prod4([add4(Y,(v,0)) for v in H[c] if v]) for c in forms}
    ok=all(2*deg4(add4(z[a],z[c]))<=D for a,c in itertools.combinations(forms,2))
    zl=(1<<(L-1),1); z2={c:prod4([add4(zl,(v,0)) for v in H[c] if v]) for c in forms}
    ok2=all(2*deg4(add4(z2[a],z2[c]))<=(s-1)*(L-1) for a,c in itertools.combinations(forms,2)) if L-1>max(deg(v) for v in V) else None
    bound=3/8*s*(s-1)*(int(math.log2(s))//2)
    key=s
    r=(D, ok, ok2, 4*sumi/s, bound, [bin(x) for x in b])
    assert ok and D>=4*sumi/s>=bound-1e-9
    if key not in res or D<res[key][0]: res[key]=r
for s in sorted(res): print("s=",s,"least D=",res[s][0],"cluster",res[s][1],"L-1 cluster?",res[s][2],"(4/s)sum iota=",res[s][3],"closed bound=",res[s][4],res[s][5])
