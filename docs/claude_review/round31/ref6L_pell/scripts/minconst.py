# Smallest limiting constant K = lim arc/R'^(1/3) over seeds for the rate (1,1,1) Pell family.
import math, sys
from itertools import combinations
from pell4 import *
c8,s8=math.cos(math.pi/8),math.sin(math.pi/8)
B=int(sys.argv[1]) if len(sys.argv)>1 else 8
V=[]
for x in range(0,B+1):
    for y in range(-B,B+1):
        if math.gcd(x,y)==1 and (x>0 or y>0): V.append((x,y))
def ab(v): return v[0]*c8+v[1]*s8, -v[0]*s8+v[1]*c8
AB=[ab(v) for v in V]
best=[]
nexact=10
LIM=float(sys.argv[2]) if len(sys.argv)>2 else 10.0
for i,j,k in combinations(range(len(V)),3):
    (a1,b1),(a2,b2),(a3,b3)=AB[i],AB[j],AB[k]
    u1,u2,u3=b1/a1,b2/a2,b3/a3
    S=[0,u1+u2,u1+u3,u2+u3]
    sp=max(S)-min(S)
    if sp<1e-12: continue
    if min(abs(S[p]-S[q]) for p,q in combinations(range(4),2))<1e-12: continue
    K1=2*abs(a1*a2*a3)**(2/3)*sp
    best.append((K1,(V[i],V[j],V[k])))
best.sort()
print("triples",len(best))
res=[]
for K1,seeds in [b for b in best if b[0]<LIM]:
    # exact gcd at large n; actual K
    bl=[Apow(v,nexact) for v in seeds]
    try:
        g,Hp,e,K,q=analyse(words(bl))
    except AssertionError:
        continue
    res.append((float(K[0]),N(g),seeds,K1))
# also triples whose g=1 constant is bigger but could have big gcd: test gcd for all with K1 < 40
res.sort()
for r in res[:15]: print("K=%.5f N(g)=%d seeds=%s K(g=1)=%.4f"%r)
