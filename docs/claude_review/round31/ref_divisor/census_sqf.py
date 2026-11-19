# Squarefree parents on first n split primes; no upper cutoff; primitive candidates with 256 t^4 m <= d.
import sys
from math import isqrt
from itertools import combinations
from collections import defaultdict, Counter
n=int(sys.argv[1])
P=[]; q=5
def isp(x): return x>1 and all(x%k for k in range(2,isqrt(x)+1))
while len(P)<n:
    if q%4==1 and isp(q): P.append(q)
    q+=1
print("primes",P[0],"..",P[-1])
full=(1<<n)-1
prod=[1]*(1<<n)
for mask in range(1,1<<n):
    low=mask&-mask; i=low.bit_length()-1
    prod[mask]=prod[mask^low]*P[i]
cands=[]  # (mask,d,s,t)
for mask in range(1,1<<n):
    d=prod[mask]; s=isqrt(d); r=d-s*s
    if r==0: continue
    t=isqrt(r)
    if t*t!=r: continue
    if 256*t**4>d: continue
    cands.append((mask,d,s,t))
print("candidate divisors",len(cands))
# for each candidate, enumerate complementary squarefree m (disjoint masks) with m<=d/(256t^4)
par=defaultdict(list)
for mask,d,s,t in cands:
    lim=d//(256*t**4)
    comp=[i for i in range(n) if not mask>>i&1]
    # DFS over subsets of comp with product<=lim
    stack=[(0,1,0)]
    while stack:
        k,pr,mm=stack.pop()
        par[mask|mm].append((mask,d,s,t))
        for j in range(k,len(comp)):
            p=P[comp[j]]
            if pr*p>lim: break
            stack.append((j+1,pr*p,mm|(1<<comp[j])))
cnt=Counter(len(v) for v in par.values())
print("parents with >=1 candidate",len(par),"distribution",sorted(cnt.items()))
inc=0; triples=set(); good=0; bad=[]
def gmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def conj(a): return (a[0],-a[1])
for N,v in par.items():
    for x,y,z in combinations(v,3):
        if x[0]^y[0]^z[0]==0:
            inc+=1; key=tuple(sorted((x[1],y[1],z[1])))
            if key in triples: continue
            triples.add(key)
            A=[(w[2],w[3]) for w in (x,y,z)]
            ok=False
            for k in range(3):  # psi_k = sum of others: Im(A_i A_j conj(A_k))=0 and Re>0
                i,j=[l for l in range(3) if l!=k]
                pr=gmul(gmul(A[i],A[j]),conj(A[k]))
                if pr[1]==0 and pr[0]>0: ok=True
            good+=ok
            if not ok: bad.append(key)
print("triangle incidences",inc,"distinct triples",len(triples),"angle relation holds",good,"fails",len(bad))
# 4-cliques with six distinct labels (rows 0,c1,c2,c3)
cl=0; cl3=0
for N,v in par.items():
    S={w[0] for w in v}
    L=sorted(S)
    for a,b,c in combinations(L,3):
        if (a^b in S) and (a^c in S) and (b^c in S):
            labs={a,b,c,a^b,a^c,b^c}
            if len(labs)==6: cl+=1
            else: cl3+=1
print("4-cliques six distinct labels",cl,"4-cliques with repeated labels (incidences)",cl3)
