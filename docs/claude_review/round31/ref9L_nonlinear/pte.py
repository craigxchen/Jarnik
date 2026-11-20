# Own regression on weighted PTE families (ref8L_attack Prop 2.1) and maximal-modulus families.
from itertools import product, combinations
from math import gcd, prod
from sympy import factorint
def lcm(a,b): return a//gcd(a,b)*b
def pp_list(q): return [p**e for p,e in factorint(q).items()]
def sigs(marks,pps):
    S=[]
    for y in marks:
        s=0
        for k,m in enumerate(pps):
            cl=sorted({z%m for z in marks}); assert len(cl)<=2,(m,cl)
            s|=cl.index(y%m)<<k
        S.append(s)
    return S
def flats(marks,pps):
    S=sigs(marks,pps); assert len(set(S))==len(S)
    idx={s:i for i,s in enumerate(S)}
    by={}
    for i,j in combinations(range(len(S)),2): by.setdefault(S[i]^S[j],[]).append((i,j))
    planes=set()
    for h,prs in by.items():
        for (i,j),(k,l) in combinations(prs,2):
            P=frozenset((i,j,k,l))
            if len(P)==4: planes.add(P)
    cubes=set()
    for P in planes:
        x=S[min(P)]
        for v in range(len(S)):
            if v in P: continue
            sh=x^S[v]
            Q=[idx.get(S[i]^sh) for i in P]
            if None not in Q: cubes.add(P|frozenset(Q))
    return S,planes,cubes
def report(label,marks,q):
    pps=pp_list(q)
    S,planes,cubes=flats(marks,pps)
    L=max(marks)-min(marks)
    diam=lambda F:max(marks[i] for i in F)-min(marks[i] for i in F)
    okp=all((2*diam(P))**2>=q for P in planes); okc=all((4*diam(C))**3>=q*q for C in cubes)
    okpg=(not planes) or (2*L)**2>=q; okcg=(not cubes) or (4*L)**3>=q*q
    import math
    print(f"{label}: B={len(marks)} log q/log L={math.log(q)/math.log(L):.3f} planes={len(planes)} 3-flats={len(cubes)} "
          f"global ineq plane/cube={okpg}/{okcg} per-flat ineq plane/cube={okp}/{okc}")
def pte(label,nodes,r,T):
    # nodes: list of (alpha,beta) forms alpha*T-beta, weight alpha^2; moments mu_0=sum eps*alpha, mu_1=sum eps*beta
    ell=[a*T-b for a,b in nodes]; tot=prod(ell)
    marks=[]
    for eps in product((0,1),repeat=len(nodes)):
        if sum(e*a for e,(a,b) in zip(eps,nodes)): continue
        if r>=2 and sum(e*b for e,(a,b) in zip(eps,nodes)): continue
        marks.append(sum(e*a*a*(tot//l) for e,(a,b),l in zip(eps,nodes,ell)))
    q=1
    for l in ell: q=lcm(q,abs(l))
    report(label,marks,q)
def q2(marks):
    B=len(marks); q=1
    for fl in product((0,1),repeat=B-1):
        bits=(0,)+fl; g=0
        for i,j in combinations(range(B),2):
            if bits[i]==bits[j]: g=gcd(g,marks[i]-marks[j])
        if g: q=lcm(q,g)
    return q
T=1009
for k in (4,5,8):
    a=[1]*((k+1)//2)+[-1]*(k//2)
    pte(f"(k={k},r=1) T={T}",[(v,v*j) for j,v in enumerate(a)],1,T)
vec=[(-3,-3),(-3,-2),(0,-1),(1,2),(1,3),(2,0),(2,1)]
pte(f"(7,2) T={T}",[(x+7*y,y) for x,y in vec],2,T)
a,b=1009,431
fams={"parallelogram":[0,a,b,a+b],
 "5 projective":[0,a*(a+b),b*(a+2*b),a*(a+b)+b*(a+2*b),(a+b)*(a+2*b)],
 "6 projective":[0,(a-b)*(a-2*b),a*(a-2*b),-2*b*(a-b),-a*b,a*a-4*a*b+2*b*b]}
nd=(1009,1013,1019)
fams["cube"]=[sum((prod(nd)//v)*e for v,e in zip(nd,bits)) for bits in product((0,1),repeat=3)]
for lab,m in fams.items(): report(lab+" (maximal q_2)",m,q2(m))
# robustness: another specialization
for T in (10007,):
    for k in (4,5):
        aa=[1]*((k+1)//2)+[-1]*(k//2)
        pte(f"(k={k},r=1) T={T}",[(v,v*j) for j,v in enumerate(aa)],1,T)
    pte(f"(7,2) T={T}",[(x+7*y,y) for x,y in vec],2,T)
