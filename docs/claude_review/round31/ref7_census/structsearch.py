import sys, itertools, collections, time
from parentmax import parent_max_fast
from witness import pmul,pconj
def dec(c):
    p=[]
    while c: c,d=divmod(c,9); p.append((d%3,d//3))
    return p
def parents(D):
    S=set()
    for line in open(f"c5_D{D}.txt"):
        _,n,*bl=line.split(); N=[(1,0)]
        for b in bl: x=dec(int(b)); N=pmul(N,pmul(x,pconj(x)))
        S.add(tuple(c[0] for c in N))
    return sorted(S)
def p3mul(a,b):
    r=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b): r[i+j]=(r[i+j]+x*y)%3
    return r
def compose(N,phi):  # N(phi(t)), low-first lists
    r=[0]; 
    for c in reversed(N):
        r=p3mul(r,phi); r[0]=(r[0]+c)%3
    while r and r[-1]==0: r.pop()
    return r
res=collections.Counter(); best=0; t0=time.time()
P={D:parents(D) for D in (12,14,16,18)}
print({D:len(v) for D,v in P.items()})
subs={'t^2':[0,0,1],'t^2+1':[1,0,1],'t^2+2':[2,0,1],'t^3+t':[0,1,0,1],'t^3+2t':[0,2,0,1]}
for D,plist in P.items():
    for name,phi in subs.items():
        L=len(phi)-1; DD=D*L
        if DD>40: continue
        for N in plist:
            M=compose(list(N),phi)
            r=parent_max_fast(M,DD)
            key=(D,name,DD, r[0] if r else 'not-sqf')
            res[key]+=1
            if r and r[0]>best: best=r[0]
print("substitutions:"); [print("  ",k,v) for k,v in sorted(res.items(),key=str)]
# products of two coprime 5-parents
res2=collections.Counter()
for D1,D2 in [(12,12),(12,14),(12,16),(14,14),(12,18),(14,16)]:
    pairs=list(itertools.combinations_with_replacement(P[D1],2)) if D1==D2 else list(itertools.product(P[D1],P[D2]))
    for A,B in pairs[:400]:
        M=p3mul(list(A),list(B)); r=parent_max_fast(M,D1+D2)
        res2[(D1,D2, r[0] if r else 'not-sqf')]+=1
print("products:"); [print("  ",k,v) for k,v in sorted(res2.items(),key=str)]
print("time",time.time()-t0)
