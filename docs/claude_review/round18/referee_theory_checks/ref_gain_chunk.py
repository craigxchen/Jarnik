# Independent recomputation of the abc-gain of unit pencils (theory.md Prop 2.6 / §2.4(2)).
# Exact: Q-factorisation by flint; nu_T of interior factors by exact symbolic substitution
# x_i -> x_{i0} + e*x_i (i in T\{i0}) and lowest power of e (no specialisation).
import itertools, sys, time, collections
from flint import fmpz_mpoly_ctx
D=int(sys.argv[1]); OI=int(sys.argv[2]); C0=int(sys.argv[3]); C1=int(sys.argv[4])
ctx=fmpz_mpoly_ctx.get(('x',8),'lex')  # x0..x6 = labels 1..7, x7 = e
G=ctx.gens(); X=G[:7]; E=G[7]
edges=[(i,j) for i in range(7) for j in range(i+1,7)]
lin={e:X[e[0]]-X[e[1]] for e in edges}
Ts=[frozenset(c) for r in range(3,7) for c in itertools.combinations(range(7),r)]
pairsT=[frozenset(e) for e in edges]
subs={}
for T in Ts:
    i0=min(T); s=[X[i] if i not in T or i==i0 else X[i0]+E*X[i] for i in range(7)]+[E]
    subs[T]=s
def order_e(p):
    return min(m[7] for m in p.monoms())
def nuT(p,T):
    q=p.compose(*subs[T])
    return order_e(q)
def canon(ms):
    verts=sorted(set(v for e in ms for v in e)); k=len(verts); best=None
    for perm in itertools.permutations(range(k)):
        mp={verts[i]:perm[i] for i in range(k)}
        t=tuple(sorted(tuple(sorted((mp[a],mp[b]))) for (a,b) in ms))
        if best is None or t<best: best=t
    return best
allms=list(itertools.combinations_with_replacement(edges,D))
import os, json as _j
_f='reps_D%d.json'%D
if os.path.exists(_f): reps=[tuple(tuple(e) for e in r) for r in _j.load(open(_f))]
else:
    reps=sorted(set(canon(ms) for ms in allms)); _j.dump(reps,open(_f,'w'))

def poly(ms):
    p=ctx.from_dict({(0,)*8:1})
    for e in ms: p=p*lin[e]
    return p
def nu_unit(ms,T): return sum(1 for (a,b) in ms if a in T and b in T)
lin_sig={}
for e in edges:
    lin_sig[str(lin[e])]=e; lin_sig[str(-lin[e])]=e
cnt=0; hist=collections.Counter(); best=[]; t0=time.time()
cacheY={}
for A in reps[OI:OI+1]:
    Aset=set(A); pa=poly(A)
    others=[e for e in edges if e not in Aset]
    for ci,C in enumerate(itertools.combinations_with_replacement(others,D)):
        if ci<C0 or ci>=C1: continue
        pc=poly(C)
        for s in (1,-1):
            B=pc - s*pa
            c,facs=B.factor()
            pairB=collections.Counter(); inter=[]
            for f,m in facs:
                if f.total_degree()==1 and str(f) in lin_sig: pairB[lin_sig[str(f)]]+=m
                else: inter.append((f,m))
            Lb=31*D; Fb=0
            cacheY={}
            for T in Ts:
                a=nu_unit(A,T); cc=nu_unit(C,T); mu=min(a,cc); Lb-=mu
                b=sum(m for (u,v),m in pairB.items() if u in T and v in T)
                for f,m in inter:
                    key=(str(f),T)
                    if key not in cacheY: cacheY[key]=nuT(f,T)
                    b+=m*cacheY[key]
                if max(a,b,cc)>mu: Fb+=1
            for e in edges:
                if e in Aset or e in C or e in pairB: Fb+=1
            for f,m in inter:
                Fb+=31*f.total_degree()-sum(cacheY[(str(f),T)] for T in Ts)
            g=Lb-Fb; cnt+=1; hist[g]+=1
            if not best or g>best[0][0]: best=[(g,A,C,s)]
            elif g==best[0][0] and len(best)<5: best.append((g,A,C,s))
import json
print(json.dumps({'D':D,'orbit':OI,'c0':C0,'c1':C1,'n':cnt,'hist':{str(k):v for k,v in hist.items()},'best':str(best[:2]) if best else ''}))
