# Independent recomputation of the abc-gain of unit pencils (theory.md Prop 2.6 / §2.4(2)).
# Exact: Q-factorisation by flint; nu_T of interior factors by exact symbolic substitution
# x_i -> x_{i0} + e*x_i (i in T\{i0}) and lowest power of e (no specialisation).
import itertools, sys, time, collections
from flint import fmpz_mpoly_ctx
D=int(sys.argv[1]); EXACT=(sys.argv[2]=='exact') if len(sys.argv)>2 else True
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
reps=sorted(set(canon(ms) for ms in allms))
print('D',D,'A-orbits',len(reps)); sys.stdout.flush()
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
for A in reps:
    Aset=set(A); pa=poly(A)
    others=[e for e in edges if e not in Aset]
    for C in itertools.combinations_with_replacement(others,D):
        pc=poly(C)
        for s in (1,-1):
            B=pc - s*pa
            open("last_B.txt","w").write(str((A,C,s)))
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
            if cnt%1000==0: import resource; print("prog",cnt,resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,round(time.time()-t0,1)); sys.stdout.flush()
            if not best or g>best[0][0]: best=[(g,A,C,s)]
            elif g==best[0][0] and len(best)<5: best.append((g,A,C,s))
    cacheY.clear(); print("A done",A,cnt,round(time.time()-t0,1),sorted(hist.items(),reverse=True)[:3]); sys.stdout.flush()
print('pencils',cnt,'time',round(time.time()-t0,1))
print('gain histogram (top 8):',sorted(hist.items(),reverse=True)[:8])
print('examples of max:',best[:3])
