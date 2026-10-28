# Referee test: unit pencils f = lam*A/C with lam in Q^*, lam != +-1 (not covered by theory.md's
# enumeration, which uses lam = s = +-1 only). For each (A,C) with deg <= Dmax, compute all rational
# "special" lam (where B = C - lam*A can be non-reduced, gain a pair factor, or vanish to higher
# order on some L_T), and evaluate the exact gain there.
import itertools, sys, time, collections, random
from fractions import Fraction
from flint import fmpz_mpoly_ctx, fmpq_poly, fmpz_poly
import sympy as sp
Dmax=int(sys.argv[1])
ctx=fmpz_mpoly_ctx.get(('x',8),'lex'); G=ctx.gens(); X=G[:7]; E=G[7]
edges=[(i,j) for i in range(7) for j in range(i+1,7)]
lin={e:X[e[0]]-X[e[1]] for e in edges}
lin_sig={}
for e in edges: lin_sig[str(lin[e])]=e; lin_sig[str(-lin[e])]=e
Ts=[frozenset(c) for r in range(3,7) for c in itertools.combinations(range(7),r)]
subs={}
for T in Ts:
    i0=min(T); subs[T]=[X[i] if i not in T or i==i0 else X[i0]+E*X[i] for i in range(7)]+[E]
def lowest(p):
    q=p.compose(*subs_cur)
def init_form(p,T):
    q=p.compose(*subs[T]); mu=min(m[7] for m in q.monoms())
    return mu,{m:c for m,c in zip(q.monoms(),q.coeffs()) if m[7]==mu}
def nuT(p,T):
    q=p.compose(*subs[T]); return min(m[7] for m in q.monoms())
def canon(ms):
    verts=sorted(set(v for e in ms for v in e)); k=len(verts); best=None
    for perm in itertools.permutations(range(k)):
        mp={verts[i]:perm[i] for i in range(k)}
        t=tuple(sorted(tuple(sorted((mp[a],mp[b]))) for (a,b) in ms))
        if best is None or t<best: best=t
    return best
def poly(ms):
    p=ctx.from_dict({(0,)*8:1})
    for e in ms: p=p*lin[e]
    return p
def nu_unit(ms,T): return sum(1 for (a,b) in ms if a in T and b in T)
def prop(d1,d2):
    # return Fraction lam with d2 = lam*d1, or None
    if set(d1)!=set(d2) or not d1: return None
    k=next(iter(d1)); lam=Fraction(int(d2[k]),int(d1[k]))
    for m in d1:
        if Fraction(int(d2[m]),int(d1[m]))!=lam: return None
    return lam
random.seed(7)
lines=[([random.randint(-50,50) for _ in range(7)],[random.randint(-50,50) for _ in range(7)]) for _ in range(2)]
t,L=sp.symbols('t L')
def restrict_line(ms,ln):
    p0,q0=ln; r=sp.Integer(1)
    for (a,b) in ms: r*= (p0[a]-p0[b]) + t*(q0[a]-q0[b])
    return sp.expand(r)
def gain(A,C,pp,qq):
    pa=poly(A); pc=poly(C); B=qq*pc-pp*pa
    if B.is_zero(): return None
    c,facs=B.factor(); pairB=collections.Counter(); inter=[]
    for f,m in facs:
        if f.total_degree()==1 and str(f) in lin_sig: pairB[lin_sig[str(f)]]+=m
        elif f.total_degree()>=1: inter.append((f,m))
    nuY={(i,T):nuT(f,T) for i,(f,m) in enumerate(inter) for T in Ts}
    D=len(A); Lb=31*D; Fb=0
    for T in Ts:
        a=nu_unit(A,T); cc=nu_unit(C,T); mu=min(a,cc); Lb-=mu
        b=sum(m for (u,v),m in pairB.items() if u in T and v in T)+sum(m*nuY[(i,T)] for i,(f,m) in enumerate(inter))
        if max(a,b,cc)>mu: Fb+=1
    for e in edges:
        if e in A or e in C or e in pairB: Fb+=1
    for i,(f,m) in enumerate(inter): Fb+=31*f.total_degree()-sum(nuY[(i,T)] for T in Ts)
    return Lb-Fb
gbest=[-10**9]; tot=0; hist=collections.Counter(); best=None; t0=time.time(); ncand=0
for D in range(1,Dmax+1):
    reps=sorted(set(canon(ms) for ms in itertools.combinations_with_replacement(edges,D)))
    for A in reps:
        others=[e for e in edges if e not in set(A)]
        pa=poly(A)
        inA={T:init_form(pa,T) for T in Ts}
        for C in itertools.combinations_with_replacement(others,D):
            pc=poly(C); cand=set()
            gen=-sum(1 for T in Ts if nu_unit(A,T)!=nu_unit(C,T))-sum(1 for e in edges if e in A or e in C)
            hist[(D,'generic',gen)]+=1
            if gen>gbest[0]: gbest[0]=gen
            for T in Ts:
                ma,da=inA[T]; mc,dc=init_form(pc,T)
                if ma==mc:
                    lam=prop(da,dc)
                    if lam is not None: cand.add(lam)
            for (i,j) in edges:
                if (i,j) in A or (i,j) in C: continue
                s=[X[k] if k!=j else X[i] for k in range(7)]+[E]
                ra=pa.compose(*s); rc=pc.compose(*s)
                lam=prop(dict(zip(ra.monoms(),ra.coeffs())),dict(zip(rc.monoms(),rc.coeffs())))
                if lam is not None: cand.add(lam)
            if D>=2:
                discs=[]
                for ln in lines:
                    a_=restrict_line(A,ln); c_=restrict_line(C,ln)
                    discs.append(sp.Poly(sp.discriminant(c_-L*a_,t),L))
                g=sp.gcd(discs[0],discs[1])
                if g.degree()>0:
                    for r in sp.Poly(g,L).ground_roots().keys() if False else sp.roots(g,filter='Q').keys():
                        if r.is_rational and r!=0: cand.add(Fraction(int(sp.numer(r)),int(sp.denom(r))))
            cand.discard(Fraction(1)); cand.discard(Fraction(-1)); cand.discard(Fraction(0))
            for lam in cand:
                g_=gain(A,C,lam.numerator,lam.denominator); ncand+=1
                if g_ is None: continue
                hist[(D,g_)]+=1
                if best is None or g_>best[0]: best=(g_,A,C,lam)
            tot+=1
    print('D',D,'(A,C) pairs so far',tot,'special lam evaluated',ncand,'time',round(time.time()-t0,1)); sys.stdout.flush()
    print('  special-lam hist top:',sorted([(k[1],c) for k,c in hist.items() if k[0]==D and len(k)==2],reverse=True)[:8]); print('  generic gain max so far',gbest[0]); sys.stdout.flush()
print('best',best)
