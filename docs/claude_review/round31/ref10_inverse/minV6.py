# Exact exhaustive column minima: owners enumerated by increasing l1 weight, exact integer solvability of the
# non-owner system via Smith normal form of D_N (no bounds on non-owner entries).
import sys, itertools, numpy as np
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_decomp
from minV import setup
def patterns(no, par, w):
    # all integer x with x_o = par_o mod 2 and sum|x| == w
    base=sum(par)
    if w<base or (w-base)%2: return
    extra=(w-base)//2  # number of +2 increments in absolute value, distributed
    for dist in itertools.combinations_with_replacement(range(no),extra):
        mag=list(par)
        for o in dist: mag[o]+=2
        nz=[o for o in range(no) if mag[o]]
        for signs in itertools.product((1,-1),repeat=len(nz)):
            x=[0]*no
            for o,s in zip(nz,signs): x[o]=s*mag[o]
            yield x
for q in map(int,sys.argv[1:]):
    I,S,i0,rows,D,Hm,Bm=setup(q); n=len(I)
    own=[t for t,k in enumerate(I) if k in S]; non=[t for t,k in enumerate(I) if k not in S]
    DN=Matrix(D[:,non].tolist())
    Sm,Pm,Qm=smith_normal_decomp(DN,domain=ZZ); assert Sm==Pm*DN*Qm
    dd=[int(Sm[i,i]) for i in range(min(Sm.shape))]; rank=sum(1 for x in dd if x)
    PD=Pm*Matrix(D.tolist())  # transformed system: (P D)[:,own] x + Sm Q^{-1} y = P b
    PDo=[[int(PD[i,j]) for j in own] for i in range(PD.rows)]
    Bf={i0:np.zeros(n,dtype=int)}
    for r_,i in enumerate(rows): Bf[i]=Bm[r_]
    Pb=Pm*Matrix(Bm.tolist())
    def ok(x,t):
        for i in range(PD.rows):
            val=int(Pb[i,t])-sum(PDo[i][o]*x[o] for o in range(len(own)))
            d=dd[i] if i<len(dd) else 0
            if d==0:
                if val!=0: return False
            elif val % d: return False
        return True
    no=len(own); tot=0; partot=0; dist={}; rowsum=[0]*no; maxe=0
    for t in range(n):
        pars=[[ (c+Bf[I[o]][t])%2 for o in own] for c in (0,1)]
        pm=min(sum(p) for p in pars); partot+=pm
        found=None; w=pm
        while found is None:
            for p in pars:
                for x in patterns(no,p,w):
                    if ok(x,t): found=x; break
                if found: break
            if found is None: w+=1
        tot+=w; dist[w-pm]=dist.get(w-pm,0)+1
        for o in range(no): rowsum[o]+=abs(found[o])
        maxe=max(maxe,max(abs(v) for v in found))
    print(dict(q=q,m=no,inv_DN=sorted(set(abs(x) for x in dd)),parity_bound=partot,exact_min=tot,
               excess_per_column=dist,owner_rows_at_one_optimum=sorted(rowsum),row_bound=(q-3)/2,max_owner_entry=maxe,
               CS_bound=round(0.5*(no*(q-1)-(no*q*(q-1))**0.5-no),2)),flush=True)
