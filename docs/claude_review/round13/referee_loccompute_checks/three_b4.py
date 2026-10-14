# For each non-width three-cut block (author's list of 126), test J(A) = J_w(R=2) + pairs of {-1,0,1} directions
# in all degrees <= rho+1. J(A) exact (fmpz); L mod p (lower bound) => equality exact.
import numpy as np, itertools, math, time, json, sys
from flint import nmod_mat
from reng import REng, PR
from rlib import J_exact, monos
regions=[(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1),(1,1,1)]
recs=[json.loads(l) for l in open('../loccompute_checks/three_mixed_out.txt')]
recs=sorted([r for r in recs if r['b0_mixed']==4],key=lambda r:(r['d'],r['rho']))
lo,hi=int(sys.argv[1]),int(sys.argv[2])
def run(rec):
    sz=rec['regions']; Ls=rec['L']
    cols=[[],[],[]]
    for rg,z in zip(regions,sz):
        for _ in range(z):
            for j in range(3): cols[j].append(rg[j])
    d=len(cols[0]); COLS=np.array(cols,dtype=np.int64)
    K=np.array(list(itertools.product([-1,0,1],repeat=d)),dtype=np.int64)
    A=K[np.all(np.abs(K@COLS.T)<=np.array(Ls),axis=1)]
    Dm=rec['rho']+1
    JA=J_exact([tuple(int(x) for x in p) for p in A],Dm,cap=2)
    dirs=[v for v in itertools.product(range(-2,3),repeat=d) if any(v) and math.gcd(*[abs(x) for x in v])==1 and next(x for x in v if x)>0]
    W={v:len(np.unique(A@np.array(v))) for v in dirs}
    d1=[v for v in dirs if max(abs(x) for x in v)<=1 and W[v]<=Dm+2]
    # pair generators (exact), computed once
    gens=[]
    for v in dirs:
        if W[v]<=Dm: gens.append((list(v),[1],[(W[v],)],W[v],1))
    npairs=0; cache={}
    for U in itertools.combinations(d1,2):
        Ua=np.array(U)
        img=np.unique(A@Ua.T,axis=0)
        if W[U[0]]*W[U[1]]==len(img): continue
        key=tuple(map(tuple,img.tolist()))
        if key not in cache: cache[key]=J_exact([tuple(int(x) for x in p) for p in img],Dm)
        J=cache[key]
        npairs+=1
        for D in range(Dm+1):
            for vec in J[D]: gens.append((list(map(list,U)),vec,monos(2,D),D,2))
    res=[]
    for DT in range(Dm+1):
        if len(JA[DT])==0: res.append((DT,0,0)); continue
        E=REng(d,DT)
        rows=[]
        for g in gens:
            U,vec,tm,D,b=g
            if D>DT: continue
            if b==1: f=E.pullback(vec,tm,[U])
            else: f=E.pullback(vec,tm,U)
            if f.any(): rows+=E.ideal_rows(f,D,DT)
        M=np.array(rows); r=nmod_mat(M.shape[0],M.shape[1],[int(x) for x in M.flatten()],PR).rank() if rows else 0
        res.append((DT,r,len(JA[DT])))
    ok=all(a==b for _,a,b in res)
    return ok,res,npairs
for i,rec in enumerate(recs[lo:hi]):
    t0=time.time()
    ok,res,npairs=run(rec)
    print(json.dumps(dict(i=lo+i,regions=rec['regions'],L=rec['L'],rho=rec['rho'],b0_col=rec['b0_col'],b0_mixed=rec['b0_mixed'],
        pairs_suffice=ok,deficits=[(D,b-a) for D,a,b in res if a!=b],nongrid_pairs=npairs,t=round(time.time()-t0,1))),flush=True)
