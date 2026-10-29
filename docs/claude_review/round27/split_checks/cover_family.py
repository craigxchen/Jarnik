# Covering families at one-sided vertices (numerical evidence).
# For a one-sided near-tight family F (15 sets B, one per cut, X_a equal), consider rational curves
#   u_j(s) = lam_j * prod_{B in F, j in B} (s - s_B)      (profile E_B.C = 1_{B in F}, B.C = 15)
# passing through e = (1:..:1) at 4 distinct parameters s_1..s_4 (mult_e C = 4, E.C~ = 4, -K_Bl.C~ = 3).
# Expected dimension of this family: 4 = dim T.  If it exists and the map (family, s) -> T is dominant, then every
# Laurent polynomial f satisfies 4 ord_e f <= sum_{B in F} pi_B(f): WTB holds at the vertex x_F.
# Unknowns: s_B (15), log-free lam_1..lam_4 (lam_5 = 1), s_1..s_4.  PGL_2 fixed by s_{B0}=0, s_{B1}=1, s_{B2}=-1 (generic).
# Equations: u_j(s_i) - u_5(s_i) = 0 normalised, j=1..4, i=1..4 (16).  Free dims: 23 - 3 - 16 = 4.
import numpy as np, itertools, sys
rng=np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 0)
def onesided_families():
    cuts=[]
    for r in (1,2):
        for c in itertools.combinations(range(5),r):
            B=frozenset(c); cuts.append((B,frozenset(range(5))-B))
    out=[]
    for mask in range(1<<15):
        F=[cuts[i][(mask>>i)&1] for i in range(15)]
        X=[sum(1 for B in F if a in B) for a in range(5)]
        if len(set(X))==1: out.append(F)
    return out
def canon(F):
    return min(tuple(sorted(tuple(sorted(p[a] for a in B)) for B in F)) for p in itertools.permutations(range(5)))
def run(F,tries=40):
    F=list(F); nB=len(F)
    inc=np.array([[1 if j in B else 0 for B in F] for j in range(5)])  # 5 x 15
    fixed={0:0.0,1:1.0,2:-1.0}
    free_idx=[i for i in range(nB) if i not in fixed]
    def unpack(v):
        sB=np.zeros(nB,dtype=complex)
        for i,val in fixed.items(): sB[i]=val
        sB[free_idx]=v[:12]; lam=np.concatenate([v[12:16],[1.0]]); si=v[16:20]
        return sB,lam,si
    def U(sB,lam,s):
        return lam*np.prod(np.where(inc==1,(s-sB)[None,:],1.0),axis=1)
    def eqs(v):
        sB,lam,si=unpack(v)
        out=[]
        for s in si:
            u=U(sB,lam,s); out.extend((u[:4]-u[4])/ (np.abs(u[4])+1e-300))
        return np.array(out)
    def jac(v,h=1e-7):
        f0=eqs(v); J=np.zeros((16,20),dtype=complex)
        for k in range(20):
            e=np.zeros(20,dtype=complex); e[k]=h; J[:,k]=(eqs(v+e)-f0)/h
        return J
    for t in range(tries):
        v=rng.normal(size=20)+1j*rng.normal(size=20)
        ok=False
        with np.errstate(all='ignore'):
            for it in range(200):
                f=eqs(v)
                if not np.all(np.isfinite(f)): break
                J=jac(v); dv=np.linalg.lstsq(J,-f,rcond=None)[0]
                v=v+dv
                if np.linalg.norm(f)<1e-12: ok=True; break
        if not ok: continue
        sB,lam,si=unpack(v)
        # nondegeneracy: distinct s_i, distinct s_B, s_i not roots
        if min(abs(a-b) for a,b in itertools.combinations(si,2))<1e-4: continue
        if min(abs(a-b) for a,b in itertools.combinations(sB,2))<1e-4: continue
        if min(abs(s-b) for s in si for b in sB)<1e-4: continue
        J=jac(v); sv=np.linalg.svd(J,compute_uv=False)
        rankJ=int(np.sum(sv>1e-8*sv[0]))
        # tangent space of the solution variety (kernel of J)
        _,_,Vh=np.linalg.svd(J); K=Vh[rankJ:].conj().T          # 20 x (20-rankJ)
        # dominance: map (v, s) -> y(s) = (u_j/u_5)_{j<4} at a generic s; differential along kernel directions and d/ds
        s0=0.37+0.21j
        def y(v,s):
            sB,lam,_=unpack(v); u=U(sB,lam,s); return u[:4]/u[4]
        h=1e-7; cols=[]
        for k in range(K.shape[1]):
            cols.append((y(v+h*K[:,k],s0)-y(v,s0))/h)
        cols.append((y(v,s0+h)-y(v,s0))/h)
        D=np.array(cols).T; dsv=np.linalg.svd(D,compute_uv=False)
        rankD=int(np.sum(dsv>1e-6*dsv[0]))
        return dict(found=True,rankJ=rankJ,family_dim=20-rankJ,rank_dominance=rankD,dsv=dsv)
    return dict(found=False)
if __name__=="__main__":
    fams=onesided_families(); seen={}
    for F in fams:
        k=canon(F)
        if k not in seen: seen[k]=F
    W=[frozenset(B) for B in [{0,3,4},{0,1,3,4},{3,4},{1,3,4},{0,2,4},{0,1,2,4},{2,4},{1,2,4},{0,2,3},{0,1,2,3},{2,3},{1,2,3},{1},{0},{0,1}]]
    kW=canon(W)
    for k,F in seen.items():
        res=run(F)
        tag="W-class" if k==kW else ""
        sizes=sorted(len(B) for B in F)
        if res['found']:
            print(sizes,tag,"family dim",res['family_dim'],"dominance rank",res['rank_dominance'],"(4 = dominant)", "sv",np.round(res['dsv'],6),flush=True)
        else:
            print(sizes,tag,"no solution found",flush=True)
