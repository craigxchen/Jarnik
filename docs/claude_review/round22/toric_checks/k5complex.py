# Complexified d=1 system at k=5: unknowns r_T, r'_T (independent), equations p_n(R_j) = p_n(R'_j), n=1..7, j=1..4,
# R_j = {r_T : T∋j}, R'_j = {r'_T : T∋j}. Gauge (complex affine): r_top = 0, r'_top = 1. 28 eqs / 28 complex unknowns.
# Real solutions of the original problem are those with r' = conj(r) after a real-affine regauge.
# Newton from random complex starts; record nondegenerate solutions (c_j = (P_j - P'_j)/2i nonzero and distinct,
# r_T != r'_T, all 30 roots distinct).
import numpy as np, itertools, sys
Ts=[T for r in range(1,5) for T in itertools.combinations(range(4),r)]
iTop=Ts.index((0,1,2,3)); idx=[i for i in range(15) if i!=iTop]
Mem=np.array([[1.0 if j in T else 0.0 for T in Ts] for j in range(4)])
def unpack(z):
    r=np.zeros(15,complex); rp=np.zeros(15,complex); r[idx]=z[:14]; rp[idx]=z[14:]; r[iTop]=0; rp[iTop]=1
    return r,rp
def F(z):
    r,rp=unpack(z)
    return np.concatenate([Mem@(r**n)-Mem@(rp**n) for n in range(1,8)])
def J(z):
    r,rp=unpack(z); rows=[]
    for n in range(1,8):
        d=n*r**(n-1); dp=n*rp**(n-1)
        for j in range(4):
            rows.append(np.concatenate([(Mem[j]*d)[idx], -(Mem[j]*dp)[idx]]))
    return np.array(rows)
def newton(z,it=200):
    for _ in range(it):
        f=F(z)
        if np.max(np.abs(f))<1e-13: return z,True
        try: dz=np.linalg.solve(J(z),-f)
        except np.linalg.LinAlgError: return z,False
        if np.max(np.abs(dz))>10: dz*=10/np.max(np.abs(dz))
        z=z+dz
    return z,np.max(np.abs(F(z)))<1e-10
def cvals(r,rp):
    return np.array([(np.prod(-r[Mem[j]>0])-np.prod(-rp[Mem[j]>0]))/2j for j in range(4)])
if __name__=="__main__":
    ntry=int(sys.argv[1]); rng=np.random.default_rng(int(sys.argv[2])); sols=[]
    conv=0
    for t in range(ntry):
        z0=(rng.normal(size=28)+1j*rng.normal(size=28))*rng.choice([0.5,1,2])
        z,ok=newton(z0)
        if not ok: continue
        conv+=1; r,rp=unpack(z); c=cvals(r,rp)
        allr=np.concatenate([r,rp]); dmin=min(abs(allr[a]-allr[b]) for a in range(30) for b in range(a+1,30))
        nd=min(np.min(np.abs(c)), min(abs(c[a]-c[b]) for a in range(4) for b in range(a+1,4)), dmin)
        sv=np.linalg.svd(J(z),compute_uv=False)
        if nd>1e-6:
            sols.append(z); print(f"try {t}: nondeg nd={nd:.2e} cond={sv[-1]/sv[0]:.1e} c={np.round(c,4)}",flush=True)
    print("converged",conv,"nondegenerate",len(sols))
    np.save(f"k5complex_{sys.argv[2]}.npy",np.array(sols) if sols else np.zeros((0,28)))
