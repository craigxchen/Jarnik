# Independent search for dominant families of rational curves through e with 4 smooth branches,
# profile = one-sided near-tight 0/1 vertex F (B.C = 15, mult_e = 4). Own formulation:
# u_a(s) = c_a prod_{B in F, a in B} (s - r_B), u_4 has c_4 = 1; branch params s_1..s_4 with s1,s2,s3 = 0,1,-1 fixed (Mobius gauge).
# Equations: c_a A_a(s_i) - D_a(s_i) = 0, A_a = prod_{B in F: a in B, 4 not in B}(s-r_B), D_a = prod_{B in F: 4 in B, a not in B}.
import numpy as np, itertools, sys, json, time
rng=np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 0)
full=frozenset(range(5))
Bs=[frozenset(s) for r in range(1,5) for s in itertools.combinations(range(5),r)]
cuts=[B for B in Bs if 0 in B]
def onesided_classes():
    perms=list(itertools.permutations(range(5)))
    seen={}
    for choice in range(1<<15):
        F=[cuts[i] if (choice>>i)&1 else full-cuts[i] for i in range(15)]
        X=[sum(1 for B in F if a in B) for a in range(5)]
        if len(set(X))!=1: continue
        best=None
        for p in perms:
            img=tuple(sorted(tuple(sorted(p[i] for i in B)) for B in F))
            if best is None or img<best: best=img
        if best not in seen: seen[best]=[frozenset(B) for B in best]
    return list(seen.values())
def setup(F):
    F=list(F)
    A=[[j for j,B in enumerate(F) if a in B and 4 not in B] for a in range(4)]
    D=[[j for j,B in enumerate(F) if 4 in B and a not in B] for a in range(4)]
    return F,A,D
def resid(z,s4,A,D):
    r=z[:15]; c=z[15:19]; S=np.array([0,1,-1,s4])
    out=[]
    for a in range(4):
        for si in S:
            out.append(c[a]*np.prod(si-r[A[a]])-np.prod(si-r[D[a]]))
    return np.array(out)
def jac(z,s4,A,D,h=1e-7):
    f0=resid(z,s4,A,D); J=np.zeros((16,19),complex)
    for k in range(19):
        dz=np.zeros(19,complex); dz[k]=h
        J[:,k]=(resid(z+dz,s4,A,D)-resid(z-dz,s4,A,D))/(2*h)
    return J
def newton(z,s4,A,D,it=200):
    for t in range(it):
        f=resid(z,s4,A,D); nf=np.linalg.norm(f)
        if nf<1e-13: return z,nf,t
        J=jac(z,s4,A,D); step=np.linalg.lstsq(J,-f,rcond=None)[0]
        lam=1.0
        while lam>1e-4:
            zn=z+lam*step
            if np.linalg.norm(resid(zn,s4,A,D))<nf*(1-0.3*lam): break
            lam/=2
        z=zn
    return z,np.linalg.norm(resid(z,s4,A,D)),it
def curve_map(z,s4,F,s):
    r=z[:15]; c=np.append(z[15:19],1.0)
    u=np.array([c[a]*np.prod([s-r[j] for j,B in enumerate(F) if a in B]) for a in range(5)])
    return u[:4]/u[4]
def dominance_rank(z,s4,F,A,D,s):
    # tangent of solution variety in (z,s4) + d/ds ; map to log(u_a/u_4)
    def G(zz,ss4): return resid(zz,ss4,A,D)
    h=1e-7
    Jz=jac(z,s4,A,D); Js4=(G(z,s4+h)-G(z,s4-h))/(2*h)
    Jfull=np.hstack([Jz,Js4[:,None]])          # 16 x 20
    _,sv,Vh=np.linalg.svd(Jfull); K=Vh[16:].conj().T   # kernel 20 x 4 (generic)
    def Phi(zz,ss4,ss): return np.log(curve_map(zz,ss4,F,ss))
    cols=[]
    for k in range(K.shape[1]):
        dz=K[:19,k]*h; ds4=K[19,k]*h
        cols.append((Phi(z+dz,s4+ds4,s)-Phi(z-dz,s4-ds4,s))/(2*h))
    cols.append((Phi(z,s4,s+h)-Phi(z,s4,s-h))/(2*h))
    M=np.array(cols).T
    return np.linalg.svd(M,compute_uv=False), sv
classes=onesided_classes()
print("classes:",len(classes)); sys.stdout.flush()
results=[]
for ci,F in enumerate(classes):
    F,A,D=setup(F)
    found=None; t0=time.time(); tries=0
    while time.time()-t0<float(sys.argv[2]) if len(sys.argv)>2 else 60:
        tries+=1
        s4=complex(*rng.normal(size=2))*1.5
        z=np.concatenate([rng.normal(size=15)+1j*rng.normal(size=15), rng.normal(size=4)+1j*rng.normal(size=4)])*1.5
        z,nf,it=newton(z,s4,A,D)
        if nf<1e-11:
            r=z[:15]; S=np.array([0,1,-1,s4])
            sep=min(min(abs(r[i]-r[j]) for i in range(15) for j in range(i)), min(abs(si-rj) for si in S for rj in r), min(abs(S[i]-S[j]) for i in range(4) for j in range(i)))
            cdeg=min(abs(z[15:19])); c1=min(abs(z[15:19]-1))
            if sep>1e-3 and cdeg>1e-3 and c1>1e-3:
                sdom,sv=dominance_rank(z,s4,F,A,D,complex(0.37,0.21))
                found=dict(tries=tries,sep=float(sep),domsv=[float(x) for x in sdom],kersv=[float(x) for x in sv[-3:]])
                break
    label=sorted(sorted(B) for B in F if len(B)<=2)
    print(ci,"class sizes",sorted(len(B) for B in F),"found:",found); sys.stdout.flush()
    results.append((ci,[sorted(B) for B in F],found))
json.dump(results,open('k_cover_out_%s.json'%(sys.argv[1] if len(sys.argv)>1 else '0'),'w'))
