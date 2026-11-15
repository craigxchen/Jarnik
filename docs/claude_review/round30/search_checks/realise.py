import sys, pickle, numpy as np
sys.path.insert(0,'.')
from basesolve import res_base, degen
from scipy.optimize import least_squares
D,M,c,starts=map(int,sys.argv[1:5]); L=D-1-(D+2*c)//2
fn=sys.argv[5] if len(sys.argv)>5 else f"types_D{D}_M{M}_c{c}.pkl"
T=pickle.load(open(fn,'rb')); rng=np.random.default_rng(1); out=[]
for idx,t in enumerate(T):
    E=np.array(t).T; Fs=(E[1:]!=E[0]).astype(float); best=(1e9,0,None)
    for s in range(starts):
        z0=np.concatenate([rng.normal(0,1,D), E[0]*np.abs(rng.normal(0.8,0.4,D))])
        meth='lm' if Fs.shape[0]*L+2>=2*D else 'trf'
        r=least_squares(res_base,z0,args=(Fs,L),method=meth,xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=4000)
        q=r.x[:D]+1j*r.x[D:]; cst=np.linalg.norm(r.fun); dg=degen(q)
        if cst<1e-10 and dg>best[1]: best=(cst,dg,q)
        if cst<1e-10 and dg>1e-2: break
    ok=best[2] is not None and best[1]>1e-4
    out.append((idx,ok,best))
    print("type",idx,"rank",np.linalg.matrix_rank(Fs),"realised" if ok else "not realised","degen %.2e"%best[1],flush=True)
pickle.dump(out,open(f"real_D{D}_M{M}_c{c}.pkl","wb"))
print("SUMMARY D",D,"M",M,"c",c,"types",len(T),"realised",sum(o[1] for o in out))
