# random types: M sign vectors with pairwise balance/distance, saturated-balance filter, generic rank M-1;
# then LM realisation. Used to test the dimension-count heuristic (c=1) and to hunt rank>=5 (c=0).
import sys, numpy as np, itertools, time
sys.path.insert(0,'.')
from basesolve import res_base, degen
from satlat import balance_violation
from scipy.optimize import least_squares
from scipy.optimize import linprog
D,c,M,ntypes,nstarts=map(int,sys.argv[1:6]); seed=int(sys.argv[6]) if len(sys.argv)>6 else 0
L=D-1-(D+2*c)//2; h=(L+1)//2; rng=np.random.default_rng(seed)
def ok_pair(x,y):
    p=np.sum((x==1)&(y==-1)); n=np.sum((x==-1)&(y==1)); return p>=h and n>=h and p+n>=L+1
def gordan(E):
    A=(E[1:]-E[0])/2; r=linprog(np.zeros(D),A_eq=A,b_eq=np.zeros(len(A)),bounds=[(1,None)]*D,method='highs'); return r.status==0
found=0; tried=0; t0=time.time(); real=0
while tried<ntypes and time.time()-t0<3000:
    E=[rng.choice([1,-1],D)]
    for it in range(2000):
        x=rng.choice([1,-1],D)
        if all(ok_pair(x,y) for y in E): E.append(x)
        if len(E)==M: break
    if len(E)<M: continue
    E=np.array(E)
    if np.linalg.matrix_rank(E[1:]-E[0])<M-1 or not gordan(E): continue
    if balance_violation((E[1:]-E[0])//2,L,h,1) is not None: continue
    tried+=1
    Fs=(E[1:]!=E[0]).astype(float)
    best=None
    for s in range(nstarts):
        z0=np.concatenate([rng.normal(0,1,D), E[0]*np.abs(rng.normal(0.8,0.4,D))])
        meth='lm' if Fs.shape[0]*L+2>=2*D else 'trf'
        r=least_squares(res_base,z0,args=(Fs,L),method=meth,xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=3000)
        q=r.x[:D]+1j*r.x[D:]; cost=np.linalg.norm(r.fun)
        if cost<1e-10 and degen(q)>1e-3: best=q; break
    if best is not None:
        real+=1
        np.save(f"rt_D{D}_c{c}_M{M}_{real}.npy",np.vstack([E[0],Fs.astype(int)]).astype(float)) 
        np.save(f"rtq_D{D}_c{c}_M{M}_{real}.npy",best)
print("D",D,"c",c,"M",M,"rank",M-1,"L",L,"expected dim",2*D-2-(M-1)*L,"types tried",tried,"realised",real,"time %.0f"%(time.time()-t0),flush=True)
