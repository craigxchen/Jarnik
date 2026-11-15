# Rank-5 hunt: for every 5-point type (rank 4) at given D, realise it (several solutions), then for every extra
# base-form set F independent of the span, solve the augmented system (type + F) from (a) the realised solutions,
# (b) random starts. Any rank-5 cluster contains such a (5-point rank-4) + 1 configuration.
import sys, pickle, numpy as np, itertools, time
sys.path.insert(0,'.')
from basesolve import res_base, degen
from satlat import balance_violation
from scipy.optimize import least_squares
D,c=int(sys.argv[1]),int(sys.argv[2]); L=D-1-(D+2*c)//2; h=(L+1)//2
T=pickle.load(open(sys.argv[3],'rb')); nrand=int(sys.argv[4]) if len(sys.argv)>4 else 3
rng=np.random.default_rng(5)
def lm(Fs,z0):
    meth='lm' if Fs.shape[0]*L+2>=2*D else 'trf'
    r=least_squares(res_base,z0,args=(Fs,L),method=meth,xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=3000)
    q=r.x[:D]+1j*r.x[D:]; return np.linalg.norm(r.fun),q
t0=time.time(); hits=0; realised=0
for idx,t in enumerate(T):
    E=np.array(t).T   # M x D normalised signs
    Fs=((E[1:]!=E[0])).astype(float)   # base-form sets relative to point 0
    y0=E[0].astype(float)              # sign of Im of base roots: base root = rho if eps=+1 else conj
    sols=[]
    for s in range(30):
        z0=np.concatenate([rng.normal(0,1,D), y0*np.abs(rng.normal(0.8,0.4,D))])
        cost,q=lm(Fs,z0)
        if cost<1e-10 and degen(q)>1e-3: sols.append(q)
        if len(sols)>=3: break
    if not sols: print("type",idx,"not realised in 30 starts",flush=True); continue
    realised+=1
    span=np.linalg.matrix_rank(Fs)
    for F in itertools.product([0,1],repeat=D):
        F=np.array(F,float)
        if F.sum()<L+1: continue
        A=np.vstack([Fs,F])
        if np.linalg.matrix_rank(A)==span: continue
        if span<4: continue   # rank-3 types cannot reach rank 5 with one more set
        enew=np.where(F>0,-E[0],E[0])
        okb=True
        for e in E:
            p=np.sum((enew==1)&(e==-1)); n=np.sum((enew==-1)&(e==1))
            if p<h or n<h or p+n<L+1: okb=False; break
        if not okb: continue
        if balance_violation(((np.vstack([E,enew])[1:]-E[0])//2).astype(int),L,h,1) is not None: continue
        starts=[np.concatenate([q.real,q.imag]) for q in sols]+[np.concatenate([rng.normal(0,1,D), y0*np.abs(rng.normal(0.8,0.4,D))]) for _ in range(nrand)]
        for z0 in starts:
            cost,q=lm(A,z0)
            if cost<1e-10 and degen(q)>1e-4:
                hits+=1; print("HIT type",idx,"F",F.astype(int),"degen",degen(q),np.round(q,6).tolist(),flush=True); break
print("D",D,"c",c,"types",len(T),"realised",realised,"hits",hits,"time %.0f"%(time.time()-t0),flush=True)
