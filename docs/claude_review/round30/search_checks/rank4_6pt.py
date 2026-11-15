# D=10 (and 14), c=0: realise the 6-point rank-4 codes found by the MILP (one multiplicity solution per S,
# plus re-solves with random objective to get several multiplicity solutions).
import sys, itertools, numpy as np, pickle
sys.path.insert(0,'.')
from codes4 import functionals
from basesolve import res_base, degen
from scipy.optimize import least_squares, milp, LinearConstraint, Bounds
D=int(sys.argv[1]); L=D-1-D//2; h=(L+1)//2; r=4; size=6
vecs=[v for v in itertools.product([0,1],repeat=r) if any(v)]
rng=np.random.default_rng(3); codes=set()
for sub in itertools.combinations(vecs,size-1):
    S=[(0,)*r]+list(sub)
    if np.linalg.matrix_rank(np.array(S))<r: continue
    cols=functionals(S,r); n=len(S); pairs=list(itertools.combinations(range(n),2)); nv=2*len(cols)
    A=[[1]*nv];lb=[D];ub=[D]
    for (x,y) in pairs:
        dist=[];pos=[];neg=[]
        for c in cols:
            dd=c[x]-c[y]
            for s in (1,-1): dist.append(abs(dd)); pos.append(1 if s*dd==1 else 0); neg.append(1 if s*dd==-1 else 0)
        A+= [dist,pos,neg]; lb+=[L+1,h,h]; ub+=[np.inf]*3
    for trial in range(6):
        res=milp(c=rng.normal(size=nv),constraints=LinearConstraint(np.array(A,float),lb,ub),integrality=np.ones(nv),bounds=Bounds(0,D))
        if res.status!=0: break
        m=np.round(res.x).astype(int); colsl=[];sg=[]
        for i in range(nv):
            colsl+=[tuple(cols[i//2])]*m[i]; sg+=[1 if i%2==0 else -1]*m[i]
        key=tuple(sorted(zip(colsl,sg))); codes.add(key)
print("D",D,"distinct (labelled) 6-point rank-4 codes sampled:",len(codes),flush=True)
real=0; out=[]
for key in codes:
    C=np.array([k[0] for k in key]).T; sg=np.array([k[1] for k in key],float)  # 6 x D sets, base signs
    Fs=C[1:].astype(float)
    for s in range(40):
        z0=np.concatenate([rng.normal(0,1,D), sg*np.abs(rng.normal(0.8,0.4,D))])
        meth='lm' if Fs.shape[0]*L+2>=2*D else 'trf'
        rr=least_squares(res_base,z0,args=(Fs,L),method=meth,xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=4000)
        q=rr.x[:D]+1j*rr.x[D:]
        if np.linalg.norm(rr.fun)<1e-10 and degen(q)>1e-3: real+=1; out.append((Fs,q)); break
print("realised",real,"of",len(codes),flush=True)
pickle.dump(out,open(f"r4_6pt_D{D}.pkl","wb"))
