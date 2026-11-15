# Rank-5 extension hunt: on a known rank-4 family (given base sets), add one more set F and solve the
# overdetermined system from starts on the family. Reports nondegenerate zeros.
import numpy as np, itertools, sys, time
sys.path.insert(0,'.'); from basesolve import *
from scipy.optimize import least_squares
def run(D,L,Fs_base,seeds_q,label,maxF=None):
    Fs_base=np.array(Fs_base); span=np.linalg.matrix_rank(Fs_base)
    hits=[]; t0=time.time(); n=0
    for F in itertools.product([0,1],repeat=D):
        F=np.array(F)
        if F.sum()<L+1 or F.sum()>D-L-1: continue
        if np.linalg.matrix_rank(np.vstack([Fs_base,F]))==span: continue
        if F[0]==1 and False: continue
        n+=1
        Fs=np.vstack([Fs_base,F]).astype(float); bestc=1e9; bestq=None
        for q0 in seeds_q:
            z0=np.concatenate([q0.real,q0.imag])
            r=least_squares(res_base,z0,args=(Fs,L),method='lm',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=3000)
            q=r.x[:D]+1j*r.x[D:]; c=np.linalg.norm(r.fun)
            if c<bestc: bestc=c; bestq=q
        if bestc<1e-10 and degen(bestq)>1e-4:
            hits.append((F,bestq)); print(label,"HIT",F,degen(bestq),np.round(bestq,5),flush=True)
    print(label,"tested",n,"extra sets; hits",len(hits),"time %.0f"%(time.time()-t0),flush=True)
    return hits
if __name__=="__main__":
    which=sys.argv[1]
    if which=="hc12":
        D=12;L=5;U=list(itertools.product([0,1],repeat=3)); cols=[]
        for u in U: cols+=[u]*((12+4*(-1)**sum(u))//8)
        cols=np.array(cols); Fs=[[1]*D]+[list(1-cols[:,i]) for i in range(3)]
        out=solve_base(np.array([[1]*D]+[list(1-cols[:,i]) for i in range(3)]+[list(cols[:,i]) for i in range(3)]),L,starts=300,seed=11)
        seeds=[o[2] for o in out if o[0]<1e-10 and o[1]>1e-3][:4]
        print("seeds",len(seeds),flush=True)
        run(D,L,Fs,seeds,"hc12")
