# reg(n Q_F(1)) for the W-class one-sided family, larger n (mod-p rank; mod-p reg >= true reg).
import sys, time; sys.path.insert(0,'/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/split_checks')
from reglib import *
W=[{0,3,4},{0,1,3,4},{3,4},{1,3,4},{0,2,4},{0,1,2,4},{2,4},{1,2,4},{0,2,3},{0,1,2,3},{2,3},{1,2,3},{1},{0},{0,1}]
Fs=set(frozenset(B) for B in W)
def reg_from(S,d0,p=2**31-1):
    S=list(S);N=len(S)
    for d in range(d0,80):
        ms=monos(d)
        if len(ms)<N: continue
        t=time.time()
        rows=[[pow(int(m[0])%p,a[0],p)*pow(int(m[1])%p,a[1],p)%p*pow(int(m[2])%p,a[2],p)%p*pow(int(m[3])%p,a[3],p)%p for a in ms] for m in S]
        rk=nmod_mat(rows,p).rank()
        print("  d=",d,"rank",rk,"of",N,f"({time.time()-t:.0f}s)",flush=True)
        if rk==N: return d
for nn in [5,6]:
    p=[nn if B in Fs else None for B in Bs]
    S=lattice_points(p)
    print("n=",nn,"|S|=",len(S),flush=True)
    r=reg_from(S,int(3.4*nn))
    print("W n=",nn,"reg=",r,"reg/n=",r/nn," vs 15/4",flush=True)
