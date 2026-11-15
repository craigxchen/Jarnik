import numpy as np, sys
sys.path.insert(0,'/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round30/search_checks')
from ffsolve import *
# cube: base Q roots r_k, -r_k (k=0..3); sets as in Thm C. Build normalised types from the exact example
y=np.array([-9/4,-2/3,1/4,8/3]); eps=np.array([1,-1,-1,1]); r=eps/y+1j*y
roots=np.concatenate([r,-r]); L=3
sgn=np.sign(roots.imag); rho=np.where(sgn>0,roots,np.conj(roots))
base=sgn.astype(int)  # Q uses root 'roots' => normalised sign = sign Im
C=cluster_through(rho,L,base); print("cluster size",len(C),"rank",rank_of(C))
best=solve(C,L,starts=40)
print("refound: cost",best[0],"degeneracy",best[2])
C2=cluster_through(best[1],L,C[0]); print("cluster size at new solution",len(C2),"rank",rank_of(C2))
