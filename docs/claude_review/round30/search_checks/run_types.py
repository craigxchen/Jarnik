import sys, pickle, numpy as np
sys.path.insert(0,'.')
from ffsolve import *
D,M,c=map(int,sys.argv[1:4]); starts=int(sys.argv[4]) if len(sys.argv)>4 else 200
L=LofD(D,c)
T=pickle.load(open(f"types_D{D}_M{M}_c{c}.pkl","rb"))
for idx,t in enumerate(T):
    E=np.array(t).T
    best=solve(E,L,starts=starts,seed=idx)
    rk=rank_of(E)
    msg="type %d rank %d cost %.2e degen %.2e"%(idx,rk,best[0],best[2])
    if best[0]<1e-10 and best[2]>1e-4:
        C=cluster_through(best[1],L,E[0]); msg+=" REALISED; full cluster %d rank %d rho=%s"%(len(C),rank_of(C),np.round(best[1],6))
    print(msg, flush=True)
    if len(T)<=10: print(E)
