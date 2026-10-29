# reg(n Q_F(1)) for one-sided families F (W-class and the smallest-volume class); compare reg/n with 15/4.
import sys, time; sys.path.insert(0,'/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/split_checks')
from reglib import *
W=[{0,3,4},{0,1,3,4},{3,4},{1,3,4},{0,2,4},{0,1,2,4},{2,4},{1,2,4},{0,2,3},{0,1,2,3},{2,3},{1,2,3},{1},{0},{0,1}]
F1=[{0}, {1}, {2}, {3}, {4}] + [set(c) for c in itertools.combinations(range(5),3)]  # sizes 1^5 3^10
for name,F in [("W",W),("1^5 3^10",F1)]:
    Fs=set(frozenset(B) for B in F)
    for nn in range(1,5):
        t=time.time()
        p=[nn if B in Fs else None for B in Bs]
        S=lattice_points(p)
        r=reg(S)
        pi=pis(S); cost=sum(pi[idx[B]] for B in Fs)
        print(name,"n=",nn,"|S|=",len(S),"reg=",r,"cost_F=",cost,"4reg/cost=",round(4*r/cost,4),"reg/n=",round(r/nn,3),f"({time.time()-t:.1f}s)",flush=True)
