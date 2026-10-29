import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/split_checks')
from reglib import *
# sanity: n R_5 (all p_B = n): reg should be floor(10n/3)
for nn in range(1,5):
    S=lattice_points([nn]*30); print("nR5 n=",nn,"|S|=",len(S),"reg=",reg(S),"floor(10n/3)=",(10*nn)//3)
# pair segment
S=[(k,-k,0,0) for k in range(0,4)]; print("segment reg",reg(S), "pis",pis(S))
