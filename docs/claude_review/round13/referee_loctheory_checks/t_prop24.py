# Prop 2.4: P = {y in [-2,2]^5 : sum y = 0 mod 3}; H_P, H_Q exact; alpha, rho.
import itertools, sys
from rlib import *
b=5; j=2
G=list(itertools.product(range(-j,j+1),repeat=b))
Pset=[y for y in G if sum(y)%3==0]; Qset=[y for y in G if sum(y)%3!=0]
print("|P|",len(Pset),"|Q|",len(Qset))
cap=[2*j+1]*b
def Hfun(X,emax):
    hs=[]
    for e in range(emax+1):
        hs.append(h_func_exact(X,e,cap)); 
        if hs[-1]==len(X): break
    H=[hs[0]]+[hs[i]-hs[i-1] for i in range(1,len(hs))]
    return H
HP=Hfun(Pset,20); print("H_P",HP,"rho(P)",len(HP)-1); sys.stdout.flush()
HQ=Hfun(Qset,20); print("H_Q",HQ,"rho(Q)",len(HQ)-1)
HG=[len(monos_eq(b,e,cap)) for e in range(21)]
print("H_Gamma",HG)
s=4*b
ok=all(HG[e]==(HP[e] if e<len(HP) else 0)+(HQ[s-e] if 0<=s-e<len(HQ) else 0) for e in range(s+1))
print("linkage identity",ok,"alpha=",s-(len(HQ)-1))
