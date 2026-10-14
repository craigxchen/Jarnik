import itertools, sys
from rlib import *
def Hf(X,b,cap):
    hs=[];e=0
    while True:
        hs.append(h_func_exact(X,e,cap)); e+=1
        if hs[-1]==len(X): break
    return [hs[0]]+[hs[i]-hs[i-1] for i in range(1,len(hs))]
for M,b,j in [(3,3,2),(3,4,2),(3,3,3),(3,4,3),(2,3,2),(2,4,2),(4,4,2),(4,3,3)]:
    cap=[2*j+1]*b
    G=list(itertools.product(range(-j,j+1),repeat=b))
    P=[y for y in G if sum(y)%M==0]; Q=[y for y in G if sum(y)%M!=0]
    HP=Hf(P,b,cap); HQ=Hf(Q,b,cap)
    print("M",M,"(b,j)",(b,j),"rho(P)",len(HP)-1,"rho(Q)",len(HQ)-1,"alpha",2*j*b-(len(HQ)-1)); sys.stdout.flush()
