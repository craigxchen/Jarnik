# Deterministic scan: all split-only N <= 10^E with prod(e+1) >= T, primes <= PMAX; keep windows with C <= CM[M].
import sys, math, pickle
from dfs_norms import gen
from enum_clusters import points
CM={3:1.0,4:2.5,5:4.5,6:6.5,7:8.0,8:9.5,9:11.0,10:12.5,11:14,12:15,13:16,14:17}
E,T,PMAX,part,nparts=map(int,sys.argv[1:6])
L=gen(10**E,T,PMAX)[part::nparts]
out=[]
for fac in L:
    ps,es,pts=points(fac); N=math.prod(p**e for p,e in fac.items())
    sR=N**0.25; n=len(pts); ext=pts+[(t+2*math.pi,z,a) for t,z,a in pts]; res=[]
    for s in range(n):
        for M in range(3,15):
            if s+M-1>=len(ext): break
            C=sR*(ext[s+M-1][0]-ext[s][0])
            if C>CM[M]: 
                if C>17: break
                continue
            res.append((M,C,[ext[s+k][1] for k in range(M)],[ext[s+k][2] for k in range(M)]))
    if res: out.append((ps,es,N,res))
pickle.dump(out,open(f"scan_{E}_{part}.pkl","wb"))
print(len(L),"norms",len(out),"with windows",sum(len(r[3]) for r in out))
