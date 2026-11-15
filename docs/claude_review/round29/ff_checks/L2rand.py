import numpy as np, sys, itertools
n=int(sys.argv[1]); B=int(sys.argv[2]); T=int(sys.argv[3]); rs=np.random.RandomState(int(sys.argv[4]) if len(sys.argv)>4 else 0)
masks=((np.arange(2**n)[:,None]>>np.arange(n))&1)
best=0
for it in range(T):
    V=np.stack([rs.randint(1,B+1,size=n),rs.randint(-B,B+1,size=n)],1)
    if len(set(map(tuple,V)))<n: continue
    s=masks@V; h=s[:,0]*1000003+s[:,1]
    m=np.bincount(h-h.min()).max()
    if m>best: best=m; print(best,V.tolist(),flush=True)
print("final",best)
