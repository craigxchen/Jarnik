import sys, time
from minsearch import search
d=5; Bs=[int(x) for x in sys.argv[1].split(',')]; Lmin,Lmax,low=int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])
for B in Bs:
    best=None; t0=time.time()
    for L in range(Lmin,Lmax+1):
        for x in range(1<<low):
            T=(1<<L)^x
            if T%2==0 and (T>>1)%2==0 and False: pass
            r=search(d,B,T,best[0]+1 if best else 10**9)
            if r[0] is not None and (best is None or r[0]<best[0]): best=(r[0],T,r[1])
    print(f"d=5 B={bin(B)} T=t^L+x, L in [{Lmin},{Lmax}], x<2^{low}: least critical D'=",(None if best is None else (best[0],bin(best[1]),best[2])),f"({time.time()-t0:.0f}s)",flush=True)
