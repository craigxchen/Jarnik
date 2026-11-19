import sys
from minsearch import search
from f4lib import gcd2
d=4; Tmax=int(sys.argv[2]); Bs=[int(x) for x in sys.argv[1].split(',')]
for B in Bs:
    best=None; cnt=0
    for T in range(Tmax):
        c=T
        for b in (1,2,4,8): c=gcd2(c,b)
        if c!=1: continue
        r=search(d,B,T,best[0]+1 if best else 10**9)
        if r[0] is not None and (best is None or r[0]<best[0]): best=(r[0],T,r[1])
    print(f"d=4 B={bin(B)} all T<2^{Tmax.bit_length()-1}: least critical D'=", (best if best is None else (best[0], bin(best[1]), best[2])), flush=True)
