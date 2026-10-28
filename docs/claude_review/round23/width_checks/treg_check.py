# Direct check: reg(T_n(a,b)) vs closed form F_n(a,b), via ranks mod p=2^61-1.
# full rank at F certifies reg <= F over Q; deficiency at F-1 mod p is evidence of equality (lower bounds over Q are
# not certified by this script).
import itertools, sys
from reglib import rank_at
def T(n,a,b):
    pts=[]
    for x in itertools.product(range(-b,a+1),repeat=n-1):
        y=list(x)+[a-b-sum(x)]
        if sum(v for v in y if v>0)<=a: pts.append(tuple(x))
    return pts
F={3:lambda a,b:a+b+min(a,b),4:lambda a,b:a+b+min((a+b)//2,2*min(a,b)),
   5:lambda a,b:a+b+min(min(a,b)+max(a,b)//3,3*min(a,b))}
n=int(sys.argv[1]); A=int(sys.argv[2])
for a in range(A+1):
    for b in range(a,A+1):
        if a+b==0: continue
        S=T(n,a,b); f=F[n](a,b)
        full=rank_at(S,f)==len(S); defic=rank_at(S,f-1)<len(S)
        print(f"n={n} a={a} b={b} |T|={len(S)} F={f} full_at_F={full} deficient_at_F-1={defic}",flush=True)
