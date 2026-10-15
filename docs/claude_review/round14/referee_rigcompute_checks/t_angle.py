import math
from reflib import *
from rverify import per_v_angle_ok
for P in ([5,13,17,29],[53,61,89,173],[97,101,109,113]):
    pis=[gprime(p) for p in P]; bad=0; tot=0
    for v in monos(4,6):
        z=Zv(pis,v); D=Dvals(z); N=norm_v(P,v)
        # genuine demand over all odd q (largest possible): odd part of |D_s|
        ms=[]
        for s in range(4):
            x=abs(D[s])
            while x%2==0: x//=2
            ms.append(x)
        a=[(1 if s%2==0 else 2**-0.5)*ms[s]/math.sqrt(N) for s in range(4)]
        tot+=1; bad+= not per_v_angle_ok(a)
    print(P,"genuine per-v angle failures",bad,"of",tot)
