# for Part V family: check minimality of L, compute delta_ij = deg(z_i-z_j)-deg gcd(z_i,z_j), and averaging identity
import itertools,sys
from f2lib import *
def fam(r,L,T=None):
    m=2**r; d=r+1; Q=2**d
    Y=((1<<L) if T is None else T,1)
    masks=[a for a in range(Q) if a&1]
    H={a:[v for v in range(Q) if bin(a&v).count('1')%2==0] for a in masks}
    z={a:prod4([add4(Y,(v,0)) for v in H[a] if v]) for a in masks}
    return masks,z
for r in (2,3,4):
    d=r+1; Lst=(d-2)*2**(d-1)
    for L in (Lst-1,Lst):
        masks,z=fam(r,L); D=deg4(z[masks[0]])
        ok=True; deltas=[]; gd=[]
        for a,b in itertools.combinations(masks,2):
            df=add4(z[a],z[b]); dd=deg4(df)
            if 2*dd>D: ok=False
            g=gcd4(z[a],z[b]); deltas.append(dd-deg4(g)); gd.append(deg4(g))
        m=len(masks)
        print(f"m={m} L={L} D={D} cluster={ok} min/avg/max delta={min(deltas)}/{sum(deltas)/len(deltas):.1f}/{max(deltas)} "
              f"bound D/(2(m-1))={D/(2*(m-1)):.1f}  sum delta={sum(deltas)} <= mD/4={m*D/4}")
