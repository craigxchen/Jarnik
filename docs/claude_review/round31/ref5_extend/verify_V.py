import itertools,sys
from f2lib import *
def run(r):
    m=2**r; d=r+1; Q=2**d; L=(d-2)*2**(d-1)
    V=list(range(Q))  # polys of degree<d
    Y=(1<<L,1)       # t^L + theta
    # functionals l(v)=parity(mask&v), need l(1)=1 -> mask odd
    masks=[a for a in range(Q) if a&1]
    assert len(masks)==m
    H={a:[v for v in V if bin(a&v).count('1')%2==0] for a in masks}
    z={a:prod4([add4(Y,(v,0)) for v in H[a] if v]) for a in masks}
    N=None
    for a in masks:
        n=norm4(z[a]); N=n if N is None else N
        assert n==N
    Nd=prod4([add4(Y,(v,0)) for v in V if v not in (0,1)])
    assert Nd[1]==0 and Nd[0]==N
    D=deg4(z[masks[0]]); assert all(deg4(z[a])==D for a in masks) and D==(m-1)*L
    mx=0; tot_delta=0; mindelta=10**9
    for a,b in itertools.combinations(masks,2):
        df=add4(z[a],z[b]); dd=deg4(df); mx=max(mx,dd)
        # gcd of z_a,z_b over F_4[t]: compute via norms? use structure: check formula
        J=[v for v in H[a] if v in set(H[b])]
        PJ1=1
        for v in J: PJ1=mul(PJ1,1^v)
        PJY=prod4([add4(Y,(v,0)) for v in J if v])
        rhs=m4((PJ1,0),PJY)
        assert rhs==df
        assert dd<=D//2 if D%2==0 else dd<=D/2
    print(f"m={m} L={L} D={D} parent={2*D} maxdiff={mx} D/2={D/2} ok")
for r in range(2,int(sys.argv[1])+1): run(r)
