# T5: moving coefficients tuned to n to create common zeros / high multiplicity
import sys, itertools, random
sys.path.insert(0,'.')
from gcdlib import *
random.seed(7)
def nullvec(M):
    # one nonzero kernel vector mod P of matrix M (rows), random combination of kernel basis
    M=[r[:] for r in M]; rows=len(M); cols=len(M[0]); piv=[]; r=0
    for c in range(cols):
        pr=next((i for i in range(r,rows) if M[i][c]%P),None)
        if pr is None: continue
        M[r],M[pr]=M[pr],M[r]; inv=pow(M[r][c],P-2,P); M[r]=[x*inv%P for x in M[r]]
        for i in range(rows):
            if i!=r and M[i][c]%P:
                f=M[i][c]; M[i]=[(x-f*y)%P for x,y in zip(M[i],M[r])]
        piv.append(c); r+=1
        if r==rows: break
    free=[c for c in range(cols) if c not in piv]
    x=[0]*cols
    for c in free: x[c]=random.randrange(1,P)
    for i,c in enumerate(piv):
        x[c]=(-sum(M[i][f]*x[f] for f in free))%P
    return x
boxpts=list(itertools.product(range(4),repeat=3))
def mult_at_1(a):
    # multiplicity of root t=1 of dense poly a mod P
    a=[int(x) for x in a]; m=0
    while True:
        s=sum(a)%P
        if s: return m
        # divide by (t-1)
        q=[0]*(len(a)-1); acc=0
        for i in range(len(a)-1,0,-1):
            acc=(acc+a[i])%P; q[i-1]=acc
        a=q; m+=1
for L in [200,400,800]:
    n=(L, round(L*0.7071), round(L*0.5773))
    for s in [4,6,8]:
        SF=random.sample(boxpts,s); SG=random.sample(boxpts,s)
        # order s-1 vanishing at t=1: sum c_u (n.u)^j = 0, j=0..s-2
        def tune(S):
            M=[[pow(sum(x*y for x,y in zip(u,n)),j,P) for u in S] for j in range(s-1)]
            c=nullvec(M); return {u:ci for u,ci in zip(S,c)}
        F=tune(SF); G=tune(SG)
        d=common([F,G],n)
        print('L',L,'s',s,'gcd deg',d,'mult at 1 of F,G',mult_at_1(spec(F,n)),mult_at_1(spec(G,n)))
