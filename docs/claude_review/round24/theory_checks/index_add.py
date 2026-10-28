"""Exact sanity check of Prop 4.2 (index additivity at the rational target (e,e)) on X x X.
For V_{m1} (x) V_{m2} and integer weights (w1,w2), the max over F != 0 of min{w1 j1 + w2 j2 : F_J != 0}
should be w1*t(m1) + w2*t(m2). Certification: full rank mod p of the conditions
'F_{(j1,j2)} = 0 for all w1 j1 + w2 j2 <= K' (random evaluations of each bihomogeneous part) proves index <= K;
the lower bound is the explicit tensor f1 (x) f2."""
import sys, random
import numpy as np
sys.path.insert(0,'.')
from isotypic import S, rank_mod, P
def check(m1,m2,w1,w2,seed=3):
    random.seed(seed)
    S1,S2=S(m1),S(m2); t1,t2=10*m1//3,10*m2//3
    target=w1*t1+w2*t2
    cols=[(a,b) for a in S1 for b in S2]; n=len(cols)
    def rows_upto(K):
        rows=[]
        for j1 in range(0,K//w1+1):
            for j2 in range(0,(K-w1*j1)//w2+1):
                nm=(j1+3)*(j1+2)*(j1+1)//6*(j2+3)*(j2+2)*(j2+1)//6
                for _ in range(min(nm,n)):
                    y=[random.randrange(P) for _ in range(4)]; z=[random.randrange(P) for _ in range(4)]
                    la=[pow(sum(a[i]*y[i] for i in range(4))%P,j1,P) for a in S1]
                    lb=[pow(sum(b[i]*z[i] for i in range(4))%P,j2,P) for b in S2]
                    rows.append([(la[i]*lb[k])%P for i in range(len(S1)) for k in range(len(S2))])
        return np.array(rows,dtype=np.int64)
    A=rows_upto(target)  # conditions for index > target impossible: weighted degree <= target all vanish
    r=rank_mod(A.T) if A.shape[0]>A.shape[1] else rank_mod(A)
    return n,target,r
for (m1,m2,w1,w2) in [(1,1,1,1),(1,1,1,2),(1,2,1,1)]:
    n,target,r=check(m1,m2,w1,w2)
    print(f"V_{m1}(x)V_{m2}, weights ({w1},{w2}): dim {n}, rank of conditions up to weighted degree {target} = {r} "
          f"-> {'max index = '+str(target)+' (certified upper; lower by tensor of extremal sections)' if r==n else 'NOT full rank'}",flush=True)
