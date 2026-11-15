# purely-imaginary-root clusters: zero-sum F iff sum_F y^m = 0 for odd m <= L ; |y| distinct
import itertools, numpy as np, sys
D=int(sys.argv[1]); L=int(sys.argv[2]); B=int(sys.argv[3])
odd=[m for m in range(1,L+1) if m%2]
subs=np.array(list(itertools.product([0,1],repeat=D)),dtype=np.int64)
signs=np.array([(1,)+s for s in itertools.product([1,-1],repeat=D-1)],dtype=np.int64)
best=0;arg=None
for absv in itertools.combinations(range(1,B+1),D):
    a=np.array(absv,dtype=np.int64)
    Y=signs*a  # S x D
    ok=np.ones((len(subs),len(signs)),dtype=bool)
    for m in odd:
        ok&=(subs@(Y**m).T==0)
    cnt=ok.sum(axis=0); i=int(cnt.argmax())
    if cnt[i]>best: best=int(cnt[i]); arg=Y[i].copy(); print(best,arg,flush=True)
print("FINAL",D,L,B,best,arg)
