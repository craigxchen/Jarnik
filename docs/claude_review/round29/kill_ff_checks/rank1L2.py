# L=2 with common real part: cluster = zero-sum subsets of signed reals y_k with distinct |y_k| (rank-1 W).
import itertools, numpy as np, sys
D=int(sys.argv[1]); B=int(sys.argv[2])
best=0;arg=None
subs=np.array(list(itertools.product([0,1],repeat=D)),dtype=np.int64)
for absv in itertools.combinations(range(1,B+1),D):
    a=np.array(absv)
    for sg in itertools.product([1,-1],repeat=D-1):
        y=a*np.array((1,)+sg)
        cnt=int((subs@y==0).sum())
        if cnt>best: best=cnt;arg=y.copy()
print(D,B,best,arg)
