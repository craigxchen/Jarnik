import sys,itertools
from collections import Counter
def pm(a,b):
    r=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): r[i+j]=(r[i+j]+x*y)%3
    return r
for D in range(1,int(sys.argv[1])+1):
    K=D//2; C=Counter()
    for a in itertools.product(range(3),repeat=D):
        A=list(a)+[1]
        for B in itertools.product(range(3),repeat=D):
            B=list(B)+[0]
            N=tuple((x+y)%3 for x,y in zip(pm(A,A),pm(B,B)))
            C[(N,tuple(A[K+1:]),tuple(B[K+1:]))]+=1
    M=max(C.values()); print(D,M,sum(1 for v in C.values() if v==M))
