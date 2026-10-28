# count n <= X with all prime factors = 1 mod 4 and <= B, and sum of d(n)=prod(e+1)
import sys
from sympy import primerange
def count(X,B):
    ps=[p for p in primerange(5,B+1) if p%4==1]
    cnt=0; sd=0
    def dfs(i,n,d):
        nonlocal cnt,sd
        cnt+=1; sd+=d
        for j in range(i,len(ps)):
            p=ps[j]; m=n*p
            if m>X: break
            e=1
            while m<=X:
                dfs(j+1,m,d*(e+1)); m*=p; e+=1
    dfs(0,1,1)
    return cnt,sd
for X in [10**12,10**14,10**16]:
    for B in [100,200]:
        print(X,B,count(X,B),flush=True)
