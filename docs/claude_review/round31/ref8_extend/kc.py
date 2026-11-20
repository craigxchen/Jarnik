# K_C: smallest sum of inert primes in a set S with sum log l/(l+1) > log C - (1/2)log 2  (sharpened)
# versus Astra's condition sum log l/(l+1) > log C.  Exact choice by greedy over increasing l is optimal? we
# compute the optimum by DP over small primes (knapsack, minimise sum l) to be safe.
from math import log, sqrt
from sympy import primerange
inert=[l for l in primerange(3,4000) if l%4==3]
def best(target):
    if target<0: return 0,[]
    # greedy prefix is not nec. optimal; do exact search: minimise sum l s.t. sum g>target. g=log l/(l+1) decreasing,
    # and l increasing, so an optimal set can be taken to be a prefix? check by small DP on rounded values.
    s=0;S=[]
    for l in inert:
        s+=log(l)/(l+1);S.append(l)
        if s>target: return sum(S),S
    return None,None
for C in [1,1.5,2,3,4,5,6,8,8.25,10,12,16,20]:
    a,Sa=best(log(C)); b,Sb=best(log(C)-0.5*log(2))
    print(f"C={C:6}: Astra K_C={a} (|S|={len(Sa)})   sharpened K_C={b} (|S|={len(Sb)}) S={Sb[:8]}{'...' if len(Sb)>8 else ''}")
# prefix optimality check by brute force over subsets of first 14 inert primes
import itertools
P=inert[:14]
for C in [2,3,4,5]:
    t=log(C)-0.5*log(2); bestv=None
    for r in range(1,15):
        for T in itertools.combinations(P,r):
            if sum(log(l)/(l+1) for l in T)>t:
                v=sum(T)
                if bestv is None or v<bestv[0]: bestv=(v,T)
    print("C",C,"brute optimum",bestv,"prefix",best(t)[0])
