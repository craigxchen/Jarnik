"""Upper bound for the Poisson intensity of gap events with w(v) > WMAX (valid since WMAX >= 10, where (4/pi)(arcsin x + arcsin(x/sqrt2)) <= 2.18 x):
   (1/2) * 2.18 * ( prod_j (1 + 2/(sqrt p_j - 1)) - 1 - sum_{0<w(v)<=WMAX} e^{-w(v)/2} ),  using Prob(Bad_v) <= 2.18 e^{-w/2}
   ((4/pi)(arcsin x + arcsin(x/sqrt2)) <= 2.18 x for x <= 5^{-1/2})."""
import sys, math
def tail(WMAX, ps):
    s=len(ps); lw=[math.log(p) for p in ps]; acc=0.0
    def dfs(j,w,nz):
        nonlocal acc
        if j==s:
            if nz: acc+=math.exp(-w/2)
            return
        m=0
        while w+m*lw[j]<=WMAX+1e-12:
            dfs(j+1,w+m*lw[j],nz or m>0)*0 if False else None
            mult=1 if m==0 else 2
            for _ in range(mult): dfs(j+1,w+m*lw[j],nz or m>0)
            m+=1
    dfs(0,0.0,False)
    full=math.prod(1+2/(math.sqrt(p)-1) for p in ps)-1
    return 0.5*2.18*(full-acc)
P=[53,61,73,89,97,101,109,113,137,149,157,173]
for s in [3,5,7,9,10,11,12]: print("first",s,"primes >=53, WMAX=22: tail <= %.3f"%tail(22,P[:s]))
print("101..157 WMAX=24: %.3f"%tail(24,[101,109,113,137,149,157]))
print("29..97 WMAX=22: %.3f"%tail(22,[29,37,41,53,61,73,89,97]))
x=5**-0.5; print("check 2.18 constant:", (4/math.pi)*(math.asin(x)+math.asin(x/math.sqrt(2)))/x)
