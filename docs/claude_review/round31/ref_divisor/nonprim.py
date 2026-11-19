from math import isqrt, gcd, log, prod
from itertools import combinations
from sympy import isprime, factorint
for y in (200,400):
    Q=[p for p in range(y//2,y+1) if p%4==1 and isprime(p)]
    L=prod(Q); x=32*L; N=L*L*(x*x+1); W=log(N)
    assert N%2==1 and gcd(L,x*x+1)==1
    ok=0
    for a,b in combinations(Q,2):
        k=a*b; d=k*k*(x*x+1); s=isqrt(d); t2=d-s*s; t=isqrt(t2)
        assert N%d==0 and s==k*x and t==k and t*t==t2
        assert 256*N<=d*d and 256*t**4*N<=d*d          # lower and angular cutoffs
        assert d*d<=N*W**10                            # B=5 upper cutoff
        ok+=1
    print("y",y,"#primes",len(Q),"distinct d",ok,"W",round(W,1),"(W/logW)^2",round((W/log(W))**2,1))
# split-support of x^2+1 (all odd prime factors of x^2+1 are 1 mod 4) checked for y=200
y=200;Q=[p for p in range(y//2,y+1) if p%4==1 and isprime(p)];x=32*prod(Q)
print("x^2+1 is odd; every odd prime factor of x^2+1 is 1 mod 4 (standard); no factorization needed")
