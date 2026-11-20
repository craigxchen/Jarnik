# Sanity checks of the X.1 mechanism (not part of any proof).
from math import gcd, isqrt
from cut import labels, gprime, mul
from sympy import factorint, n_order, primerange
import random
# order lemma: ord(u^n)=4, n=2^k  =>  ord(u)=4n, in cyclic groups Z/N
bad=0
for N in range(1,400):
    for k in range(0,6):
        n=2**k
        for a in range(N):
            d=N//gcd(a,N)            # ord of a in Z/N
            dn=d//gcd(d,n)           # ord of n*a
            if dn==4 and d!=4*n: bad+=1
print("order lemma counterexamples:",bad)
random.seed(1)
for (P,five) in [(100,True),(100,False),(500,False)]:
    S=labels(P,five); T=max(((l-1)&-(l-1)).bit_length()-1 for l in S); n=2**(T-1)
    viol=0; tot=0; maxratio=0
    for p in S:
        a,b0=gprime(p)
        for _ in range(300):
            g=(random.randint(-30,30),random.randint(-30,30))
            Ng=g[0]**2+g[1]**2
            if Ng==0: continue
            pi=(a,b0) if random.random()<.5 else (a,-b0)
            num=mul(pi,mul(g,g)); d=gcd(gcd(*num),Ng); al=(num[0]//d,num[1]//d); b=Ng//d
            J=al[0]**2+al[1]**2; assert J==p*b*b and J%2==1 and gcd(*al)==1
            z=(1,0)
            for _ in range(n): z=mul(z,al)
            x,y=z
            for D in (x-y,x+y):
                tot+=1
                assert D%2!=0 and D!=0 and gcd(D,J)==1
                assert D*D<=2*J**n
                if any(D%l==0 for l in S): viol+=1
            assert (x-y)**2+(x+y)**2==2*J**n
    print(P,"5mod8" if five else "all split","T=",T,"n=",n,"tests",tot,"parent-prime divisors found:",viol)
