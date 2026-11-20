# log L(x)/x, L(x)=lcm_{l<=x, l=3 mod 4 prime}(l+1); exact via max prime-power exponents
import sys, math
from sympy import primerange, factorint
N=int(sys.argv[1]) if len(sys.argv)>1 else 10**6
# smallest prime factor sieve up to N+1
M=N+2
spf=list(range(M))
for i in range(2,int(M**0.5)+1):
    if spf[i]==i:
        for j in range(i*i,M,i):
            if spf[j]==j: spf[j]=i
isprime=lambda n: n>1 and spf[n]==n
best={}  # prime r -> max exponent
logL=0.0; theta3=0.0
checkpoints=[10**k for k in range(2,7)]+[3*10**k for k in range(2,6)]
checkpoints=sorted(c for c in checkpoints if c<=N)
ci=0
# split by size of prime r relative to x: track contributions from large primes r > x/B at checkpoints later
out=[]
for l in range(3,N+1,4):
    if not isprime(l): 
        pass
    else:
        theta3+=math.log(l)
        n=l+1
        while n>1:
            r=spf[n]; e=0
            while n%r==0: n//=r; e+=1
            if e>best.get(r,0):
                logL+=(e-best.get(r,0))*math.log(r); best[r]=e
    while ci<len(checkpoints) and l+4>checkpoints[ci]:
        x=checkpoints[ci]
        # contribution of primes r > x/B for B=10,100
        big10=sum(e*math.log(r) for r,e in best.items() if r> x/10)
        big100=sum(e*math.log(r) for r,e in best.items() if r> x/100)
        out.append((x,logL,logL/x,theta3/x,big10/x,big100/x))
        ci+=1
print("x  logL  logL/x  theta3/x  [r>x/10]/x  [r>x/100]/x")
for row in out: print("%8d %12.1f %.4f %.4f %.4f %.4f"%row)
