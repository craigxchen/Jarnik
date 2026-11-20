import math,sys
def isprime(n): return n>1 and all(n%d for d in range(2,math.isqrt(n)+1))
def two_sq(p):
    for a in range(1,math.isqrt(p)+1):
        b=math.isqrt(p-a*a)
        if a*a+b*b==p: return (a,b)
def gm(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
for P in [30,50,100]:
    primes=[p for p in range(P,2*P+1) if isprime(p) and p%4==1]; q=math.prod(primes)
    pts=[(1,0)]
    for p in primes:
        a,b=two_sq(p); pts=[gm(z,f) for z in pts for f in ((a,b),(a,-b))]
    allp=[gm(u,z) for z in pts for u in [(1,0),(0,1),(-1,0),(0,-1)]]
    allp.sort(key=lambda z: math.atan2(z[1],z[0]))
    n=len(allp); best=1
    for i in range(n):
        k=1
        while k<n and ((allp[i][0]-allp[(i+k)%n][0])**2+(allp[i][1]-allp[(i+k)%n][1])**2)**2<=8**4*q: k+=1
        best=max(best,k)
    print(f"P={P} genuine max points within chord 8 q^(1/4) of the first (all units) = {best}")
