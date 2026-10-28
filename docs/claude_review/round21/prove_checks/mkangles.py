# emit split primes <= B, a Gaussian prime a+bi (a>b>0), and theta_p as u64 with 2^64 = pi/2
import sys
from mpmath import mp, atan2, pi, floor, mpf
from sympy import primerange
mp.dps=50
B=int(sys.argv[1])
out=[]
for p in primerange(5,B+1):
    if p%4!=1: continue
    a=next(a for a in range(1,p) if int((p-a*a)**0.5+0.5)**2==p-a*a)
    b=int((p-a*a)**0.5+0.5)
    if a<b: a,b=b,a
    th=atan2(b,a)
    u=int(floor(th/(pi/2)*mpf(2)**64+mpf('0.5')))
    out.append((p,a,b,u))
print(len(out))
for p,a,b,u in out: print(p,a,b,u)
