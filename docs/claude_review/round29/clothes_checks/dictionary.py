# exact checks of the dictionaries (a), (d) and of the twin cut inequality
import math, sympy, random, itertools
from math import gcd, isqrt
from common import gpoints
random.seed(7)
def mul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def conj(a): return (a[0],-a[1])
split=[p for p in sympy.primerange(5,300) if p%4==1]
# ---- (a): m=|z-z1|^2/2, m(2N-m)=Im(z conj z1)^2, g=gcd(m,2N), m=g a^2, 2N=g(a^2+b^2)
bad=0; cnt=0
for _ in range(40):
    N=math.prod(random.sample(split,random.randint(2,5)))*random.choice([1,5,13])
    P=gpoints(N)
    if not P: continue
    z1=random.choice(P)
    for z in P:
        if z==z1: continue
        d2=(z[0]-z1[0])**2+(z[1]-z1[1])**2
        if d2%2: bad+=1; continue
        m=d2//2; w=mul(z,conj(z1)); v=w[1]
        ok = (m*(2*N-m)==v*v) and (N-m==w[0])
        g=gcd(m,2*N); a2=m//g; b2=(2*N-m)//g
        ok = ok and m%g==0 and isqrt(a2)**2==a2 and isqrt(b2)**2==b2 and g*(a2+b2)==2*N
        cnt+=1; bad+= (not ok)
print("(a) checked",cnt,"pairs, failures",bad)
# ---- (d): z == z' mod q  <=>  (z/z') == 1 in (Z[i]/q)^*, and q | z-z' => q | a where m=g a^2 (q odd, coprime to N)
bad=0;cnt=0;hits=0
for _ in range(30):
    N=math.prod(random.sample(split[:12],random.randint(3,6)))
    P=gpoints(N)
    for q in random.sample([l for l in sympy.primerange(3,200) if N%l],5):
        for z1 in random.sample(P,min(6,len(P))):
            for z in P:
                if z==z1: continue
                cong = (z[0]-z1[0])%q==0 and (z[1]-z1[1])%q==0
                # z conj z1 == N mod q  <=> z/z1 == 1 mod q
                w=mul(z,conj(z1)); cong2 = (w[0]-N)%q==0 and w[1]%q==0
                d2=(z[0]-z1[0])**2+(z[1]-z1[1])**2; m=d2//2; g=gcd(m,2*N); a=isqrt(m//g)
                cong3 = (a%q==0)
                cnt+=1; hits+=cong; bad+= not (cong==cong2==cong3)
print("(d) checked",cnt,"pairs (",hits,"congruent ), failures",bad)
