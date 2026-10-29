# Exact check of the mnn31 "C = 7/2" five-divisor Pell family (Lean file ErdosLean/Erdos887.lean, FiveDivisors)
# n = u(u+2)(u+3)(u+5), (u+2)(u+3) = 3(w+1)(w+2); divisors d1..d5 as in the Lean 'key' lemma.
# Computes: torus height H = (max/min) lcm(d/g), dist = (max-min)/min, gamma = log H / (-log dist), C_actual.
from math import gcd, log, isqrt
from fractions import Fraction
def lcm(a,b): return a//gcd(a,b)*b
u,w = 1843,1064
print("k  log10(n)  C_actual  gcd  log(lcm)/log(n)  gamma  log(dist^{7/2}H)")
for k in range(40):
    assert (u+2)*(u+3)==3*(w+1)*(w+2)
    n=u*(u+2)*(u+3)*(u+5)
    d=[(u+2)*(u+3),3*(w+1)*(w+3),(u+2)*(u+5),(u+3)*(u+5),3*(w+2)*(w+3)]
    assert all(n%x==0 for x in d) and len(set(d))==5
    r=isqrt(n); assert all(x*x>n for x in d)
    g=0
    for x in d: g=gcd(g,x)
    v=[x//g for x in d]
    L=1
    for x in v: L=lcm(L,x)
    lH=log(L)+log(max(v)/min(v))
    D=max(v)-min(v); ldist=log(D)-log(min(v))
    gam=lH/(-ldist)
    Cact=(max(d)-r)/n**0.25 if k<25 else float('nan')
    if k%3==0 or k<4:
        print(k, round(log(n,10),1), round(Cact,5) if k<25 else '-', g, round(log(L)/log(n),4), round(gam,5), round(3.5*ldist+lH,3))
    u,w = 2*u+3*w+7, u+2*w+4
