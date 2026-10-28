"""Independent exact check of the m-czech member (a,b,c)=(1,2,4): along integer solutions of z^2=56x^2+65,
compute the 5 divisors, primitive tuple u, H=(max/min)*lcm(u) (exact toric height), dist=(max-min)/min,
gamma=log H/(-log dist), log(dist^{15/4} H), log lcm/log n, and the actual window constant C."""
from fractions import Fraction
from math import gcd, isqrt, log
def lcm(p,q): return p//gcd(p,q)*q
# all solutions with 0<=x<1e6 by brute force, then multiply by the unit 15+2sqrt56
base=[(X,isqrt(56*X*X+65)) for X in range(0,10**5) if isqrt(56*X*X+65)**2==56*X*X+65]
print("base solutions:",base[:8])
sols=set()
for X,Z in base[:6]:
    for sgn in (1,-1):
        x_,z_=X,sgn*Z
        for _ in range(30):
            if x_>0 and z_>0: sols.add((x_,z_))
            x_,z_=15*x_+2*z_,15*z_+112*x_
sols=sorted(sols)
rows=[]
for X,Z in sols:
    if X<20: continue
    n=(X*X-1)*(X*X-4)*(X*X-16)*(X*X-49)
    F=[(X-1)*(X-2)*(X-4)*(X+7), X**4-35*X*X+56+2*X*Z, (X*X-4)*(X*X-16), (X*X-1)*(X*X-16), (X*X-1)*(X*X-4)]
    assert all(n%f==0 for f in F)
    s=isqrt(n); q4=isqrt(s)
    assert all(f> s for f in F)
    C=(max(F)-s)*10**6//q4/10**6
    g=0
    for f in F: g=gcd(g,f)
    u=sorted(f//g for f in F)
    L=1
    for v in u: L=lcm(L,v)
    lH=log(L)+log(u[-1])-log(u[0])
    ld=log(u[-1]-u[0])-log(u[0])
    rows.append((X,C,g,log(L)/log(n),lH/(-ld),15/4*ld+lH,lH))
for r in rows[::max(1,len(rows)//25)]+[rows[-1]]:
    print("x=%.3e C=%.4f gcd=%d logL/logn=%.5f gamma=%.5f log(dist^3.75 H)=%.4f logH=%.1f"%r)
print("range of log(dist^{15/4}H) over",len(rows),"points:",min(r[5] for r in rows),max(r[5] for r in rows))
first=[r for r in rows if r[4]<4]
print("first integral Pell point with gamma<4: x=%.3e logH=%.1f gamma=%.4f"%(first[0][0],first[0][6],first[0][4]))
print("max log H reached:",max(r[6] for r in rows))
