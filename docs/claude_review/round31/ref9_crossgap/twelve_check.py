# measure of good set G_p = {y in R/Z : |sin((pi/2)y - arg w)|>=1/p for all 0!=w, Norm w<=p}
import math, sys
from math import gcd
def good_measure(p):
    R=math.isqrt(p); cs=[]
    for a in range(1,R+1):          # primitive directions mod units: a>0,b>=0
        for b in range(0,R+1):
            if a*a+b*b<=p and gcd(a,b)==1:
                cs.append((2/math.pi)*math.atan2(b,a))  # in [0,1)
    half=(2/math.pi)*math.asin(1/p)
    iv=sorted(((c-half)%1.0, ((c-half)%1.0)+2*half) for c in cs)
    # union length on circle
    tot=0; cur_s,cur_e=None,None
    for s,e in iv:
        if cur_e is None or s>cur_e: 
            if cur_e is not None: tot+=cur_e-cur_s
            cur_s,cur_e=s,e
        else: cur_e=max(cur_e,e)
    tot+=cur_e-cur_s
    # wrap overlap
    if cur_e>1: 
        first=iv[0][0]; tot-=max(0,min(cur_e-1,iv[0][1])-first) if cur_e-1>first else 0
    return len(cs), 1-tot, len(cs)*2*half
for p in [101,1009,10007,100003,1000003]:
    n,g,ub=good_measure(p)
    print(p, "dirs",n, "pred",round(6/math.pi**2*math.pi*p/4,1), "good measure",round(g,5),"1-unionbound",round(1-ub,5))
r=math.sin(.38*math.pi)/(.38*math.pi); print("r",r,"r^10/.38",r**10/.38)
# good log fraction in finite test
print("b",(1-math.log(64)/math.log(100))/2)
