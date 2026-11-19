from fractions import Fraction as F
from itertools import combinations
from math import isqrt
pr=(5,13,17); base={5:(2,1),13:(3,2),17:(4,1)}
rows=((0,0,0),(0,1,0),(1,0,0),(1,1,1)); vals=((0,0,0),(0,1,0),(1,0,0),(-1,1,1))
def cm(a,b): return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def unitw(v):  # prod (pi/conj pi)^v as complex Fraction
    num=(F(1),F(0))
    for p,e in zip(pr,v):
        a,b=base[p]; q=(F(a*a-b*b,p),F(2*a*b,p))  # pi/conj(pi)
        if e<0: q=(q[0],-q[1])
        for _ in range(abs(e)): num=cm(num,q)
    return num
w=[unitw(v) for v in vals]
for i,j in combinations(range(4),2):
    d=1
    for p,a,b in zip(pr,rows[i],rows[j]): d*=p**abs(a-b)
    cosang=w[i][0]*w[j][0]+w[i][1]*w[j][1]
    t2=d*(1-cosang)/2
    print((i,j),d,t2, "t =",F(isqrt(t2.numerator),isqrt(t2.denominator)) if isqrt(t2.numerator)**2==t2.numerator and isqrt(t2.denominator)**2==t2.denominator else "irrational")
