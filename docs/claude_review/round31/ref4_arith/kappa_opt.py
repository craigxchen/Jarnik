# For kappa e_orb >= (K+B)^2, search section weight a and chain weight b so that tangency chains
# (first exceptional weight b, later ones weight b) cost <= unit per intersection; report best m bound.
from fractions import Fraction as Fr
def F(K,a,k,b0,b): return K*(1-2*b+b*b0)+K*a*k*(b-b0)+(1+b-b0-a*k)**2-2*a*a*k*(k-1)/2
def best(K,steps=200):
    res=None
    for i in range(1,steps+1):
        a=Fr(i,steps)
        u=(K-2)*a*a; coef=a*(K*a-4)/8
        if coef>=0: continue
        for j in range(0,steps+1):
            b=Fr(j,steps)
            f1=F(K,a,2,0,b); f2=F(K,a,2,b,b)
            if f2<=u and f1<=u+ (u-f2) and f1<=u:   # e=1 not a chain; e>=2: f1+(e-1)f2<=e u
                pass
            ok = all(f1+(e-1)*f2<=e*u for e in range(2,60))
            if ok:
                mb=(K+1)/(-coef)
                if res is None or mb<res[0]: res=(mb,a,b)
    return res
for K in [3,Fr(7,2),4,5,8]:
    print(K,best(K))
