# [J] check: truncated square root Q of F(Y)=prod(Y+k_nu), deg(F-Q^2)<=2, 16Q integral
import sympy as sp, itertools, random
Y=sp.symbols('Y'); k=sp.symbols('k1:7')
F=sp.prod([Y+ki for ki in k])
e=[sp.Poly(F,Y).coeff_monomial(Y**(6-j)) for j in range(7)]
e1,e2,e3=e[1],e[2],e[3]
Q=Y**3+e1/2*Y**2+(e2/2-e1**2/8)*Y+e3/2-e1*e2/4+e1**3/16
D=sp.Poly(sp.expand(F-Q**2),Y)
print("deg(F-Q^2) =",D.degree())
# brute force small search: X up to 4000, window K=floor(X^(1/6))*3, min |t|/sqrt(X) for products = s^2+t^2, t!=0
import math
def sumsq_small_t(n,tmax):
    out=[]
    for t in range(1,tmax+1):
        r=n-t*t
        if r<0: break
        s=math.isqrt(r)
        if s*s==r: out.append(t)
    return out
best=1e9
random.seed(1)
for X in range(64,6000,7):
    K=max(1,int(round(2*X**(1/6))))
    for trial in range(30):
        ks=[0]+[random.randint(0,K) for _ in range(5)]
        n=1
        for kk in ks: n*=X+kk
        ts=sumsq_small_t(n,int(2*math.sqrt(X)))
        for t in ts: best=min(best,t/math.sqrt(X))
print("min |t|/sqrt X over sample (window 2X^(1/6)) =",best)
