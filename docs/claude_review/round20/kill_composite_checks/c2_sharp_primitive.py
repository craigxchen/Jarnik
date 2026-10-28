# F2 check: sharp large sieve on PRIMITIVE support <=> spacing c/X of all angles (two-point argument / Montgomery-Vaughan).
# Compare normalised min spacing  X * min |sin(angle_i - angle_j)|  over primitive ideals of norm <= X under
#   genuine angles, the MINIMAL graft (only S-smooth ideals moved), and the MULTIPLICATIVE graft (every pi^u*beta rotated).
# One fake prime p=101, pi=10+i; fake angle psi with tan psi = sqrt(7)-2 (quadratic irrational, so every j*psi is BA).
import numpy as np, math, sys
p=101; a,b=10,1; theta=math.atan2(b,a)
psi=math.atan(math.sqrt(7)-2)
def prim_ideals(X):
    r=math.isqrt(X); xs=[];ys=[]
    for x in range(1,r+1):
        ymax=math.isqrt(X-x*x); y=np.arange(0,ymax+1,dtype=np.int64)
        xs.append(np.full(len(y),x,dtype=np.int64)); ys.append(y)
    x=np.concatenate(xs); y=np.concatenate(ys)
    keep=np.gcd(x,y)==1
    return x[keep],y[keep]
def val(x,y,a,b):
    v=np.zeros(len(x),dtype=np.int64); act=np.ones(len(x),bool); x=x.copy(); y=y.copy()
    while act.any():
        re=x*a+y*b; im=y*a-x*b
        ok=act&(re%p==0)&(im%p==0)
        if not ok.any(): break
        x=np.where(ok,re//p,x); y=np.where(ok,im//p,y); v+=ok; act=ok
    return v,x,y
def minsp(ang):
    s=np.sort(np.mod(ang,math.pi/2)); d=np.diff(np.concatenate([s,[s[0]+math.pi/2]]))
    i=np.argmin(d); return d[i]
for X in [2000,10**4,5*10**4,2*10**5,10**6,4*10**6]:
    x,y=prim_ideals(X); th=np.arctan2(y,x)
    u,_,_=val(x,y,a,b); w,_,_=val(x,y,a,-b)
    D=u-w                     # primitive: at most one of u,w nonzero
    n=x*x+y*y; smooth=(n==p**(u+w))
    mult=th+D*(psi-theta)
    minimal=np.where(smooth,mult,th)
    g=X*np.sin(minsp(th)); m1=X*np.sin(minsp(minimal)); m2=X*np.sin(minsp(mult))
    print(f"X={X:8d} #prim={len(x):8d} #S-div={int((D!=0).sum()):6d}  X*minspacing: genuine={g:.4f}  minimal={m1:.4f}  multiplicative={m2:.3e}",flush=True)
