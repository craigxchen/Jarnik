import math
import numpy as np
from math import gcd
th=math.atan2(3,2)
def ideals(X):
    out=[]
    r=int(math.isqrt(X))+1
    for a in range(1,r+1):
        for b in range(0,r+1):
            n=a*a+b*b
            if n>X: break
            if gcd(a,b)!=1: continue
            # signed 13-exponent wrt pi=2+3i
            e=0; x,y=a,b
            # divide by 2+3i : (x+iy)/(2+3i) = (x+iy)(2-3i)/13
            while True:
                u,v=2*x+3*y, 2*y-3*x
                if u%13==0 and v%13==0: x,y=u//13,v//13; e+=1
                else: break
            while True:
                u,v=2*x-3*y, 2*y+3*x
                if u%13==0 and v%13==0: x,y=u//13,v//13; e-=1
                else: break
            out.append((math.atan2(b,a),e,n))
    return out
U=math.pi/2
for X in [10**3,10**4,10**5,4*10**5]:
    L=ideals(X)
    for name,f in [('genuine',lambda t,e: t),('mult graft',lambda t,e: t-e*th/2)]:
        A=np.array([f(t,e)%U for t,e,n in L]); idx=np.argsort(A); A=A[idx]
        d=np.diff(np.concatenate([A,[A[0]+U]]))
        print(X,name,'#',len(L),'X*minspacing=',X*d.min())
# Gram (dual large sieve) norm vs 2*(2K+pi X/2)
for X,K in [(300,10),(300,100),(1000,30),(1000,300)]:
    L=ideals(X)
    for name,f in [('genuine',lambda t,e: t),('mult graft',lambda t,e: t-e*th/2)]:
        A=np.array([4*f(t,e) for t,e,n in L])
        D=A[:,None]-A[None,:]
        G=np.ones_like(D)
        for k in range(1,K+1): G+=2*np.cos(k*D)
        lam=np.linalg.eigvalsh(G).max()
        print(X,K,name,'norm=',round(lam,1),' 2K+piX/2=',round(2*K+math.pi*X/2,1),' 4K+piX=',round(4*K+math.pi*X,1))
