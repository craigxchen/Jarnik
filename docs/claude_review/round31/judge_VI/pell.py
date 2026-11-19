# Integer analogue of VI.1 fails: z1 = gamma^2 (2+i), z0 = conj(z1), |z1-z0| = 2.
from math import log
def gm(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def vp(n,p):
    k=0
    while n%p==0: n//=p; k+=1
    return k
# solutions of u^2-5y^2=+-1, x=u-2y
sols=[]; u,y=1,0
# fundamental unit 2+sqrt5 (norm -1): (u,y)->(2u+5y, u+2y)
for _ in range(12):
    u,y=2*u+5*y,u+2*y; sols.append((u-2*y,y))
for x,y in sols:
    g=(x,y); z1=gm(gm(g,g),(2,1)); z0=(z1[0],-z1[1])
    N=z1[0]**2+z1[1]**2
    assert N==(x*x+y*y)**2*5
    # parity word of z1/z0 = (gamma/gammabar)^2 * (2+i)/(2-i): odd valuations only at 5
    print("x,y=",x,y," |z1-z0| =",abs(2*z1[1]),"  W=log N =%.2f"%log(N), " v5(N)=",vp(N,5))
