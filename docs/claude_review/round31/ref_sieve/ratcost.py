# Independent search: rational graft direction psi = arg b, height h; min over primitive alpha,beta (norms coprime to p)
# of sin||psi+arg beta-arg alpha||_u * max(h N beta, N alpha); compare with sqrt(h/Q).
import math
from math import gcd
from fractions import Fraction
def mul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def egcd(a,b):
    if b==0: return (a,1,0)
    g,x,y=egcd(b,a%b); return (g,y,x-(a//b)*y)
def run(base,k,p,B=40):
    b=(1,0)
    for _ in range(k): b=mul(b,base)
    Q=b[0]**2+b[1]**2; h=p**k
    best=None
    for x in range(-B,B+1):
        for y in range(0,B+1):
            if (x,y)==(0,0) or gcd(x,y)!=1 or (x*x+y*y)%p==0: continue
            v=mul(b,(x,y)); g=gcd(abs(v[0]),abs(v[1])); vp=(v[0]//g,v[1]//g)
            Nb=x*x+y*y; nv2=v[0]**2+v[1]**2
            # alpha with det(vp,alpha)=1: vp0*a1 - vp1*a0 = 1
            gg,s,t=egcd(vp[0],-vp[1])  # vp0*s + (-vp1)*t = gg = +-1
            if gg<0: s,t=-s,-t
            a0=(t,s)   # det(vp,(t,s)) = vp0*s - vp1*t = 1
            # alpha = a0 + j vp ; choose j so that |alpha|^2 ~ h*Nb*... target |alpha|^2 = sqrt(h/Q)*|v|*|beta|... use r
            target=math.sqrt(h/Q)*math.sqrt(nv2)   # |alpha| target = r*|v|
            # alpha component along vp: (a0.vp)/|vp| + j|vp| ; perpendicular comp = 1/|vp|
            nvp=math.sqrt(vp[0]**2+vp[1]**2)
            c0=(a0[0]*vp[0]+a0[1]*vp[1])/nvp
            for sgn in (1,-1):
                j0=round((sgn*target-c0)/nvp)
                for j in range(j0-2,j0+3):
                    al=(a0[0]+j*vp[0],a0[1]+j*vp[1]); Na=al[0]**2+al[1]**2
                    if Na%p==0 or gcd(*al)!=1: continue
                    det=abs(v[0]*al[1]-v[1]*al[0]); dot=abs(v[0]*al[0]+v[1]*al[1])
                    num=min(det,dot)
                    if num==0: continue
                    # normalized^2 = num^2 * max(h Nb, Na)^2 / (nv2*Na)
                    val=Fraction(num*num*max(h*Nb,Na)**2, nv2*Na)
                    if best is None or val<best[0]: best=(val,(x,y),al)
    ratio=best[0]/Fraction(h,Q)
    print(f"base={base} k={k}: Q={Q} h={h}  min normalized^2/(h/Q) = {float(ratio):.5f}  (beta,alpha)={best[1]},{best[2]}")
for k in (2,4,8):
    run((3,4),k,5)
for k in (2,4):
    run((3,2),k,5)
