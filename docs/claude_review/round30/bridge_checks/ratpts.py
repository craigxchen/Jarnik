# rational points on the cube surface a b c d = 1, a+b+c+d = 0, height <= H, nondegenerate for some eps-assignment
import sys, math, itertools
from fractions import Fraction as Fr
H=int(sys.argv[1])
sols=set()
for q in range(1,H+1):
  for p in range(-H,H+1):
    if p==0 or math.gcd(p,q)!=1: continue
    for s in range(1,H+1):
      for r in range(-H,H+1):
        if r==0 or math.gcd(r,s)!=1: continue
        if (p*s,q*r)>(r*q,s*p) : pass
        X=p*r*(p*r*(p*s+r*q)**2-4*q**3*s**3)
        if X<0: continue
        t=math.isqrt(X)
        if t*t!=X: continue
        a=Fr(p,q); b=Fr(r,s); A=a*b; B=a*b*(a+b)
        sq=Fr(t, (q*s)**2)   # sqrt(Delta)
        for c in ((-B+sq)/(2*A),(-B-sq)/(2*A)):
            y=tuple(sorted((a,b,c,-(a+b+c))))
            if 0 in y: continue
            assert y[0]*y[1]*y[2]*y[3]==1
            sols.add(y)
def nondeg(y):
    ok=[]
    for plus in itertools.combinations(range(4),2):
        minus=[k for k in range(4) if k not in plus]
        yy=[y[plus[0]],y[minus[0]],y[minus[1]],y[plus[1]]]
        R=[]
        for e,v in zip([1,-1,-1,1],yy):
            x=(Fr(e)/v, v); R+= [x,(-x[0],-x[1]),(x[0],-x[1]),(-x[0],x[1])]
        if len(set(R))==16: ok.append(yy)
    return ok
out=[]
for y in sols:
    nd=nondeg(y)
    if nd: out.append((max(max(abs(v.numerator),v.denominator) for v in y),y))
out.sort()
for h,y in out: print(h,' '.join(str(v) for v in y))
