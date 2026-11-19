from fractions import Fraction as Fr
from math import isqrt
def issq(n):
    if n<0: return False
    r=isqrt(n); return r*r==n
def sqfr(q):
    q=Fr(q)
    if q<0: return None
    n,d=q.numerator,q.denominator
    if issq(n) and issq(d): return Fr(isqrt(n),isqrt(d))
    return None
res=[]
for r in range(1,6):
  for s in range(r+1,8):
    pts=[]
    for v in range(1,13):
      for u in range(-60*v*v,60*v*v+1):
        x=Fr(u,v*v)
        if x.denominator!=v*v: continue
        y=sqfr(x*(x+r*r)*(x+s*s))
        if y is not None and y!=0: pts.append((x,y))
    tors={Fr(r*s),Fr(-r*s)}
    nt=[p for p in pts if p[0] not in tors]
    if nt: res.append((r,s,nt[:4]))
for z in res[:6]: print(z)
