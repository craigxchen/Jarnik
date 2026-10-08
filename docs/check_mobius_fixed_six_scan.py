from fractions import Fraction
from itertools import product
from math import gcd

def central(ts):
 ws=[]
 for i,t in enumerate(ts):
  D=Fraction(1)
  for j,s in enumerate(ts):
   if i!=j:D*=t-s
  ws.append(Fraction((t*t+1)**2,D))
 den=1
 for w in ws: den=den*w.denominator//gcd(den,w.denominator)
 ints=[w.numerator*(den//w.denominator) for w in ws]
 g=0
 for x in ints:g=gcd(g,abs(x))
 return [x//g for x in ints]
def evalmap(ts,a,b,c,d):
 out=[]
 for t in ts:
  den=c*t+d
  if den==0:return None
  out.append((a*t+b)/den)
 return out
for base0 in [[2,3,5,12,17,23],[1,2,3,4,7,11],[2,4,7,12,19,31]]:
 ts0=[Fraction(x) for x in base0]; best=None
 for a,b,c,d in product(range(-7,8),repeat=4):
  if a*d-b*c==0:continue
  ts=evalmap(ts0,Fraction(a),Fraction(b),Fraction(c),Fraction(d))
  if ts is None or len(set(ts))<6 or min(ts)<=0:continue
  u=central(ts); ratio=Fraction(max(map(abs,u)),1)/min(ts)
  if best is None or ratio<best[0]:best=(ratio,(a,b,c,d),ts,u)
 print(base0, best)
