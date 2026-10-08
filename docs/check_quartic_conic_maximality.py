#!/usr/bin/env python3
"""Exact certificate for quartic_conic_family_maximality.md.

Uses only the Python standard library. Polynomial lists are in increasing
coefficient order. R is Q(v); C is Q(v)(i). No floating point arithmetic,
parameter sampling, or external symbolic algebra system is used.

The rank test checks gcds of coefficient minors, after removing explicitly
listed factors whose real zeros are already excluded degenerations.
"""
from fractions import Fraction as Q

def tr(p):
 p=list(p)
 while len(p)>1 and p[-1]==0:p.pop()
 return p
def pa(a,b):return tr([(a[j] if j<len(a) else 0)+(b[j] if j<len(b) else 0) for j in range(max(len(a),len(b)))])
def pm(a,b):
 o=[Q(0)]*(len(a)+len(b)-1)
 for j,x in enumerate(a):
  for k,y in enumerate(b):o[j+k]+=x*y
 return tr(o)
def pd(a,b):
 a=tr(a);o=[Q(0)]*max(1,len(a)-len(b)+1)
 while a!=[0] and len(a)>=len(b):
  k=len(a)-len(b);c=a[-1]/b[-1];o[k]=c
  for j in range(len(b)):a[j+k]-=c*b[j]
  a=tr(a)
 return tr(o),a
def pg(a,b):
 while b!=[0]:a,b=b,pd(a,b)[1]
 return [x/a[-1] for x in a]
class R:
 def __init__(self,a=0,b=None):
  if isinstance(a,R):self.a,self.b=a.a,a.b;return
  if not isinstance(a,list):a=[Q(a)]
  if b is None:b=[Q(1)]
  a=tr(a);b=tr(b);assert b!=[0]
  if a==[0]:self.a,self.b=[Q(0)],[Q(1)];return
  g=pg(a,b);a=pd(a,g)[0];b=pd(b,g)[0];lead=b[-1];self.a=[x/lead for x in a];self.b=[x/lead for x in b]
 def __add__(self,o):o=R(o);return R(pa(pm(self.a,o.b),pm(o.a,self.b)),pm(self.b,o.b))
 __radd__=__add__
 def __neg__(self):return R([-x for x in self.a],self.b)
 def __sub__(self,o):return self+-R(o)
 def __rsub__(self,o):return R(o)+-self
 def __mul__(self,o):o=R(o);return R(pm(self.a,o.a),pm(self.b,o.b))
 __rmul__=__mul__
 def __truediv__(self,o):o=R(o);return R(pm(self.a,o.b),pm(self.b,o.a))
 def __rtruediv__(self,o):return R(o)/self
 def __pow__(self,n):
  o=R(1)
  for j in range(n):o=o*self
  return o
 def __repr__(self):return str(self.a)+' / '+str(self.b)
 def __eq__(self,o):o=R(o);return self.a==o.a and self.b==o.b
class C:
 def __init__(self,a=0,b=0):
  if isinstance(a,C):self.x,self.y=a.x,a.y
  else:self.x,self.y=R(a),R(b)
 def __add__(self,o):o=C(o);return C(self.x+o.x,self.y+o.y)
 __radd__=__add__
 def __neg__(self):return C(-self.x,-self.y)
 def __sub__(self,o):return self+-C(o)
 def __mul__(self,o):o=C(o);return C(self.x*o.x-self.y*o.y,self.x*o.y+self.y*o.x)
 __rmul__=__mul__
 def conj(self):return C(self.x,-self.y)
 def __truediv__(self,o):o=C(o);n=o.x**2+o.y**2;r=self*o.conj();return C(r.x/n,r.y/n)
 def __repr__(self):return '('+str(self.x)+') + i('+str(self.y)+')'
def norm(r):return [r.x**2+r.y**2,-2*r.x,R(1)]
def mulpoly(a,b):
 o=[R(0)]*(len(a)+len(b)-1)
 for j,x in enumerate(a):
  for k,y in enumerate(b):o[j+k]=o[j+k]+x*y
 return o

from itertools import combinations
import math,json,time
# primitive rational polynomial helper
def primitive(p):
 p=tr(p)
 if p==[0]:return p
 dl=1
 for x in p:dl=math.lcm(dl,x.denominator)
 ints=[int(x*dl) for x in p];g=0
 for x in ints:g=math.gcd(g,abs(x))
 if ints[-1]<0:g=-g
 return [Q(x//g) for x in ints]
def scalar(a,c):return tr([x*c for x in a])
def ps(a,b):return pa(a,scalar(b,-1))
def normvec(r):
 D=pd(pm(r.x.b,r.y.b),pg(r.x.b,r.y.b))[0]
 X=pm(r.x.a,pd(D,r.x.b)[0]);Y=pm(r.y.a,pd(D,r.y.b)[0])
 out=[pa(pm(X,X),pm(Y,Y)),scalar(pm(D,X),-2),pm(D,D)]
 g=pg(pg(out[0],out[1]),out[2])
 out=[pd(p,g)[0] for p in out]
 # one common numerical denominator; scalar normalization irrelevant
 dl=1
 for p in out:
  for x in p:dl=math.lcm(dl,x.denominator)
 return [scalar(p,dl) for p in out]
def vvprod(a,b):
 o=[[Q(0)] for j in range(len(a)+len(b)-1)]
 for j,p in enumerate(a):
  for k,q in enumerate(b):o[j+k]=pa(o[j+k],pm(p,q))
 return o
def pad(v):return v+[[Q(0)]]*(9-len(v))
allowed=[[Q(0),Q(1)],[Q(-4),0,Q(1)],[Q(-16),0,Q(1)]]+[[Q(c),0,Q(1)] for c in [2,4,8,16,32]]
def strip(p):
 p=primitive(p)
 for f in allowed:
  while len(p)>1:
   quo,rem=pd(p,f)
   if rem!=[0]:break
   p=primitive(quo)
 return p

def rankcondition(p,q,d):
 # all 3x3 minors vanish iff the three coefficient vectors have rank <=2
 g=None;nonzero=0
 for i,j,k in combinations(range(9),3):
  if all(p[t]==q[t]==d[t]==[0] for t in [i,j,k]):continue
  det=pa(ps(pm(p[i],ps(pm(q[j],d[k]),pm(q[k],d[j]))),pm(q[i],ps(pm(p[j],d[k]),pm(p[k],d[j])))),pm(d[i],ps(pm(p[j],q[k]),pm(p[k],q[j]))))
  if det==[0]:continue
  nonzero+=1
  g=strip(det) if g is None else primitive(pg(g,det))
  if len(g)==1:return None
 return [] if g is None else g


v=R([Q(0),Q(1)]);s=6*v/(v*v+8)
for branch in [-1,1]:
 a,b,c,d=C(0,1),C(3*s,2),C(4*s,-1),C(s,-2)
 if branch==-1:
  e=C(8*(v*v-4)/(v*(v*v+8)),-3*v*v/(v*v+8))
  f=C(v*(v*v+32)/(2*(v*v+8)),-24/(v*v+8))
  g=C(16*(v*v+2)/(v*(v*v+8)),3*v*v/(v*v+8))
 else:
  e=C(v*(16-v*v)/(2*(v*v+8)),-24/(v*v+8))
  f=C(16*(v*v+2)/(v*(v*v+8)),-3*v*v/(v*v+8))
  g=C(v*(v*v+32)/(2*(v*v+8)),24/(v*v+8))
 names=['ABCD','ABEF','ACEG'];blocks=dict(zip('ABCDEFG',[a,b,c,d,e,f,g]));nv={key:normvec(z) for key,z in blocks.items()}
 dirs=[]
 for name in names:
  dd=[]
  for size in range(5):
   for sub in combinations(name,size):
    p=[[Q(1)]]
    for key in sub:p=vvprod(p,nv[key])
    dd.append((''.join(sub),pad(p)))
  dirs.append(dd)
 print('Checking branch',branch,flush=True)
 generic=[];exceptions=[]
 for aa,bb in combinations(range(3),2):
  sh=set(names[aa])&set(names[bb]);delta=[[Q(1)]]
  for key in sorted(sh):delta=vvprod(delta,nv[key])
  delta=pad(delta)
  for namep,p in dirs[aa]:
   for nameq,q in dirs[bb]:
    result=rankcondition(p,q,delta)
    if result==[]:generic.append([aa,bb,namep,nameq])
    elif result is not None:exceptions.append([aa,bb,namep,nameq,result])

 # Classify all identically rank-two direction triples.
 expected=set()
 nontrivial={(0,1):[('AC','AE'),('BD','BF')],(0,2):[('AB','AE'),('CD','CG')],(1,2):[('AB','AC'),('EF','EG')]}
 for aa,bb in combinations(range(3),2):
  shared=''.join(sorted(set(names[aa])&set(names[bb])))
  for np,_ in dirs[aa]:
   for nq,_ in dirs[bb]:
    if np==shared or nq==shared or (np==nq and set(np).issubset(set(shared))):expected.add((aa,bb,np,nq))
  for np,nq in nontrivial[(aa,bb)]:expected.add((aa,bb,np,nq))
 assert set(map(tuple,generic))==expected
 assert len(generic)==108 and len(exceptions)==77
 # Every exceptional residual gcd either has positive even coefficients,
 # or is v^2-8, where roots F and G are conjugates.
 exceptional_counts={}
 for aa,bb,np,nq,p in exceptions:
  pp=tuple(p);exceptional_counts[pp]=exceptional_counts.get(pp,0)+1
  if p==[Q(-8),Q(0),Q(1)]:continue
  assert p[0]>0 and all(x>=0 for x in p) and all(p[j]==0 for j in range(1,len(p),2))
 assert exceptional_counts[(Q(-8),Q(0),Q(1))]==3
 expected_exception_polynomials={
  (Q(-8),Q(0),Q(1)),
  (Q(16),Q(0),Q(5)),
  (Q(20),Q(0),Q(1)),
  (Q(64),Q(0),Q(52),Q(0),Q(1)),
  tuple(map(Q,[1024,0,1152,0,276,0,5] if branch==-1 else [1280,0,1104,0,72,0,1])),
 }
 assert set(exceptional_counts)==expected_exception_polynomials
 print('  768 direction pairs: 108 generic, 3 degenerate real exceptions; all other exceptional roots nonreal.',flush=True)
 # Check the only nontrivial generic candidates against the remaining row.
 tau=12*s*(1+s*s)
 la=(e-c)*(e-d)/((e-a.conj())*(e-b.conj()))
 mu=(e-b)*(e-d)/((e-a.conj())*(e-c.conj()))
 assert la.y==mu.y==0
 la,mu=la.x,mu.x
 Fs=[[(19*s*s+4)/tau,-(12*s**3+20*s)/tau,(19*s*s+5)/tau,-8*s/tau,1/tau]]
 for lam,qq in [(la,mulpoly(norm(a),norm(b))),(mu,mulpoly(norm(a),norm(c)))]:Fs.append([Fs[0][j]-lam*qq[j]/tau for j in range(5)])
 fv=C(-Fs[1][3]/Fs[1][4])-a-b-e
 gv=C(-Fs[2][3]/Fs[2][4])-a-c-e
 assert (fv.x,fv.y)==(f.x,f.y) and (gv.x,gv.y)==(g.x,g.y)
 triples=[(0,1,2,b,d,f),(0,2,1,c,d,g),(1,2,0,e,f,g)]
 expected_common=pm([Q(-16 if branch==-1 else -4),0,Q(1)],[Q(4 if branch==-1 else 16),0,Q(1)])
 for aa,bb,cc,common,r1,r2 in triples:
  p=mulpoly(norm(common),norm(r1));q=mulpoly(norm(common),norm(r2));delta=[Fs[bb][j]-Fs[aa][j] for j in range(5)]
  lam=None
  for j in range(4):
   det=-p[j]*q[4]+p[4]*q[j]
   if det!=0:
    lam=(-delta[j]*q[4]+delta[4]*q[j])/det;mm=(p[j]*delta[4]-p[4]*delta[j])/det;break
  assert lam is not None and all(lam*p[j]-mm*q[j]==delta[j] for j in range(5))
  lead=Fs[aa][4]+lam*p[4]-Fs[cc][4]
  cub=Fs[aa][3]+lam*p[3]-Fs[cc][3]
  joint=primitive(pg(lead.a,cub.a))
  expected_joint=pm(expected_common,[Q(8),0,Q(1)])
  if (aa,bb)!=(1,2):expected_joint=pm(expected_joint,[Q(8),0,Q(1)])
  assert joint==primitive(expected_joint),(branch,aa,bb,joint)
 print('  All three companion cancellation gcds have only degenerate real roots.',flush=True)
print('PASS: both real conic branches are maximal whenever the seven factors remain nondegenerate.',flush=True)
