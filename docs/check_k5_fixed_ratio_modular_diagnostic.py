"""Finite-field K5 fixed-ratio diagnostic.

This is an exact computation over F_1000003 for x01=1, x02=2.  It is not a
characteristic-zero Groebner certificate.  The reduced basis makes the partial
product P9=prod(j=1..9)(x0^2-xj^2) reduce to zero.
"""

from itertools import combinations
from fractions import Fraction
from collections import deque
P=1000003
n=4
# polys mod p dict mon tuple -> int
def norm(f):
 f={m:c%P for m,c in f.items() if c%P};return f
def add(*ps):
 d={}
 for f in ps:
  for m,c in f.items():d[m]=(d.get(m,0)+c)%P
 return norm(d)
def scale(f,a):return norm({m:c*a for m,c in f.items()})
def mul(f,g):
 d={}
 for a,x in f.items():
  for b,y in g.items():
   m=tuple(i+j for i,j in zip(a,b));d[m]=(d.get(m,0)+x*y)%P
 return norm(d)
def powp(f,k):
 r={(0,)*n:1}
 for _ in range(k):r=mul(r,f)
 return r
def lin(c,co):
 d={(0,)*n:c}
 for j,v in co.items():
  m=tuple(1 if i==j else 0 for i in range(n));d[m]=v
 return norm(d)
def lm(f):return max(f)
def divisible(a,b):return all(x>=y for x,y in zip(a,b))
def submono(a,b):return tuple(x-y for x,y in zip(a,b))
def mulmono(f,m,a=1):return {tuple(x+y for x,y in zip(k,m)):c*a%P for k,c in f.items()}
def reduce_poly(f,G):
 f=dict(f);r={}
 while f:
  m=max(f); c=f[m]; found=False
  for g in G:
   mg=lm(g)
   if divisible(m,mg):
    q=submono(m,mg); inv=pow(g[mg],P-2,P)
    for k,v in mulmono(g,q,c*inv).items():
     f[k]=(f.get(k,0)-v)%P
     if not f[k]:del f[k]
    found=True;break
  if not found:r[m]=c;del f[m]
 return norm(r)
def monic(f):
 z=f[lm(f)];return scale(f,pow(z,P-2,P))
def spol(f,g):
 a,b=lm(f),lm(g); l=tuple(max(x,y) for x,y in zip(a,b));
 return add(mulmono(f,submono(l,a),1),scale(mulmono(g,submono(l,b),1),-1))
# construct adjacent equations
Z=(0,)*n
def cube(f):return mul(mul(f,f),f)
vals={(0,1):lin(1,{}),(0,2):lin(2,{})}
vals[(0,3)]=lin(-3,{0:1,2:1,3:1});vals[(0,4)]=lin(-2,{0:-1,1:2,2:1,3:1})
vals[(1,2)]=lin(-4,{1:1,2:1,3:2});vals[(1,3)]=lin(1,{0:-1,1:1,2:1})
vals[(1,4)]=lin(0,{0:1});vals[(2,3)]=lin(0,{1:1});vals[(2,4)]=lin(0,{2:1});vals[(3,4)]=lin(0,{3:1})
edges=list(combinations(range(5),2)); F=[]
for v in range(1,5):F.append(add(*[cube(vals[e]) for e in edges if v in e],scale(add(*[cube(vals[e]) for e in edges if 0 in e]),-1)))
G=[monic(f) for f in F]
pairs=deque((i,j) for i in range(len(G)) for j in range(i))
steps=0
while pairs:
 i,j=pairs.popleft(); h=reduce_poly(spol(G[i],G[j]),G)
 if h:
  h=monic(h); k=len(G);G.append(h)
  for i in range(k):pairs.append((k,i))
 steps+=1
 if steps%1000==0:print('step',steps,'basis',len(G),'queue',len(pairs),flush=True)
 if steps>50000:print('stop');break
print('done steps',steps,'basis',len(G))
assert steps == 8128 and len(G) == 128

# Interreduce the completed modular basis: retain minimal leading monomials,
# then reduce every non-leading term by the other retained polynomials.
minimal=[]
for i,g in enumerate(G):
    if not any(j!=i and divisible(lm(g),lm(h)) for j,h in enumerate(G)):
        minimal.append(g)
R=[]
for i,g in enumerate(minimal):
    h=reduce_poly(g,[q for j,q in enumerate(minimal) if j!=i])
    R.append(monic(h))
# A second pass removes any accidental leading divisibility after reduction.
R=[g for i,g in enumerate(R) if not any(j!=i and divisible(lm(g),lm(h)) for j,h in enumerate(R))]
R=[monic(reduce_poly(g,[q for j,q in enumerate(R) if j!=i])) for i,g in enumerate(R)]
print('interreduced basis',len(R),'lms',sorted(lm(g) for g in R))
assert len(R) == 16
# Count standard monomials by closure under multiplication by variables.
standard={(0,)*n}; frontier=[(0,)*n]
while frontier:
    m=frontier.pop()
    for c in range(n):
        q=list(m);q[c]+=1;q=tuple(q)
        if any(divisible(q,lm(g)) for g in R): continue
        if q not in standard: standard.add(q);frontier.append(q)
print('standard monomials',len(standard))
assert len(standard) == 54
# reduce product using reduced basis
prod={(0,)*n:1}
for i,j in combinations(range(10),2):
 q=add(mul(vals[edges[i]],vals[edges[i]]),scale(mul(vals[edges[j]],vals[edges[j]]),-1))
 prod=reduce_poly(mul(prod,q),R)
 if not prod: print('P9 zero after factor',i,j);break
print('P9 remainder terms',len(prod))
assert not prod
