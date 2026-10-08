"""Exact all-integer-character barrier for H12 x H2, columns 0,1 deleted."""
from itertools import product
from collections import Counter

def h12():
 q=11
 chi=lambda x: 0 if x%q==0 else (1 if pow(x%q,5,q)==1 else -1)
 return [[1]*12]+[[1]+[-1 if a==b else chi(a-b) for b in range(q)] for a in range(q)]
H=h12(); rr=(0,1,2,4)
S=[[H[i][j]*(1 if e*d==0 else -1) for j in range(1,12) for d in range(2)] for i in rr for e in range(2)]
cols=list(zip(*S))
groups=[(0 if abs(sum(col))==4 else 1) if k%2==0 else (2 if abs(sum(col[::2]))==2 else 3) for k,col in enumerate(cols)]
w=[(5,4,6,4)[g] for g in groups]
assert Counter(groups)=={0:8,1:3,2:8,3:3}
assert sum(w)==112 and sum(w[j] for j,g in enumerate(groups) if g==0)==40
heavy=[j for j,g in enumerate(groups) if g in (0,2)]
# The sixteen heavy columns are a Hadamard order eight with every column doubled,
# allowing column sign reversals. T T^t=16 I supplies the inverse norm bound.
assert all(sum(S[i][j]*S[k][j] for j in heavy)==(16 if i==k else 0) for i in range(8) for k in range(8))
mins={}; count=0; slopes=set()
for c in product((-1,0,1),repeat=8):
 if sum(c) or not any(c):continue
 count+=1
 nums=[sum(c[i]*S[i][j] for i in range(8)) for j in range(22)]
 assert all(x%2==0 for x in nums)
 v=[abs(x)//2 for x in nums]
 V=sum(x*y for x,y in zip(v,w)); h=sum(abs(x) for x in c)
 mins[h]=min(mins.get(h,V),V)
 # Parameter weights (a,t,a+t/4,t): V=A*a+B*t.
 A=sum(v[j] for j in heavy)
 B4=4*sum(v[j] for j,g in enumerate(groups) if g in (1,3))+sum(v[j] for j,g in enumerate(groups) if g==2)
 slopes.add((A,B4))
 assert A>=8 and B4>=16
 assert V>=56
# For max|c|>=2, Hadamard inversion gives sum_heavy |cS/2| >= 8 max|c|.
# All heavy weights >=5, so V>=40 max|c|>=80>56. This exhausts every integer c.
print('nonzero ternary zero-sum characters:',count)
print('minimum V by support:',sorted(mins.items()))
print('parameter coefficient pairs (A,4B):',sorted(slopes))
print('PASS: all nonzero integer zero-sum c have V>=56; W=112, W_pair=40, W_bal=72.')

# Gaussian arithmetic and actual independent prime blocks.
def mul(z,w):return (z[0]*w[0]-z[1]*w[1],z[0]*w[1]+z[1]*w[0])
def conj(z):return (z[0],-z[1])
def norm(z):return z[0]*z[0]+z[1]*z[1]
def associates(z):return {z,(-z[0],-z[1]),(-z[1],z[0]),(z[1],-z[0])}
def orbit(z):return associates(z)|associates(conj(z))
def isprime(n):return n>=2 and all(n%d for d in range(2,int(n**.5)+1))
blocks=[]; norms=set()
for x in range(1,40):
 for y in range(1,x):
  n=x*x+y*y
  if isprime(n) and n not in norms:
   blocks.append((x,y));norms.add(n)
  if len(blocks)==24:break
 if len(blocks)==24:break
assert len(blocks)==24 and len(norms)==24
T=[[H[i][j]*(1 if e*d==0 else -1) for j in range(12) for d in range(2)] for i in rr for e in range(2)]
balanced=[k for k in range(24) if sum(row[k] for row in T)==0]
assert len(balanced)==15
for deleted in balanced:
 keep=[k for k in range(24) if k not in (0,deleted)]
 assert Counter(min(sum(T[i][k]==1 for i in range(8)),sum(T[i][k]==-1 for i in range(8))) for k in keep)=={2:8,4:14}
 anti=[k for k in keep if k%2]
 base=[k for k in keep if not k%2]
 Avec=[];Bvec=[]
 for ii in range(4):
  A=B=(1,0)
  for k in anti:A=mul(A,blocks[k] if T[2*ii][k]==1 else conj(blocks[k]))
  for k in base:B=mul(B,blocks[k] if T[2*ii][k]==1 else conj(blocks[k]))
  Avec.append(A);Bvec.append(B)
  direct=[]
  for e in range(2):
   z=(1,0)
   for k in keep:z=mul(z,blocks[k] if T[2*ii+e][k]==1 else conj(blocks[k]))
   direct.append(z)
  assert direct==[mul(B,A),mul(B,conj(A))]
 assert len({norm(A) for A in Avec})==len({norm(B) for B in Bvec})==1
 assert all(Avec[i] not in orbit(Avec[j]) for i in range(4) for j in range(i))
print('PASS: all 15 balanced deletions have the profile and actual conjugation-injective paired products.')

# Exact projection-level checks on actual integer circles. Each representative
# is chosen modulo units and conjugation; no endpoint assertion is made.
from math import atan2,sqrt,cos,pi,floor
for N in (5,25,65,85,325,1105,5525,32045):
 points=[]
 for x in range(1,int(sqrt(N))+1):
  for y in range(x+1):
   if x*x+y*y==N:points.append((x,y))
 assert all(points[i] not in orbit(points[j]) for i in range(len(points)) for j in range(i))
 classes=[[],[]];errors=[]
 for x,y in points:
  t=atan2(y,x);diagonal=int(t>pi/8)
  classes[diagonal].append((x,y));errors.append(abs(t-diagonal*pi/4))
 for kind,ps in enumerate(classes):
  levels=[x if kind==0 else x+y for x,y in ps]
  assert len(levels)==len(set(levels))
  if kind:assert all(v%2==N%2 for v in levels)
 if points:
  Delta=2*max(errors)+1e-9
  s=sqrt(N)*(1-cos(Delta/2))
  assert len(points)<=2+floor(s+1e-7)+floor(s/sqrt(2)+1e-7)
print('PASS: actual axis/diagonal levels are injective, have the asserted spacing, and satisfy capacity.')
