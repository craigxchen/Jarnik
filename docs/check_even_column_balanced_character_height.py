"""Exact balanced-character averaging and literal Gaussian identity checks."""
import json
from collections import Counter
from itertools import product
from math import gcd
from pathlib import Path
from random import Random

L=[v for v in product((-1,1),repeat=8) if sum(v)==0 and v[0]==1]
assert len(L)==35
for col in product((-1,1),repeat=8):
 h=sum(x==1 for x in col)
 if h not in (2,4,6):continue
 nums=[sum(a*b for a,b in zip(lam,col)) for lam in L]
 assert all(n%4==0 for n in nums)
 values=[abs(n)//4 for n in nums]
 assert Counter(values)==({0:20,1:15} if h in (2,6) else {0:18,1:16,2:1})
print('PASS: all 126 even nonconstant columns have exact sums 15 (pair) or 18 (balanced).')

f=json.loads(Path(__file__).with_name('paley24_even_column_moment_fixture.json').read_text())
H=f['hadamard'];count=0
for rows in f['positive_row_sets']:
 T=[H[i] for i in rows];cols=list(zip(*T))
 assert all(sum(a*b for a,b in zip(T[i],T[j]))==(24 if i==j else 0) for i in range(8) for j in range(8))
 assert cols[0]==(1,)*8
 assert Counter(min(col.count(1),col.count(-1)) for col in cols)=={0:1,2:8,4:15}
 for deleted,col in enumerate(cols):
  if sum(col):continue
  count+=1
  # SS^t=24I-11^t-bb^t. Since 1 and b are orthogonal and
  # both have squared length 8, its eigenvalues are 16,16,24,...,24.
  assert sum(col)==0 and sum(x*x for x in col)==8
  for lam in L:
   c1=sum(lam);cb=sum(a*b for a,b in zip(lam,col))
   energy=24*sum(x*x for x in lam)-c1*c1-cb*cb
   assert energy>=16*sum(x*x for x in lam)>0
assert count==11385
print('PASS: all 11385 positive deleted codes have positive Gram rank and no zero balanced character.')

def mul(z,w):return (z[0]*w[0]-z[1]*w[1],z[0]*w[1]+z[1]*w[0])
def conj(z):return (z[0],-z[1])
def norm(z):return z[0]*z[0]+z[1]*z[1]
def gp(z,n):
 r=(1,0)
 for _ in range(n):r=mul(r,z)
 return r
def isprime(n):return n>=2 and all(n%d for d in range(2,int(n**.5)+1))
blocks=[];seen=set()
for x in range(1,40):
 for y in range(1,x):
  n=x*x+y*y
  if isprime(n) and n not in seen:blocks.append((x,y));seen.add(n)
  if len(blocks)==22:break
 if len(blocks)==22:break
assert len(blocks)==22
rep=f['representative'];S=[[H[i][j] for j in rep['retained_columns']] for i in rep['rows']]
units=[(1,0),(0,1),(-1,0),(0,-1)]
for odd in (False,True):
 exps=[i%4 for i in range(8)];exps[0]+=int(odd)
 eps=[units[k%4] for k in exps]
 Z=[]
 for row,u in zip(S,eps):
  z=mul((2,1),u)
  for s,h in zip(row,blocks):z=mul(z,h if s==1 else conj(h))
  Z.append(z)
 assert len({norm(z) for z in Z})==1 and len(set(Z))==8
 for lam in L:
  nums=[sum(lam[i]*S[i][j] for i in range(8)) for j in range(22)]
  assert all(n%4==0 for n in nums)
  v=[n//4 for n in nums];assert any(v)
  beta=(1,0);normheight=1
  for e,h in zip(v,blocks):
   beta=mul(beta,gp(h if e>=0 else conj(h),abs(e)))
   normheight*=norm(h)**abs(e)
  assert norm(beta)==normheight and gcd(*beta)==1 and norm(beta)>1 and norm(beta)%2==1
  plus=minus=(1,0)
  for l,z in zip(lam,Z):
   if l==1:plus=mul(plus,z)
   else:minus=mul(minus,z)
  k=sum(l*e for l,e in zip(lam,exps))%4;u=units[k]
  assert k%2==int(odd)
  assert mul(plus,gp(conj(beta),2))==mul(mul(u,minus),gp(beta,2))
print('PASS: all 35 literal Gaussian identities, exact heights, nonvanishing and both unit-parity branches.')

# Repeated threshold blocks are not independent.  Compare their signed
# products, after cancellation, with the direct prime-allocation formula.
rng = Random(20260921)
cases = cancellation_cases = 0
for trial in range(24):
 allocations = []
 for _ in range(7):
  a = [0, 0, 1, 1, 3, 3, 6, 6]
  rng.shuffle(a)
  allocations.append(a)
 eps = [units[rng.randrange(4)] for _ in range(8)]
 Z = []
 for i in range(8):
  z = mul((2, 1), eps[i])
  for h, a in zip(blocks, allocations):
   z = mul(z, mul(gp(h, a[i]), gp(conj(h), 6-a[i])))
  Z.append(z)
 assert len({norm(z) for z in Z}) == 1
 for lam in L:
  beta = raw = (1, 0)
  height_bound = 1
  for h, a in zip(blocks, allocations):
   numerator = sum(lam[i]*a[i] for i in range(8))
   assert numerator % 2 == 0
   e = numerator // 2
   beta = mul(beta, gp(h if e >= 0 else conj(h), abs(e)))
   layer_sum = 0
   for k in range(1, 7):
    col = [1 if x >= k else -1 for x in a]
    assert sum(x == 1 for x in col) in (2, 4, 6)
    n = sum(lam[i]*col[i] for i in range(8))
    assert n % 4 == 0
    v = n // 4
    layer_sum += v
    raw = mul(raw, gp(h if v >= 0 else conj(h), abs(v)))
    height_bound *= norm(h)**abs(v)
   assert layer_sum == e
  content = gcd(abs(raw[0]), abs(raw[1]))
  assert content > 0 and raw == (content*beta[0], content*beta[1])
  assert gcd(*beta) == 1 and norm(beta) % 2 == 1
  assert norm(raw) == height_bound and norm(beta) <= height_bound
  cancellation_cases += content > 1
  plus = minus = u = (1, 0)
  for l, z, epsilon in zip(lam, Z, eps):
   if l == 1:
    plus = mul(plus, z)
    u = mul(u, epsilon)
   else:
    minus = mul(minus, z)
    u = mul(u, conj(epsilon))
  assert mul(plus, gp(conj(beta), 2)) == mul(mul(u, minus), gp(beta, 2))
  cases += 1
assert cases == 840 and cancellation_cases > 0
print(f'PASS: {cases} nested-layer identities, including {cancellation_cases} nontrivial content cancellations.')
