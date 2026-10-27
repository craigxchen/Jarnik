# (1) vertex divisor bound: #{d | n : sqrt n <= d <= sqrt n + Y} <= floor(Y^2/sqrt n)+2 (via s = d + n/d integer)
# (2) near-diagonal circle count: #{(a,b): a^2+b^2=m, a,b in [N0,N0+H]} <= 2*(floor(H^2/(sqrt2*N0))+2)
import math, random
from sympy import divisors, factorint
random.seed(1)
worst1=0
# highly composite-ish n
cands=[]
for _ in range(4000):
    n=1
    for p in random.sample([2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53],random.randint(4,9)):
        n*=p**random.randint(1,3)
    if n<10**30: cands.append(n)
for n in cands:
    r=math.isqrt(n); r0=r if r*r==n else r+1
    ds=[d for d in divisors(n) if d>=r0]
    for Y in [int(n**0.25*c) for c in (0.3,1,2,4)]:
        cnt=sum(1 for d in ds if d<=r0+Y)
        bound=(Y+1)**2/math.isqrt(n)+2
        worst1=max(worst1,cnt-bound)
print('vertex divisor bound violations (max cnt-bound, must be <=0):',worst1)
# window n^{1/4-delta}: count <=2
mx=0
for n in cands:
    r=math.isqrt(n); r0=r if r*r==n else r+1
    Y=int(n**(0.25-0.05))
    cnt=sum(1 for d in divisors(n) if r0<=d<=r0+Y)
    mx=max(mx,cnt)
print('max divisors in [sqrt n, sqrt n + n^{0.20}] over',len(cands),'smooth n:',mx)
# (2) near-diagonal
def reps_box(N0,H):
    from collections import defaultdict
    D=defaultdict(int)
    for a in range(N0,N0+H+1):
        for b in range(N0,N0+H+1):
            D[a*a+b*b]+=1
    return max(D.values())
for N0,H in [(10**4,150),(10**4,400),(10**5,700),(2*10**5,1500)]:
    mxr=reps_box(N0,H)
    print('N0=%d H=%d gamma=%.3f max reps on one circle in box: %d ; sagitta bound 2*(H^2/(sqrt2 N0)+2)=%.1f'%(N0,H,math.log(H)/math.log(N0),mxr,2*(H*H/(math.sqrt(2)*N0)+2)))
# (3) exponent map: general-position arc R^lam  <->  near-diagonal window N0^gamma with gamma=(1+lam)/2
for g in [0.5,(5**0.5-1)/2,2/3,0.75]:
    print('gamma=%.4f -> general-position arc exponent lambda=2gamma-1=%.4f'%(g,2*g-1))
