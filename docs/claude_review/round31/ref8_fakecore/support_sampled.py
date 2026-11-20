# Sampled and structured balanced primitive v with larger ||v||_1; min support of Kv mod l (odd l<=h) and mod 4.
import numpy as np, math, random, sys
from math import gcd
from functools import reduce
from sympy import primerange, primitive_root
def K_of(q):
    sq=set((i*i)%q for i in range(1,q)); c=np.array([0]+[1 if i in sq else -1 for i in range(1,q)])
    idx=np.arange(q); return c[(idx[None,:]-idx[:,None])%q]
rng=random.Random(1)
def rand_v(q,h):
    v=[0]*q
    half=h//2
    for _ in range(half): v[rng.randrange(q)]+=1
    for _ in range(half): v[rng.randrange(q)]-=1
    return v
def check(K,q,v,stats):
    v=np.array(v); h=int(np.abs(v).sum())
    if h==0 or v.sum()!=0 or reduce(gcd,[int(x) for x in v])!=1: return
    kv=K@v; bound=math.sqrt(q)-h-1
    for l in [4]+[int(p) for p in primerange(3,h+1)]:
        s=int((kv%l!=0).sum())
        key=(h,l); stats[key]=min(stats.get(key,10**9),s)
        if s<bound: print("VIOLATION",q,h,l,s,v.nonzero())
for q in [int(x) for x in sys.argv[1:]]:
    K=K_of(q); stats={}
    for h in [6,8,10,12,16,20]:
        for _ in range(3000): check(K,q,rand_v(q,h),stats)
    # structured: v = 1_H - 1_{aH} for multiplicative subgroups H, and translates/arith progressions
    g=primitive_root(q)
    for k in [d for d in range(2,q) if (q-1)%d==0 and d<=12]:
        H=[pow(g,(q-1)//k*t,q) for t in range(k)]
        for a in [g, 2, 3]:
            v=[0]*q
            for x in H: v[x]+=1
            for x in H: v[a*x%q]-=1
            check(K,q,v,stats)
    for k in range(2,11):  # arithmetic progressions: 1_{[0,k)} - 1_{[b,b+k)}
        for b in range(1,q):
            v=[0]*q
            for x in range(k): v[x]+=1; v[(b+x)%q]-=1
            check(K,q,v,stats)
    for k in range(1,6):  # squares interval stuff: 1_{x^2: x<=k} - shift
        for b in range(1,q,7):
            v=[0]*q
            for x in range(1,k+1): v[x*x%q]+=1; v[(x*x+b)%q]-=1
            check(K,q,v,stats)
    print("q=%d sqrt=%.2f"%(q,math.sqrt(q)))
    for key in sorted(stats): print("  h=%d mod %d: min supp %d (bound %.2f)"%(key[0],key[1],stats[key],math.sqrt(q)-key[0]-1))
