# Range cost h'(v) = sum_P |sigma_P(v)| of a zonotope generator v in sdc_q (exact integers).
# Proved lower bound (survive.md Prop 4.1): h'(v) >= (q+1)||v||_2^2/(2||v||_1) >= (q+1)/2.
import numpy as np, itertools, random, sys
from slib import chi
def tables(q):
    Dp=np.zeros(q,dtype=np.int64); Dm=np.zeros(q,dtype=np.int64)
    Dp[0]=1; Dm[0]=1
    for t in range(1,q):
        if chi(t,q)==1: Dp[t]=1
        else: Dm[t]=1
    return Dp,Dm
def hprime(q,v,Dp,Dm):
    # S_a = sum_{t in D'} v_{a+t}, D'={0}uQR ; S^-_a = sum_{t in {0}uNQR} v_{a+t}
    S=np.zeros(q,dtype=np.int64); Sm=np.zeros(q,dtype=np.int64)
    a=np.arange(q)
    for x,val in v.items():
        S+=val*Dp[(x-a)%q]; Sm+=val*Dm[(x-a)%q]
    return int(np.abs(S).sum()+np.abs(Sm).sum()), S, Sm
def check_parseval(q,v,S,Sm):
    K=sum(v.values()); n2=sum(x*x for x in v.values())
    return int((S*S).sum())*4==(q+1)*(K*K+n2) and int((Sm*Sm).sum())*4==(q+1)*(K*K+n2)
for q in [int(a) for a in sys.argv[1:]] or [3251]:
    Dp,Dm=tables(q); best={}
    # supports of size 1..3 up to AGL (first two points 0,1), entries in {+-1,+-2}
    vals=[1,-1,2,-2]
    cands=[{0:1},{0:2}]
    for b in vals: cands.append({0:1,1:b})
    for y in range(2,q):
        for b,c in itertools.product(vals,vals): cands.append({0:1,1:b,y:c})
    random.seed(1)
    for s in range(4,9):
        for _ in range(3000):
            pts=random.sample(range(q),s); cands.append({x:random.choice(vals) for x in pts})
    allpar=True; mins={}
    for v in cands:
        h,S,Sm=hprime(q,v,Dp,Dm)
        allpar&=check_parseval(q,v,S,Sm)
        s=len(v)
        if s not in mins or h<mins[s][0]: mins[s]=(h,dict(v))
    print(f"q={q}: Parseval identity sum S_a^2 = (q+1)(K^2+||v||^2)/4 holds on all {len(cands)} candidates: {allpar}")
    for s in sorted(mins):
        h,v=mins[s]
        print(f"  support {s}: min h'(v) = {h} = {h/(q+1):.4f}(q+1)  at v={v if s<=3 else '(random sample)'}")
