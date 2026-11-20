# Own enumeration: ordered pairs (P,N) of disjoint t-multisets over a label set L in F_q,
# counted by the number of initial power sums p_1..p_r that agree mod q.
# Also: vector-level count vs Astra's bound binom(|L|+t-1,t)*binom(|L|+t-r-1,t-r)
# (her k = r+1), and the recovery map (P, smallest t-r of N) -> N is injective.
from itertools import combinations_with_replacement as cwr
from math import comb
from collections import defaultdict
def run(q,t,labels):
    ms=list(cwr(labels,t))
    ps=[tuple(sum(pow(x,s,q) for x in a)%q for s in range(1,t+1)) for a in ms]
    sup=[set(a) for a in ms]
    cnt=[0]*(t+1); inj_ok=True
    for r in range(t+1):
        seen=defaultdict(set)
    keys=[defaultdict(set) for _ in range(t+1)]
    for i,a in enumerate(ms):
        for j,b in enumerate(ms):
            if sup[i]&sup[j]: continue
            r=0
            while r<t and ps[i][r]==ps[j][r]: r+=1
            for rr in range(r+1):
                cnt[rr]+=1
                key=(a,b[:t-rr])   # b is sorted: first t-rr smallest labels
                keys[rr][key].add(b)
    inj=all(len(v)==1 for k in keys for v in k.values())
    n=len(labels)
    bounds=[comb(n+t-1,t)*comb(n+t-r-1,t-r) for r in range(t+1)]
    return cnt,bounds,inj
print("F_11, all 11 labels, t=3:",run(11,3,list(range(11))))
print("F_11, labels 1..10, t=3:",run(11,3,list(range(1,11))))
for q,t in ((7,1),(7,2),(7,3),(7,4),(13,3),(13,4),(11,4)):
    c,b,inj=run(q,t,list(range(1,q)))
    assert all(x<=y for x,y in zip(c,b)) and inj and c[t]==0
    print(q,t,c,"bounds",b,"injective",inj)
