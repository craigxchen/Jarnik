# check: for a boundary valuation vector a = sum_k c_k 1_{A_k} (nested chain, min 0),
# ord J = sum_{i<j} min(a_i,a_j) - 2 sum a_i  equals  sum_k c_k (C(|A_k|,2) - 2|A_k|);
# pair-cut layers (|A|=2 or 3) give -3, singleton-cut layers (|A|=1 or 4) give -2.
import itertools, random
from math import comb
random.seed(0)
for trial in range(2000):
    perm=list(range(5)); random.shuffle(perm)
    sizes=sorted(random.sample(range(1,5),random.randint(1,4)),reverse=True)
    a=[0]*5; pred=0
    for s in sizes:
        c=random.randint(1,5); A=perm[:s]
        for i in A: a[i]+=c
        pred+=c*(comb(s,2)-2*s)
    val=sum(min(a[i],a[j]) for i,j in itertools.combinations(range(5),2))-2*sum(a)
    assert val==pred,(a,sizes)
print("ordJ layer formula OK; per-layer values:",{s:comb(s,2)-2*s for s in range(1,5)})
