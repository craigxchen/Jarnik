# -K-degree of Huang's "general lines" (p.54) through e on the permutohedral fan (split, over Q(i)).
# Relation: sum_{rho in sigma(1)} rho + sum_{rho in tau(1)} c_rho rho = 0, tau the cone containing -sum in its relative interior.
import itertools
for k in (3,4,5):
    s=list(range(k)); cone=[frozenset(s[:i]) for i in range(1,k)]
    v=[0]*k
    for A in cone:
        for a in A: v[a]+=1
    w=[-x for x in v]; m=min(w); w=[x-m for x in w]   # -sum, normalised representative in Z^k/Z1
    # decompose w (nonneg integer vector, min 0) as sum of 1_{level sets}: flag decomposition
    c=0; levels=sorted(set(w))
    for lo,hi in zip(levels,levels[1:]): c+=hi-lo
    print(f"k={k}: n={k-1}, -K.line = n + sum c = {k-1}+{c} = {k-1+c}; -sum = {w}")
