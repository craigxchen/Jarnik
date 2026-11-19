"""Check possibilities for four distinct sixth moment labels modulo the inert prime 7."""
from itertools import combinations, product
from collections import defaultdict, Counter
from tenpair_modular import power_moments

p=7
roots=[(x,y) for x in range(p) for y in range(1,p)]
opp=lambda r: (r[0],-r[1]%p)
M={r:power_moments(p,*r) for r in roots}
def sig(s):
    return tuple(sum(c)%p for c in zip(*(M[r] for r in s)))
A=defaultdict(list); B=defaultdict(list)
for s in combinations(roots,4):
    if not set(s).isdisjoint(opp(r) for r in s):continue
    A[sig(s)].append(s)
for b in roots:
    for c,d in combinations(roots,2):
        s=(b,b,c,d)
        if len(set(s))!=3 or not set(s).isdisjoint(opp(r) for r in s):continue
        B[sig(s)].append((b,frozenset((c,d))))
match=0
for k, rows in B.items():
    for triple in combinations(rows,3):
        bs=[b for b,cd in triple]
        if len(set(bs))<3:continue
        cds=[cd for b,cd in triple]
        if any(len(a&b)!=1 for a,b in combinations(cds,2)) or len(set.union(*(set(q) for q in cds)))!=3:continue
        if any(not set(bs).isdisjoint(cd) for cd in cds):continue
        for s in A[k]:
            if set(bs).issubset(s) and not ((set(s)-set(bs))&set.union(*(set(q) for q in cds))):
                match+=1
print("Possible four-label residue profiles (even before testing sixth moments):",match)
