# Lemma 2.1 checks: AGL(1,q) preserves the sdc_q pattern set (exact); l_s(tau k)=zeta^s l_s(k), l_1(mu_s k)=l_s(k)
# on random integer k (exact in Z[x]/(x^q-1)); distinctness of 2i-m mod q for m<q.
import random
from slib import patterns
def elt(k,q,s):  # l_s(k) as coefficient vector in Z[x]/(x^q-1)
    v=[0]*q
    for x in range(q): v[(s*x)%q]+=k[x]
    return v
ok=True
for q in [7,11,19,23,31,43,47,59,67,71,79,83,3251]:
    if q<=83:
        pats=set(patterns(q))
        aff=all(frozenset((s*x+b)%q for x in S) in pats for S in pats for s in range(1,q) for b in (0,1))
    else: aff=None
    random.seed(q); good=True
    for _ in range(20):
        k=[random.randint(-3,3) for _ in range(q)]
        tk=[k[(x-1)%q] for x in range(q)]
        s=random.randint(1,q-1); si=pow(s,-1,q)
        mk=[k[(si*x)%q] for x in range(q)]
        # l_s(tau k) = zeta^s l_s(k): shift coefficient vector by s
        a=elt(tk,q,s); b=elt(k,q,s); b=[b[(t-s)%q] for t in range(q)]
        good&= (a==b) and (elt(mk,q,1)==elt(k,q,s))
    dist=all(len({(2*i-m)%q for i in range(m+1)})==m+1 for m in range(1,q))
    ok&=good and dist and (aff is not False)
    print(f"q={q}: AGL(1,q) preserves pattern set: {aff}; translation/multiplier identities: {good}; 2i-m distinct mod q for all m<q: {dist}")
print("ALL OK" if ok else "FAILURE")
