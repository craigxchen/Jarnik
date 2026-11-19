"""Degree-one K5 residue test over F49 with a 4+1 partition of the sixth moments."""
from collections import Counter, defaultdict
from itertools import combinations_with_replacement
from tenpair_modular import power_moments

p=7
roots=[(x,y) for x in range(p) for y in range(1,p)]
def opp(r):
    return (r[0],(-r[1])%p)
imag6={}
data={}
for r in roots:
    x,y=r
    a,b=1,0
    for n in range(6):a,b=(a*x-b*y)%p,(a*y+b*x)%p
    data[r]=power_moments(p,x,y)
    imag6[r]=b

def sig(s):
    return tuple(sum(c)%p for c in zip(*(data[r] for r in s)))
def netok(s):
    return set(s).isdisjoint(opp(r) for r in s)

base=(0,1)
target=sig((base,)*3)
triples=[s for s in combinations_with_replacement(roots,3) if base not in s and netok(s+(base,)) and sig(s)==target and sum(imag6[r] for r in s)%p!=3*imag6[base]%p]
print('4+1 with no star cancellation, candidate triples',triples)
# Cancellation in each star. Find possible common residual two-root multisets u,v by signature.
pairs=defaultdict(list)
for uv in combinations_with_replacement(roots,2):
    if netok(uv):
        pairs[sig(uv)].append(uv)
cand=[]
for b23 in combinations_with_replacement(roots,3):
    bs=(base,)+b23
    if not netok(bs):continue
    mult=Counter(bs)
    for uv in pairs.get(sig(bs),[]):
        if not set(uv).isdisjoint(bs):continue
        if any(n%2 and opp(r) not in uv for r,n in mult.items()):continue
        if sum(imag6[r] for r in uv)%p == sum(imag6[r] for r in bs)%p:continue
        cand.append((bs,uv))
print('4+1 with star cancellation, candidate multisets',cand)
