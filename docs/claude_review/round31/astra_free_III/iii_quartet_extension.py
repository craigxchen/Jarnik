"""All matching partitions at a parent with a planted passing quartet."""
from collections import Counter
from itertools import combinations, permutations, product
from math import gcd, isqrt, prod

from ii_paley import double_core
from iii_quartet import modular_minimum, quartet_types
from iii_quartet_exhaust import small_relation


PLANTED = (
    (1000000873621, 1000000976641, 1000001041141),
    (1000000874941, 1000000973341, 1000001043121),
    (1000000876981, 1000000968661, 1000001045761),
)
EXTRA = (1000000963813, 1000000963837, 1000000963877, 1000000963937, 1000000964021)
NAMES = tuple(permutations(range(3)))
ALL = (1<<14)-1


def fingerprint(left, right):
    """Orbit of the 4 by 4 intersection table under renaming the two triples."""
    def orders(triple):
        complement = ALL ^ (triple[0]|triple[1]|triple[2])
        return [tuple(triple[i] for i in perm)+(complement,) for perm in NAMES]
    return min(tuple((a&b).bit_count() for a,b in product(x,y))
               for x,y in product(orders(left),orders(right)))


def main():
    ps=sorted([p for triple in PLANTED for p in triple]+list(EXTRA))
    assert len(set(ps))==14
    assert all(p%4==1 and all(p%q for q in range(3,isqrt(p)+1,2)) for p in ps)
    denom=16*isqrt(prod(ps))
    u=(prod(ps[-8:])+denom-1)//denom
    bound=isqrt(u)
    planted_products=list(map(prod,PLANTED))
    a,b,c=planted_products
    coeff=(b-c,c-a,a-b)
    coeff=tuple(v//gcd(*coeff) for v in coeff)
    assert sum(x*y for x,y in zip(coeff,planted_products))==0
    for v in coeff:
        # Each of the three terms splits into two permitted even cofactors after scaling by 4.
        assert any(abs(v)%d==0 and abs(v)//d<=bound//2 for d in range(1,min(isqrt(abs(v)),bound//2)+1))
    print('U, edge bound, planted primitive coefficient triple:',u,bound,coeff,flush=True)
    groups=[(sum(1<<i for i in s),prod(ps[i] for i in s)) for s in combinations(range(14),3)]
    passing=[]
    tested=0
    for i,(ma,a) in enumerate(groups):
        for j in range(i+1,len(groups)):
            mb,b=groups[j]
            if ma&mb:
                continue
            used=ma|mb
            for mc,c in groups[j+1:]:
                if used&mc:
                    continue
                tested+=1
                if modular_minimum(a,b,c)<=2*u*u and small_relation(a,b,c,u) is not None:
                    passing.append((ma,mb,mc))
    assert tested==560560 and len(passing)==1
    print('total disjoint triples of products, passing:',tested,len(passing),flush=True)
    cols=[c[1:] for family in double_core(7) for c in family]
    five_types=Counter()
    for five in combinations(range(7),5):
        partitions=[]
        for four in combinations(five,4):
            classes=quartet_types(cols,four)
            if len(classes[0])==3:
                partitions.append(tuple(sum(1<<i for i in s) for s in classes))
        assert len(partitions)==2
        five_types[fingerprint(*partitions)]+=1
    compatible=sum(fingerprint(x,y) in five_types for x,y in product(passing, repeat=2))
    print('five-mark restrictions by overlap fingerprint:',sorted(five_types.values()),flush=True)
    print('ordered pairs of passing prime partitions compatible with any five marks:',compatible,flush=True)
    assert compatible==0


if __name__=='__main__':
    main()
