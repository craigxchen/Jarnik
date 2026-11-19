"""Optional exact test of the one-anchor-pair valuation projection in the double core."""
from collections import Counter
from itertools import combinations
from sympy import Matrix
from ii_paley import double_core


def ranks(q):
    columns=sum(double_core(q),[])
    counts=Counter()
    for a,b in combinations(range(q+1),2):
        active=[c for c in columns if c[a]!=c[b]]
        others=[i for i in range(q+1) if i not in (a,b)]
        h=others[0]
        values=Matrix([[c[i]-c[h] for c in active] for i in others[1:]])
        counts[(len(active),values.rank())]+=1
    assert list(counts)==[(q+1,q-2)]
    print('m',q+1,'(active factors, affine rank): counts',dict(counts))


if __name__=='__main__':
    for q in (7,11):
        ranks(q)
