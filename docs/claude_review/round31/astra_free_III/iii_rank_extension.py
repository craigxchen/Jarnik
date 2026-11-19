"""Exact check of the projected row information before cofactor errors."""
from itertools import combinations
from collections import Counter
from ii_paley import double_core


def rank_mod(rows, p=1009):
    pivots = {}
    for source in rows:
        row=[t%p for t in source]
        for i in sorted(pivots):
            if row[i]:
                a = row[i]
                row = [(b-a*c)%p for b,c in zip(row, pivots[i])]
        pivot = next((i for i,t in enumerate(row) if t),None)
        if pivot is not None:
            a = pow(row[pivot],-1,p)
            pivots[pivot] = [a*t%p for t in row]
    return len(pivots)


def run(q):
    cols = [c[1:] for family in double_core(q) for c in family]
    ranks=Counter()
    for i,j,k in combinations(range(q),3):
        mask = [c for c in cols if len({c[i],c[j],c[k]}) == 2]
        vertices = [u for u in range(q) if u not in (i,j,k)]
        for c in mask:
            l = i if c[i]!=c[k] else j
            for x in vertices:
                assert min(c[x],c[l])-min(c[x],c[k]) == c[x]*(c[l]-c[k])
        base=vertices[0]
        ranks[rank_mod([[c[x]-c[base] for c in mask] for x in vertices[1:]])]+=1
    print('ordinary marks, projected affine ranks across anchors:',q,sorted(ranks.items()))
    assert set(ranks) == {q-4}  # The number of independent differences is at most q-4.


if __name__=='__main__':
    for q in (7,11,19):
        run(q)
