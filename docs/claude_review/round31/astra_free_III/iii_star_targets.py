"""Exact gain tests over every real center for model stars and four-edge matchings."""
from collections import Counter
from fractions import Fraction
from itertools import product
from ii_paley import double_core


def optimal_mean_gain(rows,edges):
    # Equal model column weights. A coordinate median minimizes total L1 cost.
    total=0
    for col in zip(*rows):
        differences=sorted(col[i]-col[j] for i,j in edges)
        med=differences[len(edges)//2]
        total+=sum(abs(d)-abs(d-med) for d in differences)
    return Fraction(total,len(edges)*len(rows[0]))


def perfect_matchings(vertices):
    if not vertices:
        yield ()
        return
    first,*tail=vertices
    for index,second in enumerate(tail):
        for pairs in perfect_matchings(tail[:index]+tail[index+1:]):
            yield ((first,second),)+pairs


def check(q):
    cols=sum(double_core(q),[])
    rows=list(zip(*cols)); m=len(rows)
    star=Fraction(1,m-1)
    for i in range(m):
        assert optimal_mean_gain(rows,[(i,j) for j in range(m) if j!=i])==star
    print('m',m,'each full outgoing star, maximum normalized mean gain',star)
    if m==8:
        hist=Counter()
        for pairs in perfect_matchings(list(range(8))):
            for orientations in product((0,1),repeat=4):
                edges=[(p[t],p[1-t]) for p,t in zip(pairs,orientations)]
                hist[optimal_mean_gain(rows,edges)]+=1
        assert sum(hist.values())==1680
        print('all 1680 oriented matchings, histogram of optimum normalized gain',
              sorted(hist.items()))


if __name__=='__main__':
    for q in (7,11,19):
        check(q)
