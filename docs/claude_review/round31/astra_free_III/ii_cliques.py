"""All subset-sum fibres and their minimum-shell cliques for the three examples."""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from ii_template import ROOTS, SELECT, directions, moments

SECOND = ((-23,4),(-11,10),(-8,16),(2,16),(5,10),(5,12),(13,6),(17,4))


def maximum_cliques(roots, eligible=None):
    h = [moments(root) for root in roots]
    buckets = defaultdict(list)
    for x in range(1 << len(roots)):
        signature = tuple(sum(h[i][j] for i in range(len(roots)) if x & (1 << i))
                          for j in range(3))
        buckets[signature].append(x)
    words = directions(roots)
    supports = {sum(1 << i for i, sign in words[k])
                for k in (range(len(words)) if eligible is None else eligible)}
    largest_fibre = max(map(len,buckets.values()))
    largest = largest_sidon = 0
    for bucket in buckets.values():
        for size in range(1,len(bucket)+1):
            for marks in combinations(bucket,size):
                labels = [a ^ b for a,b in combinations(marks,2)]
                if all(v in supports for v in labels):
                    largest = max(largest,size)
                    if len(set(labels)) == len(labels):
                        largest_sidon = max(largest_sidon,size)
    return largest_fibre,largest,largest_sidon


if __name__ == '__main__':
    fractions = (Fraction(-9,4),Fraction(-2,3),Fraction(1,4),Fraction(8,3))
    cube = [(x,abs(t)) for t,e in zip(fractions,(1,-1,-1,1)) for x in (e/t,-e/t)]
    for label, roots, chosen in [('cube',cube,None),('seven',ROOTS,SELECT),
                                 ('nine',SECOND,(0,1,2,3,5,7,8))]:
        print(label,maximum_cliques(roots), 'selected', maximum_cliques(roots,chosen))
