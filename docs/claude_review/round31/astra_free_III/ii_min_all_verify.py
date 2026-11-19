"""All integral vectors of l1 norm four, including repeated factors.

Run ./ii_min_flat 12 all > ii_min_allcases.txt before this verifier.
The ternary-only, faster check is ii_min_verify.py.
"""
from collections import Counter
from itertools import permutations

import sympy as sp

from exact_plane import plane_certificate
from ii_min_class import coordinates, read, true_count
from ii_min_plane import a, b, interpolation


def equal(x, y):
    assert sp.expand(x-y) == 0


def midpoint_same_height(h, heights, indices, substitution):
    i, j, k = indices
    assert h.row(i)+h.row(k) == 2*h.row(j)
    for pair in [(i,j),(j,k)]:
        equal((heights[pair[0]] - heights[pair[1]]).subs(substitution), 0)


def triangle_test(ids):
    h, heights, eq = interpolation(ids)
    if ids == [0,8,19,23,33]:
        # 2a+b-1=y_4+b>0. Unless a=1, the two nonnegative squared
        # factors force a=b and a+b=2, a contradiction.
        equal(eq[0], 3*(a-1)*(a-b)**2*(2*a+b-1))
        equal(eq[1], 3*(a-1)*(a+b-2)**2*(2*a+b-1))
        equal(2*a+b-1, heights[4]+b)
        midpoint_same_height(h, heights, (4,5,7), {a:1})
    elif ids == [6,13,20,24,27]:
        # a+b-1>0 since 2(a+b-1)=y_0+1. The indicated midpoint
        # collapses for each of the remaining two branches.
        equal(eq[1], 48*(a-1)**2*(b-1)*(a+b-1))
        equal(2*(a+b-1), heights[0]+1)
        midpoint_same_height(h, heights, (3,5,7), {a:1})
        midpoint_same_height(h, heights, (4,6,7), {b:1})
    else:
        assert ids == [3,14,21,29,35]
        # a-5b+2=-(y_0+1)<0. On either possible branch b=1:
        # 2b-1=y_2>0 and 3b-1=b+y_2>0.
        equal(eq[0], -12*(a-b)*(a-5*b+2)*(a-2*b+1)**2)
        equal(a-5*b+2, -heights[0]-1)
        equal(2*b-1, heights[2])
        equal(eq[1].subs(a,b), 48*(b-1)**3*(2*b-1))
        equal(eq[1].subs(a,2*b-1), 12*(b-1)**3*(3*b-1))
        midpoint_same_height(h, heights, (2,6,7), {b:1})


def centered_hexagon(h):
    points = {tuple(row) for row in h.tolist()}
    def add(a,b):return tuple(x+y for x,y in zip(a,b))
    def sub(a,b):return tuple(x-y for x,y in zip(a,b))
    for center in points:
        displacements = {sub(p,center) for p in points if p != center}
        for u, v in permutations(displacements, 2):
            w = add(u,v)
            triple = {u,v,w}
            if len(triple) != 3:continue
            with_negatives = triple | {tuple(-x for x in t) for t in triple}
            if len(with_negatives) == 6 and with_negatives <= displacements:
                return True
    return False


if __name__ == '__main__':
    groups, exterior = read('ii_min_allcases.txt')
    assert len(groups) == 59 and len(exterior) == 10
    special = [[0,8,19,23,33], [6,13,20,24,27], [3,14,21,29,35]]
    counts = Counter()
    for group in groups.values():
        # Representatives with >=13 observed ordered patterns come first
        # when available, for agreement with the threshold-13 census.
        ids = next((ids for n,ids in group if n>=13), group[0][1])
        h = coordinates(ids)
        counts[true_count(h, True)] += 1
        if ids in special:
            triangle_test(ids)
            print('Factored equations exclude', ids)
        else:
            reason, _ = plane_certificate(h)
            assert reason in {'unit ideal','vanishing height product','vanishing distance product'}
    assert counts == {12:30,13:15,14:9,15:1,17:2,18:1,19:1}
    cube = {(x,y,z,1-x-y-z) for x in (0,1) for y in (0,1) for z in (0,1)}
    n_cube = 0
    for n, ids, h in exterior:
        assert n == true_count(h, True) == 12
        if set(map(tuple, h.tolist())) == cube:
            n_cube += 1
        else:
            assert centered_hexagon(h)
    assert n_cube == 2
    print('Ruled out all 59 planar shapes. Rank-four survivors equal the ternary cube.')
