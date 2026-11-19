"""Finite certificates for the squarefree eight-root divisors in Part II.

First compile and run ii_min_flat.cpp (see the report), producing
ii_min_cases.txt. All remaining checks use exact rational arithmetic.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations

import sympy as sp

from ii_min_class import coordinates, read, true_count
from ii_min_plane import interpolation
from exact_plane import plane_certificate


def verify_ternary():
    groups, exterior = read()
    assert len(groups) == 3 and len(exterior) == 2
    reps = [g[0][1] for g in groups.values()]
    first = [0, 5, 12, 14, 16]
    second = [0, 5, 12, 14, 18]
    last = [2, 7, 14, 23, 27]
    assert reps == [first, second, last]
    for ids in (first, second):
        h = coordinates(ids)
        reason, pair = plane_certificate(h)
        assert reason == 'vanishing distance product'
        print('Planar', ids, reason, pair)
    h, heights, equations = interpolation(last)
    a, b = sp.symbols('a b')
    assert equations[0] == 12*(a-b)**2*(b-1)*(2*a+b-1)
    assert sp.expand(equations[1] - 12*(b-1)*(2*a+b-1)*(2*a*a-2*a*b-2*a+b*b+1)) == 0
    assert h[3, :] == 2*h[6, :]-h[7, :]
    assert [heights[i].subs(b, 1) for i in (3,6,7)] == [1,1,1]
    print('Third planar case forces coincident roots by a constant-y midpoint.')
    cube = {(a, b, c, 1-a-b-c) for a in (0,1) for b in (0,1) for c in (0,1)}
    for n, ids, h in exterior:
        assert set(map(tuple, h.tolist())) == cube and true_count(h, False) == n == 12
    # Distinct minimum support vectors in the known eight-row extremal code.
    code = [15,51,60,85,170,195,204,240]
    minimum = Counter(i ^ j for i,j in combinations(code,2) if (i ^ j).bit_count() == 4)
    assert len(minimum) == 12 and set(minimum.values()) == {2}
    print('Cube code:',len(minimum),'four-supports, each occurring twice as an edge.')
    q = [Fraction(-9,4), Fraction(-2,3), Fraction(1,4), Fraction(8,3)]
    points = []
    for h, e in zip(q, (1,-1,-1,1)):
        for x in (e/h, -e/h):
            y = abs(h)
            points.append((y,x*y,3*x*x*y-y**3))
    assert true_count(sp.Matrix(points), False) == 12
    print('Explicit rational template: 12 primitive ternary directions up to sign.')


def verify_gap():
    groups, exterior = read('ii_min_gap.txt', threshold=11)
    assert len(exterior) == 5 and all(true_count(h, False) == 12 for _, _, h in exterior)
    eleven = 0
    for group in groups.values():
        h = coordinates(group[0][1])
        count = true_count(h, False)
        if count == 11:
            reason, pair = plane_certificate(h)
            print('Eleven-direction class', group[0][1], reason, pair)
            assert reason == 'vanishing distance product'
            eleven += 1
        else:
            assert count >= 12
    assert eleven == 2


if __name__ == '__main__':
    import sys
    if '--gap' in sys.argv:
        verify_gap()
    else:
        verify_ternary()
