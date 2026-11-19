"""Exact affine classification for the pair-sum flat enumeration."""
from collections import defaultdict
from functools import reduce
from itertools import combinations, permutations
from math import gcd
from pathlib import Path

import sympy as sp

QUADRUPLES = list(combinations(range(8), 4))
TRIPLES = list(combinations(range(8), 3))


def coordinates(ids):
    matrix = sp.zeros(len(ids), 8)
    for r, idx in enumerate(ids):
        if idx < len(QUADRUPLES):
            i, j, k, l = QUADRUPLES[idx]
            matrix[r, i] = matrix[r, l] = 1
            matrix[r, j] = matrix[r, k] = -1
        else:
            i, j, k = TRIPLES[idx - len(QUADRUPLES)]
            matrix[r, i] = matrix[r, k] = 1
            matrix[r, j] = -2
    return sp.Matrix.hstack(*matrix.nullspace())


def plane_key(h):
    # Choose an ordered triangle of maximal (absolute) area and map its
    # vertices to (0,0), (1,0), (0,1). Minimize the resulting point list.
    translated = h - sp.ones(8, 1) * h[0, :]
    coords = translated[:, translated.rref()[1]]
    assert coords.cols == 2
    scale = sp.ilcm(*(z.q for z in coords))
    points = [tuple(int(z*scale) for z in row) for row in coords.tolist()]
    diff = [[(a[0]-b[0], a[1]-b[1]) for b in points] for a in points]

    def det(a, b):
        return a[0]*b[1] - a[1]*b[0]

    areas = {indices: abs(det(diff[j][i], diff[k][i]))
             for indices in TRIPLES for i, j, k in [indices]}
    largest = max(areas.values())
    keys = []
    for triple, area in areas.items():
        if area != largest:
            continue
        for i, j, k in permutations(triple):
            u, v = diff[j][i], diff[k][i]
            sign = 1 if det(u, v) > 0 else -1
            coords = sorted((sign*det(diff[t][i], v), sign*det(u, diff[t][i]))
                            for t in range(8))
            divisor = reduce(gcd, [area, *(z for p in coords for z in p)])
            keys.append((area // divisor, tuple((x//divisor, y//divisor) for x, y in coords)))
    return min(keys)


def read(filename='ii_min_cases.txt', threshold=12):
    groups = defaultdict(list)
    other = []
    for line in Path(filename).read_text().splitlines():
        r, n, *ids = map(int, line.split())
        if n < threshold:
            continue
        h = coordinates(ids)
        if r < 5:
            other.append((n, ids, h))
        else:
            groups[plane_key(h)].append((n, ids))
    return groups, other


def true_count(h, with_midpoints):
    from collections import Counter
    pairs = combinations(range(8), 2)
    sums = Counter(tuple(h[j, l]+h[k, l] for l in range(h.cols)) for j, k in pairs)
    result = sum(v*(v-1)//2 for v in sums.values())
    if with_midpoints:
        result += sum(sums[tuple(2*z for z in h.row(i))] for i in range(8))
    return result


if __name__ == '__main__':
    import sys
    from collections import Counter
    args=sys.argv[1:]
    groups, other = read(*([args[0]] if args else []),
                        threshold=int(args[1]) if len(args)>1 else 12)
    for n,ids,h in other:
        print('3-d shape', ids, 'exact count', true_count(h, bool(args)))
    stats = Counter()
    for group in groups.values():
        count = true_count(coordinates(group[0][1]), bool(args)); stats[count]+=1
        print('planar:', count, 'copies', len(group), 'representative', group[0][1])
    print('class distribution',sorted(stats.items()))
