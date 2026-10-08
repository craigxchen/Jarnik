"""Exact simple-root checks at six-point two-line conic degenerations."""

from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import gcd, lcm, prod
from random import Random

from check_segre_gradient_arithmetic import rref


def determinant(a, b, c):
    return (a[0]*(b[1]*c[2]-b[2]*c[1])
            -a[1]*(b[0]*c[2]-b[2]*c[0])
            +a[2]*(b[0]*c[1]-b[1]*c[0]))


def check(xs, ys):
    raw = []
    for values, sign in ((xs, 1), (ys, -1)):
        raw.extend(F(sign*prod(values), value*prod(value-other
                   for j, other in enumerate(values) if j != i))
                   for i, value in enumerate(values))
    denominator = lcm(*(u.denominator for u in raw))
    u = [int(v*denominator) for v in raw]
    common = reduce(gcd, u)
    u = [v//common for v in u]
    e = [[F(1)]*6, list(map(F, xs))+[F(0)]*3,
         [F(0)]*3+list(map(F, ys))]
    beta = lambda a,b: sum(c*x*y for c,x,y in zip(u,a,b))
    assert all(beta(a,b) == 0 for a in e for b in e)
    indices = (0,1,3)
    matrix = [[e[i][j]*u[j] for j in indices]+[F(i==k) for k in range(3)]
              for i in range(3)]
    inverse, pivots = rref(matrix, 3)
    assert pivots == [0,1,2]
    temp = [[F(0)]*6 for _ in range(3)]
    for i in range(3):
        for k,j in enumerate(indices):
            temp[i][j] = inverse[k][3+i]
    gram = [[beta(a,b) for b in temp] for a in temp]
    dual = [[temp[i][j]-sum(gram[i][k]*e[k][j]/2 for k in range(3))
             for j in range(6)] for i in range(3)]
    assert all(beta(e[i],dual[j]) == int(i==j) for i in range(3) for j in range(3))
    assert all(beta(a,b) == 0 for a in dual for b in dual)
    dx,dy = dual[2],[-v for v in dual[1]]
    assert beta(e[1],dy) == -1 and beta(e[2],dx) == 1

    zeros = []
    for inds in combinations(range(6),3):
        one = [F(1)]*3
        x,y = [[row[i] for i in inds] for row in e[1:]]
        value = determinant(one,x,y)
        if value == 0:
            derivative = (determinant(one,[dx[i] for i in inds],y)
                          +determinant(one,x,[dy[i] for i in inds]))
            assert derivative != 0
            zeros.append(inds)
    assert zeros == [(0,1,2),(3,4,5)]
    return all(sum(u[i] for i in inds) != 0
               for size in (1,2,3) for inds in combinations(range(6),size))


def main():
    rng = Random(20260913)
    nonresonant = 0
    for _ in range(120):
        xs = sorted(rng.sample(range(1,31),3))
        ys = sorted(rng.sample(range(1,31),3))
        nonresonant += check(xs,ys)
    assert nonresonant
    print(f'PASS: 120 exact two-line degenerations, both simple complementary '
          f'minors; {nonresonant} nonresonant weight vectors.')


if __name__ == '__main__':
    main()
