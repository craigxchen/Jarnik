"""Integral certificates for five-mark restrictions of the skew double core."""
from itertools import combinations
from math import comb

from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form


def double_core(q):
    # The border is row 0, and field element t has row t+1.
    squares = {x*x % q for x in range(1, q)}
    def chi(x):
        if x % q == 0:
            return 0
        return 1 if x % q in squares else -1
    matrix = [[0] * (q+1) for _ in range(q+1)]
    matrix[0] = [1] + [-1]*q
    for x in range(q):
        matrix[x+1][0] = 1
        for y in range(q):
            matrix[x+1][y+1] = chi(y-x) + int(y == x)
    first = [[(matrix[x][a]+1)//2 for x in range(q+1)] for a in range(1,q+1)]
    second = [[(matrix[a][x]*matrix[0][x]+1)//2 for x in range(q+1)] for a in range(1,q+1)]
    return first, second


def verify(q):
    first, second = double_core(q)
    m = q+1
    columns = first + second
    differences = [[a-b for a,b in zip(c,d)] for c,d in zip(first,second)]
    for a, d in enumerate(differences, start=1):
        assert d == [int(i == a)-int(i == 0) for i in range(m)]
    assert all(sum(c) == m//2 for c in columns)
    anchored = Matrix([[c[i]-c[0] for c in columns] for i in range(1,m)])
    full_index = hermite_normal_form(anchored).det()
    assert full_index == m//2
    count = 0
    for rows in combinations(range(m),5):
        # Left multiplication by the anchoring map sends Z^5 onto Z^4.
        anchor = rows[0]
        sub = Matrix([[c[i]-c[anchor] for c in columns] for i in rows[1:]])
        assert hermite_normal_form(sub) == Matrix.eye(4)
        count += 1
    assert count == comb(m,5)
    print(f"m={m}: whole index {full_index}; all {count} five-subsets have index 1")


if __name__ == '__main__':
    for field in (7,11):
        verify(field)
