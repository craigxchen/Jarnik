"""Evaluate the determinant of quadrics on the ten projected pair sums."""
from itertools import combinations, combinations_with_replacement
from sympy import Matrix, eye, ones
v=list(eye(4).columnspace())+[-ones(4,1)]
p=[v[i]+v[j] for i,j in combinations(range(5),2)]
M=Matrix([[a[j]*a[k] for j,k in combinations_with_replacement(range(4),2)] for a in p])
print(M.det())
assert M.det()==-3
