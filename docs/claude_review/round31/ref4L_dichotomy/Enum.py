import sympy
from run_tests import iv2
x = (-27, -270, 135, 120, 40)
from math import prod
from itertools import combinations
edges = list(combinations(range(5), 2))
for q in [30, 60]:
  for P in [10**20+39, 10**30+57, 10**40+3]:
    V, J = iv2(x, P, q)
    n = [prod(J[e] for e in edges if i in e) for i in range(5)]
    M = sympy.Matrix([[v[0]**2, v[0]*v[1], v[1]**2] for v in V[:3]])
    s = M.LUsolve(sympy.Matrix(n[:3]))
    E = [1-s[0], -s[1], 1-s[2]]
    print(q, P, [sympy.nsimplify(e) if abs(e) < 10**60 else float(e) for e in E], 'detE', E[0]*E[2]-E[1]**2/4)
