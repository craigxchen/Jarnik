"""Exact checks of signed fixed-cofactor incidence and private-prime converse."""
from itertools import combinations
from math import prod

from sympy import I, Matrix, Rational, S, cancel, conjugate, expand, im, re

marks = range(5)
edges = list(combinations(marks, 2))
blocks = dict(zip(edges, [2+I, 3+2*I, 4+I, 5+2*I, 6+I, 5+4*I,
                          7+2*I, 6+5*I, 8+3*I, 8+5*I]))
def simp(a):
    return cancel(a, extension=I)

n = {e: int(expand(b*conjugate(b))) for e, b in blocks.items()}
H = prod(n.values())
P = [expand(prod(b for e,b in blocks.items() if i in e)) for i in marks]
q = [simp(v/conjugate(v)) for v in P]
m = Matrix(5, 5, lambda i,j: im(expand(conjugate(P[i])*P[j])))
t = Matrix(5,5,lambda i,j: m[i,j]/n[tuple(sorted((i,j)))] if i!=j else 0)
assert all(t[e].is_Integer and t[e] != 0 for e in edges)
c = Matrix([1/prod(t[i,j] for j in marks if i!=j) for i in marks])
assert m*c == Matrix.zeros(5,1)
D = {e:simp((q[e[0]]-q[e[1]])/(-2*I*t[e])) for e in edges}
K = [simp(prod(d for e,d in D.items() if i in e)) for i in marks]
B = simp(H**2/prod(q))
assert simp(B**3-prod(q)/prod(K)) == 0
for i,j in edges:
    assert simp(D[i,j]*q[i]*q[j]/(K[i]*K[j]*B)-n[i,j]) == 0
    word = expand(conjugate(P[i])*P[j])/n[i,j]
    s = re(word)
    assert s.is_Integer and simp((s+I*t[i,j])/(s-I*t[i,j])-q[j]/q[i]) == 0

rows = Matrix([[S.One, z, 1/z] for z in q]).T
gale = Matrix(5,5,lambda i,j: c[i]*c[j]*m[i,j])
assert rows*gale == Matrix.zeros(3,5) or (rows*gale).applyfunc(simp) == Matrix.zeros(3,5)
print("Signed cubic inversion, unique oriented ratios and the Grassmann relations verified.")
