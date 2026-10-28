# Exact: Vandermonde prod_{i<j}(x_i-x_j) (5 vars), shifted by -(2,2,2,2,2): support in 3P cap M, ord at e = 10 exactly.
import itertools, sympy as sp
from fractions import Fraction
x = sp.symbols('x0:5')
V = sp.Poly(sp.prod([x[i]-x[j] for i in range(5) for j in range(i+1,5)]), *x)
terms = [(tuple(e[t]-2 for t in range(5)), int(c)) for e, c in V.terms()]
A_list = [A for r in range(1,5) for A in itertools.combinations(range(5), r)]
inside = all(sum(u[i] for i in A) >= -3 for u,_ in terms for A in A_list)
print("support size", len(terms), "inside 3P:", inside, " inside 2P:", all(sum(u[i] for i in A) >= -2 for u,_ in terms for A in A_list))
def moment(alpha):
    return sum(c*sp.prod([u[i]**alpha[i] for i in range(5)]) for u,c in terms)
from itertools import product
zero_below = all(moment(a)==0 for d in range(10) for a in product(range(d+1), repeat=5) if sum(a)==d)
nonzero10 = any(moment(a)!=0 for a in product(range(11), repeat=5) if sum(a)==10)
print("all moments of degree <10 vanish:", zero_below, "; some degree-10 moment nonzero:", nonzero10)
