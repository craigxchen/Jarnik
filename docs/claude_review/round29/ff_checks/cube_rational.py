# Symmetric cube ansatz: roots +-r0,+-r1,+-r4,+-r5, r_k=x_k+i y_k, x_k y_k = eps_k (eps=(+1,-1,-1,+1)),
# conditions: sum y_k = 0, sum (3/y_k - y_k^3) = 0  (k in {0,1,4,5}).  Search rational y.
from fractions import Fraction as Fr
import itertools, sympy as sp
a,b,c=sp.symbols('a b c')
d=-(a+b+c)
F=sp.together(sum(3/y-y**3 for y in (a,b,c,d)))
num=sp.factor(sp.numer(F)); print("numerator factored:",num)
