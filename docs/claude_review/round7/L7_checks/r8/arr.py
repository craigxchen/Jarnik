# Fixed explicit star arrangement over Q used for the write-up.
from sympy import Matrix
X3 = Matrix([[0,0,0,1]])
m1 = Matrix([[1,2,-3,0]]); m2 = Matrix([[2,-1,5,0]]); m3 = Matrix([[-3,4,1,0]])
PLANES = [X3, m1 + 2*X3, m1 - 3*X3, m2 + 5*X3, m2 - X3, m3 + 7*X3, m3 - 4*X3]  # O,a1,a2,b1,b2,c1,c2
CHOSEN = [(1,3,5),(1,3,6),(1,4,5),(1,4,6),(2,3,5),(2,4,5)]
