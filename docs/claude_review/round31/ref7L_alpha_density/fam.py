# Own exact construction of the five-mark polynomial Pell family.
import sympy as sp
from sympy import I, Rational as Q, symbols, expand, Poly, Matrix
t = symbols('t', real=True)
S = symbols('s1:5', real=True)
EDGES = [(0,1),(0,2),(1,2),(0,3)]

def conjp(e):
    # coefficientwise complex conjugation; t and s real
    return sp.conjugate(e)

def marks(n, seeds, tt=t):
    A = Matrix([[tt**2+1, tt],[tt,1]])
    An = A**n
    G = []
    for s in seeds:
        v = An*Matrix([1, s])
        G.append(expand(v[0] + I*v[1]))
    Gb = [expand(conjp(g)) for g in G]
    beta = tt + 2*I; betab = tt - 2*I
    z0 = beta
    for g in Gb: z0 *= g
    zs = [expand(z0)]
    for (j,k) in EDGES:
        z = betab*G[j]*G[k]
        for l in range(4):
            if l not in (j,k): z *= Gb[l]
        zs.append(expand(z))
    return G, Gb, zs
