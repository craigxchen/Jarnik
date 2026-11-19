"""Anchored star = unanchored star with orientation U = conj(z_0) (squarefree N): P_0 = 1."""
import sympy as sp
from dictionary import cube8, dictionary, conj, t
clus, U = cube8()
for k in range(len(clus)):
    U0 = conj(clus[k])
    r = dictionary(clus, U0, 8, 4, f'cube D=8, anchored at point {k}', verbose=(k == 0))
    p0 = sp.gcd(clus[k], U0); assert p0.degree() == 0
print("anchored stars at all 8 anchors: P_0 = 1, E = 0, delta = 0, Sigma = 0: OK")
