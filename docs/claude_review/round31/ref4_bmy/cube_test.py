import sympy as sp
from fractions import Fraction as Fr
from genuine import cluster
from bmy_global import analyse
Y = [sp.Rational(-9,4), sp.Rational(-2,3), sp.Rational(1,4), sp.Rational(8,3)]
roots = [e/y + sp.I*y for y in Y for e in (1,-1)]
cols, N, nl = cluster(roots, 4)
print('cluster size', len(cols))
res = analyse(cols, N, nl)
print({k: (str(v) if isinstance(v, Fr) else v) for k, v in res.items()})
