import sys
from fractions import Fraction as Fr
from beta7lift import setup, direct_nullity
from etdata import etdata
p = int(sys.argv[1])
xs = [Fr(2), Fr(3), Fr(-5), Fr(1)]
E, D = setup(xs, p)
labs = {}
for lab, A, j, P in E['pts']:
    labs[tuple((c.numerator % p) * pow(c.denominator % p, -1, p) % p for c in P)] = (lab, P)
for i, P in enumerate(E['aux']):
    labs[tuple((c.numerator % p) * pow(c.denominator % p, -1, p) % p for c in P)] = ('r%d' % (i+1), P)
for q, (lab, P) in labs.items():
    base = {r: 1 for r in D['R']}; base[q] = base.get(q, 0) + 1
    n = direct_nullity(D, base)
    if n: print(lab, [str(c) for c in P], 'nullity', n)
