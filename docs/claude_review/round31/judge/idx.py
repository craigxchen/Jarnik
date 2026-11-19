import itertools, random
from fractions import Fraction
from math import gcd
import sympy as sp
for d in range(1,9): print(d, float(sp.Rational(4,5)*(70*d)**sp.Rational(1,4)))
print('124/33=',124/33, '15/4=',3.75)
# Paley/Hadamard-type profile: rows = points, columns = singleton primes, entries 0/1
def paley(q):
    QR={(x*x)%q for x in range(1,q)}
    chi=lambda a: 0 if a%q==0 else (1 if a%q in QR else -1)
    # skew Paley tournament q=3 mod 4: rows x, columns y: entry 1 if chi(y-x)=1
    return [[1 if chi(y-x)==1 else 0 for y in range(q)] for x in range(q)]
def index(rows5):
    cols=[[r[c]-rows5[4][c] for r in rows5[:4]] for c in range(len(rows5[0]))]
    M=sp.Matrix(cols).T
    if M.rank()<4: return 0
    g=0
    for S in itertools.combinations(range(M.shape[1]),4):
        g=gcd(g,int(M[:,list(S)].det()))
        if g==1: return 1
    return g
from collections import Counter
for q in (11,19,23):
    P=paley(q); random.seed(1)
    c=Counter(index(random.sample(P,5)) for _ in range(200))
    print('paley',q,c)
