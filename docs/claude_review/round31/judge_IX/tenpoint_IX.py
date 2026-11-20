# Judge: the record 10-point cluster N=1176852625 under (IX.1)/(IX.2).
import math
from math import isqrt, log
from sympy import factorint, legendre_symbol, primerange
N = 1176852625
pts = []
for s in range(48400, isqrt(2*N)+1):
    t2 = 2*N - s*s
    t = isqrt(t2)
    if t*t == t2 and (s+t) % 2 == 0:
        pts += [((s+t)//2, (s-t)//2), ((s-t)//2, (s+t)//2)]
pts = sorted(set(pts), key=lambda z: math.atan2(z[1], z[0]))
print('points', len(pts), pts[0], pts[-1])
R = math.sqrt(N)
ang = [math.atan2(y, x) for x, y in pts]
print('arc constant C =', R*(ang[-1]-ang[0])/math.sqrt(R))
xs = sorted(set(x+y for x, y in pts)); print('projections b=1+i (x+y):', xs, 'X=', xs[-1]-xs[0])
M = len(pts); n = len(xs); D = 2*N
def nu(p):
    if p == 2: return 2
    if D % p == 0: return p if p % 4 == 1 else 1
    return (p - (1 if p % 4 == 1 else -1) + 1 + legendre_symbol(D % p, p))//2
lhs = sum(log(p)*(n*n/nu(p)-n) for p in primerange(2, M+1))
rhs = sum(log(abs(a-b)) for a in xs for b in xs if a != b)
print('sieve chain', lhs, '<=', rhs, '<=', n*(n-1)*log(xs[-1]-xs[0]))
print('factor N', factorint(N))
# is the cluster symmetric under (x,y)->(y,x)?
print('diag-symmetric:', set(pts) == set((y, x) for x, y in pts))
