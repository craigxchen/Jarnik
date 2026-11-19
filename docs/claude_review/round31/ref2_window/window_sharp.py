"""Sharpness of the exponent 1/6 in the six-norm lemma.
If F(Y)=prod(Y+k_v) (min k = 0), Q its cubic square-root part, R=F-Q^2=r2 Y^2+r1 Y+r0, r0=-Q(0)^2.
If r2=0 and X=(Q(0)^2+t^2)/r1 is a positive integer with Q(X) integral, then F(X)=Q(X)^2+t^2.
Search such k-vectors and report the window ratio K/X^(1/6) and t/sqrt(X)."""
from fractions import Fraction as Fr
from itertools import combinations_with_replacement
import math

def esym(ks):
    e = [1] + [0]*6
    for k in ks:
        for j in range(6, 0, -1): e[j] += e[j-1]*k
    return e

best = []
for K in range(1, 61):
    for rest in combinations_with_replacement(range(0, K+1), 4):
        ks = (0,) + rest + (K,)
        e = esym(ks)
        a = Fr(e[1], 2); b = Fr(e[2], 2) - Fr(e[1]**2, 8)
        c = Fr(e[3], 2) - Fr(e[1]*e[2], 4) + Fr(e[1]**3, 16)
        r2 = e[4] - (b*b + 2*a*c); r1 = e[5] - 2*b*c; r0 = e[6] - c*c
        assert r0 == -c*c
        if r2 != 0 or r1 <= 0: continue
        for t in range(1, 200):
            Xf = (c*c + t*t) / r1
            if Xf.denominator != 1 or Xf <= K: continue
            X = int(Xf)
            QX = X**3 + a*X**2 + b*X + c
            if QX.denominator != 1: continue
            Fv = math.prod(X+k for k in ks)
            assert Fv == int(QX)**2 + t*t
            best.append((K / X**(1/6), t / math.sqrt(X), X, ks, t, r1))
            break
best.sort()
print('number of exact examples found:', len(best))
for row in best[:12]:
    print('K/X^(1/6)=%.3f  t/sqrtX=%.2e  X=%d  k=%s  t=%d  r1=%s' % row)
