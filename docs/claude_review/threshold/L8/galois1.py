from c5ansatz import *
import random
from math import isqrt
random.seed(3)
def interp(xs, ys):
    n = len(xs); res = []
    for i in range(n):
        num = [Fr(1)]; den = Fr(1)
        for j in range(n):
            if j != i:
                num = mul(num, linpoly(1, -xs[j])); den *= xs[i] - xs[j]
        res = add(res, scal(ys[i] / den, num))
    return res
def issq(fr):
    if fr < 0: return False
    a, b = fr.numerator, fr.denominator
    return isqrt(a)**2 == a and isqrt(b)**2 == b
cnt = 0
for trial in range(200):
    vals = random.sample(range(-30, 31), 9)
    rho = {e: Fr(v) for e, v in zip(EDGES[1:], vals[:4])}
    sigma = {T: Fr(v) for T, v in zip(TRIPLES, vals[4:])}
    def N(x):
        r = dict(rho); r[(1, 2)] = x
        D = prod([linpoly(-1, sigma[T]) for T in triples_of(1)] + [linpoly(-1, sigma[T]) for T in triples_of(2)])
        return det_frac(matrix(r, sigma)) * ev(D, x)
    xs = [Fr(10**6 + 17 * k, 1) for k in range(5)]
    p = interp(xs, [N(x) for x in xs])
    q, rem = divmod_(p, linpoly(1, -sigma[(1, 2, 4)]))
    if deg(q) != 3: continue
    d, c, b, a = q
    disc = 18*a*b*c*d - 4*b**3*d + b**2*c**2 - 4*a*c**3 - 27*a**2*d**2
    cnt += issq(disc)
print('square discriminants:', cnt, 'of 200')
