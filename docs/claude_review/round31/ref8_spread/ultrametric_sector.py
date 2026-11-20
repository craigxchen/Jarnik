"""Ultrametric analogue of the sector count used in VIII.11 (referee check).
Model: O = F_9[t] = F_3[t][i], i^2 = -1 (p = 3 = 3 mod 4, infinity inert).  For monic z of degree n
(one representative per associate class), the 'argument' is u(z) = z / zbar in F_9((1/t)); its norm is 1.
A 'sector' of depth k is the set of u with prescribed coefficients of t^0..t^-k (a ball of radius 3^-(k+1)).
We count monic z of degree n in each nonempty sector and compare with the 'area' 9^n / 3^k.
"""
import itertools
from collections import Counter

# F_9 = F_3[i]/(i^2+1): elements (a, b) = a + b i
def add(x, y): return ((x[0] + y[0]) % 3, (x[1] + y[1]) % 3)
def sub(x, y): return ((x[0] - y[0]) % 3, (x[1] - y[1]) % 3)
def mul(x, y): return ((x[0] * y[0] - x[1] * y[1]) % 3, (x[0] * y[1] + x[1] * y[0]) % 3)
def conj(x): return (x[0], (-x[1]) % 3)
F9 = [(a, b) for a in range(3) for b in range(3)]
ZERO, ONE = (0, 0), (1, 0)

def series_ratio(num, den, prec):
    """num, den: coefficient lists in powers t^0, t^-1, ... with den[0] = 1. Return num/den to prec terms."""
    out = []
    rem = list(num) + [ZERO] * prec
    for k in range(prec):
        c = rem[k]
        out.append(c)
        for m, d in enumerate(den):
            if k + m < len(rem):
                rem[k + m] = sub(rem[k + m], mul(c, d))
    return tuple(out)

for n in (2, 3):
    zs = list(itertools.product(F9, repeat=n))  # z = t^n + a1 t^(n-1) + ... ; written as 1 + a1 t^-1 + ... times t^n
    for k in range(0, 2 * n + 2):
        cnt = Counter()
        for a in zs:
            num = [ONE] + list(a)
            den = [ONE] + [conj(x) for x in a]
            u = series_ratio(num, den, k + 1)
            cnt[u] += 1
        mx = max(cnt.values()); nonempty = len(cnt)
        print(f"n={n} depth k={k}: classes={len(zs)}, nonempty sectors={nonempty}, max per sector={mx}, area 9^n/3^k={9**n/3**k:.2f}")
