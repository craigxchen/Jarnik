# Exact Pr(|k(S)| > j), k uniform on {-1,0,1}^s, via trinomial coefficients (big integers); Prop 7.1 of reduction.md.
import math
from mpmath import mp, mpf, log10, sqrt, log, ceil, floor
mp.dps = 50

def trinomial_row(s):
    # coefficients of (1 + x + x^2)^s, index t+s for k(S) = t
    row = [1]
    for _ in range(s):
        new = [0] * (len(row) + 2)
        for i, c in enumerate(row):
            new[i] += c; new[i + 1] += c; new[i + 2] += c
        row = new
    return row

def params(n, r):
    d = n - 1
    w = sqrt(8 * d * log(4 * r))
    j = int(floor(w / 2)); mprime = int(ceil(4 * w))
    return d, w, j, mprime

cache = {}
for (n, r, lab) in [(2901, 2901, 'r=n'), (5000, 5000, 'r=n'), (6000, 6000 ** 2, 'r=n^2')]:
    d, w, j, mp_ = params(n, r)
    print('n=%d %s: d=%d w=%s j=%d deg=2j+1=%d m\'=%d  316 ln4r=%s (regime d>=316 ln4r? %s)' % (
        n, lab, d, mp.nstr(w, 8), j, 2 * j + 1, mp_, mp.nstr(316 * log(4 * r), 6), d >= 316 * log(4 * r)))
    for s in sorted(set([d // 2, d])):
        if s not in cache: cache[s] = trinomial_row(s)
        row = cache[s]
        tail = sum(row[t + s] for t in range(-s, s + 1) if abs(t) > j)
        p = mpf(tail) / mpf(3) ** s
        hoeff = 2 * mp.e ** (-mpf((j + 1) ** 2) / (2 * s))
        print('   s=%d: Pr(Q_S != 0) = 10^%s   Hoeffding bound 10^%s   max normalised influence <= (2j+1)/s = %s' % (
            s, mp.nstr(log10(p), 5), mp.nstr(log10(hoeff), 5), mp.nstr(mpf(2 * j + 1) / s, 4)))
