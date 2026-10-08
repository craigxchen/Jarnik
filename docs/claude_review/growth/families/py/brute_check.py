# Brute-force check of the Rust scanner on ALL n<=N (no restriction to primes 1 mod 4).
import math, sys
from collections import defaultdict
N = int(sys.argv[1])
pts = defaultdict(list)
r = int(math.isqrt(N))
for x in range(-r, r+1):
    for y in range(-r, r+1):
        n = x*x+y*y
        if 0 < n <= N:
            pts[n].append(math.atan2(y, x))
KM = 8
first = {}
blockmin = {}
for n in sorted(pts):
    a = sorted(pts[n]); d = len(a)
    n14 = n ** 0.25
    for k in range(3, min(d, KM)+1):
        best = min(((a[(i+k-1) % d] - a[i]) % (2*math.pi)) for i in range(d))
        c = best * n14
        j = n.bit_length()-1
        if c < blockmin.get((k, j), (9e9,))[0]:
            blockmin[(k, j)] = (c, n)
        for ct in (0.5, 1.0, 1.4142):
            if c <= ct and (k, ct) not in first:
                first[(k, ct)] = n
for key in sorted(first): print("first k=%d C<=%s n=%d" % (key[0], key[1], first[key]))
for key in sorted(blockmin): print("blockmin k=%d j=%d min=%.6f n=%d" % (key[0], key[1], blockmin[key][0], blockmin[key][1]))
