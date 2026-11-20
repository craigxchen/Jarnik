# Data on CG Conjecture 9: max #roots of x^2 = a (mod q) in a half-open real interval of length sqrt(q),
# split into unit a (gcd(a,q)=1) and all a.  Exact: interval [x, x+sqrt q) with x a root; compare via squares.
import sys
from math import isqrt, gcd
Qmax = int(sys.argv[1])
best_all = (0,None); best_unit = (0,None)
for q in range(2, Qmax+1):
    roots = {}
    for x in range(q):
        roots.setdefault(x*x % q, []).append(x)
    for a, rs in roots.items():
        if len(rs) <= best_all[0] and (gcd(a,q) != 1 or len(rs) <= best_unit[0]): continue
        ext = rs + [r + q for r in rs]
        # count roots y with x <= y < x + sqrt q  <=> (y-x)^2 < q
        m = 0; j = 0
        for i in range(len(rs)):
            if j < i: j = i
            while j + 1 < len(ext) and (ext[j+1]-ext[i])**2 < q: j += 1
            m = max(m, j - i + 1)
        if m > best_all[0]: best_all = (m, (q, a))
        if gcd(a, q) == 1 and m > best_unit[0]: best_unit = (m, (q, a))
print("q <=", Qmax, "max roots in sqrt(q)-interval: all a", best_all, "; unit a", best_unit)
