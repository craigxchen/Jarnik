"""(B) near-miss: n = prod_{i=1}^t (1 + i q), q divisible by t!, has >= 2^t divisors = 1 mod q
(the subset products; factors pairwise coprime since gcd(1+iq,1+jq) | (j-i) and q is divisible by
every prime <= t).  alpha = log q / log n -> 1/t from below.  Exact check, t = 4, 5."""
from math import gcd, log, prod
from itertools import combinations
for t in (4, 5):
    for e in (10, 40, 200):
        q = 120 * (10 ** e + 1)          # divisible by 5! = 120
        f = [1 + i * q for i in range(1, t + 1)]
        assert all(gcd(a, b) == 1 for a, b in combinations(f, 2))
        n = prod(f)
        ds = {prod(c) for k in range(t + 1) for c in combinations(f, k)}
        assert all(d % q == 1 for d in ds) and all(n % d == 0 for d in ds)
        print(f"t={t} q~10^{e+2}: {len(ds)} divisors = 1 mod q, alpha = log q/log n = {log(q)/log(n):.6f} (1/t = {1/t:.6f})")
