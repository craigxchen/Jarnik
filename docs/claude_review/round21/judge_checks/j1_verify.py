# Exact re-check of each J1 window with 50-digit endpoints: count split primes in [q, q*exp(1/(4r))],
# and confirm that no smaller split prime q' gives >= r (re-run in mpmath only near the boundary).
import mpmath, re
from sympy import isprime, primerange
mpmath.mp.dps = 50
for line in open('j1_out.txt'):
    r = int(re.search(r'r=(\d+)', line).group(1)); q = int(re.search(r'minP=(\d+)', line).group(1))
    hi = mpmath.mpf(q)*mpmath.e**(mpmath.mpf(1)/(4*r))
    H = int(mpmath.floor(hi))
    cnt = sum(1 for p in primerange(q, H+1) if p % 4 == 1)
    # distance from endpoint to nearest integer (sanity: no rounding issue possible if > 1e-6)
    d = min(hi-H, H+1-hi)
    # previous split prime must fail
    prev = max(p for p in primerange(5, q) if p % 4 == 1)
    hp = mpmath.mpf(prev)*mpmath.e**(mpmath.mpf(1)/(4*r))
    cp = sum(1 for p in primerange(prev, int(mpmath.floor(hp))+1) if p % 4 == 1)
    print(f"r={r} q={q} count={cnt} (>= r: {cnt>=r}) prev={prev} prev_count={cp} (< r: {cp<r}) endpoint_dist={mpmath.nstr(d,3)}")
