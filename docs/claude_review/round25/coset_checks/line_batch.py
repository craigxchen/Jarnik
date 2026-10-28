"""Section 8.2 batch: 60 random lines (b random primitive in [-3,3]^8, g uniform) and 30 lines through the
transplant point xb (random b).  Primes 13..73, T = 22.  Reports character status and surviving length."""
import sys, math, json
import numpy as np
sys.path.insert(0, '.')
from common import *
from line_cosets import line_survivors
allsp = split_primes(200)
ps = [13, 17, 29, 37, 41, 53, 61, 73]
qs = [allsp[allsp.index(p) - 1] for p in ps]
xb = np.array([math.atan2(*gaussian_prime(q)[::-1]) / (math.pi / 2) for q in qs])
V, W = monomials([math.log(p) for p in ps], 22.0)
rng = np.random.default_rng(11)
def rand_b():
    while True:
        b = rng.integers(-3, 4, size=8)
        if np.gcd.reduce(np.abs(b[b != 0])) == 1 and (b != 0).sum() >= 2:
            return b
stats = {'random_g': [], 'through_transplant': []}
for i in range(60):
    b = rand_b(); g = rng.uniform(0, 1, 8)
    r = line_survivors(V, W, b, g)
    stats['random_g'].append((r['char_margin'], r['surv_len']))
for i in range(30):
    b = rand_b()
    r = line_survivors(V, W, b, xb)
    stats['through_transplant'].append((r['char_margin'], r['surv_len']))
rg = stats['random_g']
summary = {
 'random_g_char_ok': sum(1 for c, s in rg if c is None or c > 0),
 'random_g_char_ok_and_survivors': sum(1 for c, s in rg if (c is None or c > 0) and s > 0),
 'random_g_char_fail_but_noncharacter_survivors': sum(1 for c, s in rg if (c is not None and c <= 0) and s > 0),
 'transplant_lines_survivors': sum(1 for c, s in stats['through_transplant'] if s > 0),
 'transplant_lines_min_len': min(s for c, s in stats['through_transplant']),
 'transplant_lines_char_ok': sum(1 for c, s in stats['through_transplant'] if c is None or c > 0)}
print(json.dumps(summary, indent=1))
json.dump({'summary': summary, 'stats': stats}, open('line_batch_out.json', 'w'), indent=1)
