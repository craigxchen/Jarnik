"""Exhaustive small-height census (Data).

All configurations of k distinct points of P^1(Q) that are represented by integers
0 = x_0 < x_1 < ... < x_{k-1} <= B with gcd 1 (infinity unmarked).  For each, the anchor x_0 is
moved to infinity (X_i = L/(x_i - x_0), L = lcm, divided by the gcd of differences), and the
119/56/25/10 visible contact integers d_T are computed exactly (census_lib).  Recorded:
  cov  = #{T : d_T > 1}            (coverage of the boundary divisors),
  eps  = (ell_max - ell_min)/(ell_max + ell_min)  (=1 if cov < V),
  w*   = (ell_max + ell_min)/2,
  adm  = admissible (two chains at every prime, anchor x_0).
Output: distribution of cov, and the frontier min{eps : w* >= W} for full-coverage points,
plus the best configurations.  Usage: python3 census_small.py k B
"""
import sys, time, json
from itertools import combinations
from math import gcd, log
from functools import reduce
from census_lib import contact_integers, visible_blocks, admissible_at_prime

def lcm(a, b):
    return a * b // gcd(a, b)

def anchored(xs):
    x0 = xs[0]
    diffs = [x - x0 for x in xs[1:]]
    L = reduce(lcm, diffs)
    X = [L // dd for dd in diffs]
    g = 0
    for i in range(len(X)):
        for j in range(i + 1, len(X)):
            g = gcd(g, X[i] - X[j])
    return [x // g for x in X]

def main(k, B):
    m = k - 1
    V = 2 ** m - m - 2
    t0 = time.time()
    covhist = {}
    full = []
    n = 0
    for rest in combinations(range(1, B + 1), k - 1):
        if rest[-1] != B and B > k:  # enumerate by max = b for incremental output; handled below
            pass
        xs = (0,) + rest
        if reduce(gcd, xs) != 1:
            continue
        n += 1
        X = anchored(xs)
        d, top, l = contact_integers(X)
        cov = sum(1 for v in d.values() if v > 1)
        covhist[cov] = covhist.get(cov, 0) + 1
        if cov == V:
            lg = [log(v) for v in d.values()]
            lo, hi = min(lg), max(lg)
            adm = all(admissible_at_prime(list(lp.keys()), m) for lp in l.values())
            full.append(((hi - lo) / (hi + lo), (hi + lo) / 2, xs, adm, lo, hi))
    full.sort()
    res = {"k": k, "B": B, "configs": n, "cov_hist": dict(sorted(covhist.items())),
           "full_coverage": len(full), "best": [(round(e, 4), round(w, 3), xs, adm, round(lo, 3), round(hi, 3))
                                                 for e, w, xs, adm, lo, hi in full[:15]],
           "seconds": round(time.time() - t0, 1)}
    # frontier: min eps among full-coverage points with w* >= W
    fr = {}
    for W in [1, 2, 3, 4, 5, 6, 8, 10]:
        c = [e for e, w, *_ in full if w >= W]
        fr[W] = round(min(c), 4) if c else None
    res["frontier_min_eps_for_wstar_ge"] = fr
    print(json.dumps(res))

if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]))
