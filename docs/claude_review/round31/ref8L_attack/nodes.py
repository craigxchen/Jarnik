"""Node families for the two-class problem (exact rational arithmetic).

A family is: degree d, interpolation nodes s_0..s_d, extra nodes r_1..r_e (r may be 'inf'),
column weights kappa_i.  For A subset of {0..d} let f_A be the polynomial of degree <= d with
f_A(s_i) = c_i [i in A], where c_i = kappa_i * w'(s_i), w = prod (t - s_i).  At an extra node r
    f_A(r) = w(r) * sum_{i in A} kappa_i/(r - s_i)       (r finite)
    lead(f_A) = sum_{i in A} kappa_i                       (r = inf)
F = { A : for every extra node m the value lambda_m(A) lies in {0, v_m} } for a choice of v_m.
Then the |F| polynomials f_A take <= 2 values at each of the k = d+1+e nodes, and homogenising and
specialising at coprime (alpha,beta) gives |F| integers with ratio log q2 / log diam -> k/d.
"""
from fractions import Fraction as Fr
from itertools import combinations, product
import random, sys

def lam_vectors(s, r, kappa):
    X = []
    for rm in r:
        if rm == 'inf': X.append(list(kappa))
        else: X.append([kappa[i] / (rm - s[i]) for i in range(len(s))])
    return X

def best_family(s, r, kappa):
    D = len(s)
    subsets = list(range(1 << D))
    X = lam_vectors(s, r, kappa)
    sums = []
    for row in X:
        sm = {}
        for A in subsets:
            sm[A] = sum((row[i] for i in range(D) if (A >> i) & 1), Fr(0))
        sums.append(sm)
    best = [1, [0]]
    def rec(m, cand):
        if len(cand) <= best[0]: return
        if m == len(sums):
            best[0] = len(cand); best[1] = cand; return
        sm = sums[m]
        zero = [A for A in cand if sm[A] == 0]
        vals = {}
        for A in cand:
            v = sm[A]
            if v != 0: vals.setdefault(v, []).append(A)
        opts = sorted(vals.values(), key=len, reverse=True)
        if not opts: rec(m + 1, zero)
        for o in opts:
            if len(zero) + len(o) <= best[0]: break
            rec(m + 1, zero + o)
    rec(0, subsets)
    return best[0], sorted(best[1])

def search(d, e, pool, kpool, trials, seed):
    random.seed(seed)
    D = d + 1; best = (0,)
    for _ in range(trials):
        s = [Fr(0), Fr(1)] + random.sample([p for p in pool if p not in (0, 1)], D - 2)
        rest = [p for p in pool if p not in s]
        r = ['inf'] + random.sample(rest, e - 1)
        kappa = [Fr(1)] + [random.choice(kpool) for _ in range(D - 1)]
        n, F = best_family(s, r, kappa)
        if n > best[0]:
            best = (n, s, r, kappa, F)
    return best

if __name__ == "__main__":
    d = int(sys.argv[1]); e = int(sys.argv[2]); trials = int(sys.argv[3])
    pool = sorted({Fr(a, b) for a in range(-6, 7) for b in range(1, 4)})
    kpool = sorted({Fr(a, b) for a in range(-6, 7) for b in range(1, 4) if a != 0})
    b = search(d, e, pool, kpool, trials, seed=7 * d + e)
    k = d + 1 + e
    print(f"d={d} e={e} k={k} E={k}/{d}={k/d:.4f} best B={b[0]}")
    print("  s=", [str(x) for x in b[1]], "r=", [str(x) for x in b[2]], "kappa=", [str(x) for x in b[3]])
    print("  F=", b[4])
