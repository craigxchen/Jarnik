"""Referee r4: stress test of Theorem 2.4 / Corollary 2.5 (gain + n(f) <= Lambda_{k(f)}.beta <= 31) on NON-unit
pencils, where it is the only bound: A, C products of powers of translation-invariant forms (differences and
interior linear/quadratic forms), with high multiplicities, at every candidate special lam.
k(f) is bounded above by (#indices occurring) + 1 (capped at 8); Lambda of that bound is used.
Usage: r4_nonunit.py chunk -> r4_nonunit.jsonl (resumable).
"""
import sys, os, json, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from r3_lib import *

def rand_form(rng, pts):
    t = rng.random()
    if t < 0.35:
        i, j = rng.sample(pts, 2); return X[i] - X[j]
    if t < 0.7:
        i, j, k = rng.sample(pts, 3); a = rng.choice([1, 2, 3]); b = rng.choice([1, 2])
        return a * X[i] + b * X[j] - (a + b) * X[k]
    if t < 0.85 and len(pts) >= 4:
        i, j, k, l = rng.sample(pts, 4); return X[i] + X[j] - X[k] - X[l]
    i, j, k, l = rng.sample(pts, 4) if len(pts) >= 4 else (pts + pts)[:4]
    return (X[i] - X[j]) * (X[k] - X[l]) + rng.choice([1, -1, 2]) * (X[i] - X[k]) * (X[j] - X[l])

def rand_poly(rng, pts, D):
    P = CTX.from_dict({(0,) * (NV + 1): 1}); deg = 0
    while deg < D:
        f = rand_form(rng, pts)
        dd = int(f.total_degree())
        if deg + dd > D:
            continue
        m = rng.choice([1, 1, 2, 3])
        while deg + m * dd > D:
            m -= 1
        P = P * f ** m; deg += m * dd
    return P

if __name__ == "__main__":
    chunk = int(sys.argv[1])
    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(here, "r4_nonunit.jsonl")
    if os.path.exists(out):
        for line in open(out):
            if json.loads(line)["chunk"] == chunk:
                print("done"); sys.exit(0)
    rng = random.Random(5000 + chunk)
    recs = []; tested = 0
    while len(recs) < 25:
        D = rng.randint(1, 5)
        pts = rng.sample(range(7), rng.randint(3, 7))
        A, C = rand_poly(rng, pts, D), rand_poly(rng, pts, D)
        g = A.gcd(C)
        if g.total_degree() > 0:
            continue
        used = set()
        for P in (A, C):
            for mon in P.monoms():
                used |= {i for i in range(NV) if mon[i]}
        kb = min(8, len(used) + 1)
        if kb < 4:
            continue
        try:
            cands = special_lambdas(A, C, rng) | {Fr(1), Fr(-1)}
        except RuntimeError:
            continue
        best = None
        for lam in cands:
            d = gain_data(A, C, lam)
            tested += 1
            ex = d['gain'] + d['n'] - LAMBDA[kb]
            row = (ex, d['gain'], d['n'], str(lam), d['L'])
            if best is None or row > best:
                best = row
            if ex > 0:
                print("VIOLATION of Cor 2.5", str(A), str(C), lam, d, flush=True)
        recs.append([D, kb, len(cands)] + list(best))
    summ = dict(chunk=chunk, pencils=tested, pairs=len(recs), max_excess=max(r[3] for r in recs),
                max_gain=max(r[4] for r in recs), recs=recs)
    with open(out, "a") as fh:
        fh.write(json.dumps(summ) + "\n")
    print(json.dumps({k: v for k, v in summ.items() if k != 'recs'}))
