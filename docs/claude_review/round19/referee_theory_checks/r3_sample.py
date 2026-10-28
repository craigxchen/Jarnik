"""Referee r3_sample: random unit pencils lam*A/C, D = 1..5, with ALL candidate special lam found by
r3_lib.special_lambdas (initial-form proportionality along every D_T, and rational roots of the
discriminant on 2 random lines), plus lam = +-1 and one generic lam.
Usage: r3_sample.py chunk  (seed = chunk; 30 (A,C) pairs per chunk) -> r3_sample.jsonl (resumable).
"""
import sys, os, json, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from r3_lib import *

def rand_pair(rng):
    D = rng.choice([1, 2, 3, 4, 4, 5, 5, 5])
    npts = rng.randint(3, 7)
    pts = rng.sample(range(7), npts)
    pp = [p for p in PAIRS if p[0] in pts and p[1] in pts]
    eA = {}
    for _ in range(D):
        p = rng.choice(pp); eA[p] = eA.get(p, 0) + 1
    rest = [p for p in pp if p not in eA] or [p for p in PAIRS if p not in eA]
    eC = {}
    for _ in range(D):
        p = rng.choice(rest); eC[p] = eC.get(p, 0) + 1
    return D, eA, eC

if __name__ == "__main__":
    chunk = int(sys.argv[1])
    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(here, "r3_sample.jsonl")
    if os.path.exists(out):
        for line in open(out):
            if json.loads(line)["chunk"] == chunk:
                print("done"); sys.exit(0)
    rng = random.Random(1000 + chunk)
    recs = []
    for it in range(30):
        D, eA, eC = rand_pair(rng)
        A, C = unit_poly(eA), unit_poly(eC)
        e = dict(eA)
        for p, v in eC.items(): e[p] = e.get(p, 0) - v
        if not any(e.values()):
            continue
        k = kf_unit(e)
        cands = special_lambdas(A, C, rng) | {Fr(1), Fr(-1)}
        gen = Fr(1000003, 7919)
        dg = gain_data(A, C, gen)
        Nu = dg['n0'] + dg['ninf']
        assert dg['gain'] == -Nu and dg['n1'] == 0, ("generic lam not generic", dg)
        best = None
        for lam in cands:
            d = gain_data(A, C, lam)
            N = d['n0'] + d['ninf']
            assert N == Nu
            ok = d['gain'] <= LAMBDA[k] - N - d['n1'] and d['gain'] + d['n'] <= LAMBDA[k] and N >= 32
            row = (d['gain'], d['gain'] + d['n'] - LAMBDA[k], str(lam), ok)
            if best is None or row > best:
                best = row
            if not ok:
                print("VIOLATION", eA, eC, lam, d, flush=True)
        recs.append(dict(D=D, A=sorted(eA.items()), C=sorted(eC.items()), k=k, N=Nu, ncand=len(cands),
                         max_gain=best[0], max_excess=best[1], lam_at_max=best[2], all_ok=best[3]))
    summ = dict(chunk=chunk, pairs=len(recs), max_gain=max(r['max_gain'] for r in recs),
                max_excess=max(r['max_excess'] for r in recs), all_ok=all(r['all_ok'] for r in recs),
                special_lams=sorted({r['lam_at_max'] for r in recs if r['lam_at_max'] not in ('1', '-1')}),
                recs=[[r['D'], r['k'], r['N'], r['ncand'], r['max_gain'], r['max_excess'], r['lam_at_max']] for r in recs])
    with open(out, "a") as fh:
        fh.write(json.dumps(summ) + "\n")
    print(json.dumps({k: v for k, v in summ.items() if k != 'recs'}))
