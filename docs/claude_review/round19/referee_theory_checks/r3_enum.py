"""Referee r3_enum: independent full enumeration of the author's c3 family, for comparison of histograms.
Family: A = product of D differences (one representative per S_7-orbit of multisets), C = product of D differences
sharing no pair with A, lam in the 14 values {+-1,+-2,+-1/2,+-3,+-1/3,+-4,+-1/4}.
Usage: r3_enum.py D orbit_index  -> appends one jsonl line to r3_enum_D{D}.jsonl (skips if done).
Records: histogram of gain, histogram of gain+n, max gain, max (gain + n - Lambda_{k(f)}).
"""
import sys, os, json, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from r3_lib import *

LAMS = [Fr(s * a, b) for s in (1, -1) for (a, b) in ((1, 1), (2, 1), (1, 2), (3, 1), (1, 3), (4, 1), (1, 4))]

def multisets(D):
    return list(itertools.combinations_with_replacement(range(21), D))

def canon(ms):
    best = None
    for perm in itertools.permutations(range(7)):
        img = tuple(sorted(PAIRS.index(tuple(sorted((perm[PAIRS[p][0]], perm[PAIRS[p][1]])))) for p in ms))
        if best is None or img < best:
            best = img
    return best

def orbit_reps(D):
    seen = {}
    for ms in multisets(D):
        c = canon(ms)
        seen.setdefault(c, ms)
    return sorted(seen)

if __name__ == "__main__":
    D = int(sys.argv[1]); oi = int(sys.argv[2])
    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(here, f"r3_enum_D{D}.jsonl")
    done = set()
    if os.path.exists(out):
        for line in open(out):
            done.add(json.loads(line)["orbit"])
    reps = orbit_reps(D)
    if oi >= len(reps):
        print("no such orbit; number of orbits", len(reps)); sys.exit(0)
    if oi in done:
        print("done already"); sys.exit(0)
    Ams = reps[oi]
    eA = {}
    for p in Ams: eA[PAIRS[p]] = eA.get(PAIRS[p], 0) + 1
    A = unit_poly(eA)
    hg, hgn = {}, {}
    maxgain = -10**9; maxexcess = -10**9; cnt = 0; worst = None
    for Cms in multisets(D):
        if set(Cms) & set(Ams):
            continue
        eC = {}
        for p in Cms: eC[PAIRS[p]] = eC.get(PAIRS[p], 0) + 1
        C = unit_poly(eC)
        e = dict(eA)
        for p, v in eC.items(): e[p] = e.get(p, 0) - v
        k = kf_unit(e)
        for lam in LAMS:
            d = gain_data(A, C, lam)
            g = d['gain']; gn = g + d['n']
            hg[g] = hg.get(g, 0) + 1; hgn[gn] = hgn.get(gn, 0) + 1
            cnt += 1
            maxgain = max(maxgain, g)
            if gn - LAMBDA[k] > maxexcess:
                maxexcess = gn - LAMBDA[k]; worst = [g, d['n'], k, [PAIRS[p] for p in Cms], str(lam)]
    rec = dict(orbit=oi, A=[PAIRS[p] for p in Ams], triples=cnt, max_gain=maxgain,
               max_gain_plus_n_minus_Lambda_k=maxexcess, worst=worst,
               hist_gain=sorted(hg.items(), reverse=True), hist_gain_plus_n=sorted(hgn.items(), reverse=True),
               norbits=len(reps))
    with open(out, "a") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(json.dumps(rec)[:400])
