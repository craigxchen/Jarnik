"""Independent referee search, one chunk: python3 rsearch.py D orbit c0 c1 OUT
For each coprime unit pair (A = rep, C = multiset of D edges outside supp A):
  generic value -N0; lambda=+-1 exact gain; all special lambda != +-1 (rational roots of the
  interpolated Disc_D on line L1), each decided exactly (sqfree on L2/L3 -> -N0, else literal formula).
Appends one json line to OUT."""
import sys, json, itertools, random, time
from fractions import Fraction
from collections import Counter
import numpy as np
W = '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round19/referee_search_checks'
sys.path.insert(0, W)
from rcore import *

D, orb, c0, c1, OUT = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
SEED = int(sys.argv[6]) if len(sys.argv) > 6 else 4242
VERIFY = float(sys.argv[7]) if len(sys.argv) > 7 else 0.0   # fraction of +-1 pencils re-done with full formula
rng = random.Random(SEED)


def rline():
    while True:
        p = [rng.randint(-1000, 1000) for _ in range(7)]
        q = [rng.randint(-1000, 1000) for _ in range(7)]
        if len(set(q)) == 7:
            return p, q


LINES = [rline() for _ in range(3)]
# random point on each L_T (|T|>=3): coordinates in T equal, all classes distinct
PT = []
for T in BIG:
    while True:
        vals = [rng.randint(-10 ** 4, 10 ** 4) for _ in range(7)]
        t0 = min(T)
        for i in T:
            vals[i] = vals[t0]
        reps_ = [vals[i] for i in range(7) if i not in T or i == t0]
        if len(set(reps_)) == len(reps_):
            break
    PT.append(vals)
VAL = np.array([[PT[k][i] - PT[k][j] for (i, j) in EDGES] for k in range(98)], dtype=object)
BIGIDX = [VIS.index(T) for T in BIG]

reps = json.load(open(W + '/rreps_D%d.json' % D))
A = [tuple(e) for e in reps[orb]]
supp = set(A)
rest = [e for e in EDGES if e not in supp]
av = vec(A)
aT = INC @ av
aT_big = aT[BIGIDX]
Apoly = poly_of(A)
Aline = [restrict(A, p, q) for (p, q) in LINES]
Aidx = [EIDX[e] for e in A]
Aval = np.prod(VAL[:, Aidx], axis=1)

hist_gen, hist_pm, hist_spec = Counter(), Counter(), Counter()
cnt = Counter()
specials, top, mism = [], [], []
t_start = time.process_time()
for k, Cm in enumerate(itertools.combinations_with_replacement(rest, D)):
    if k < c0:
        continue
    if k >= c1:
        break
    cnt['pairs'] += 1
    Cl = list(Cm)
    cv = vec(Cl)
    cT = INC @ cv
    n0 = int(np.count_nonzero(aT - cT))
    hist_gen[-n0] += 1
    if -n0 > -40:
        top.append([-n0, A, Cl, 'generic'])
    Cline = [restrict(Cl, p, q) for (p, q) in LINES]
    Cidx = [EIDX[e] for e in Cl]
    Cval = np.prod(VAL[:, Cidx], axis=1)
    cT_big = cT[BIGIDX]
    # ---- lambda = s = +-1
    for s in (1, -1):
        ok = any(sqfree_form(Cline[i] - s * Aline[i], D) for i in range(3))
        if ok:
            exc = 0
            zero = (aT_big == 0) & (cT_big == 0)
            for idx in np.nonzero(zero)[0]:
                if Cval[idx] == s * Aval[idx]:
                    cnt['S_candidates'] += 1
                    B = poly_of(Cl) - s * Apoly
                    nb = nu_T(B, BIG[idx])
                    if nb > 0:
                        exc += nb - 1
                        cnt['S_true'] += 1
            g = -n0 + exc
            if VERIFY and rng.random() < VERIFY:
                g2 = full_gain(A, Cl, s, 1)
                cnt['verified_pm'] += 1
                if g2 != g:
                    mism.append([A, Cl, s, g, g2])
        else:
            cnt['pm_fallback'] += 1
            g = full_gain(A, Cl, s, 1)
        hist_pm[g] += 1
        if g > -40:
            top.append([g, A, Cl, s])
    # ---- special lambda != +-1
    a1, c1_ = Aline[0], Cline[0]
    assert a1.degree() == D and c1_.degree() == D
    Dl = disc_poly(a1, c1_, D)
    if Dl.is_zero():
        cnt['DISC_ZERO'] += 1
        raise SystemExit('disc zero')
    for lam in rational_roots(Dl):
        if lam == 0 or lam == 1 or lam == -1:
            continue
        cnt['cand'] += 1
        p_, q_ = lam.numerator, lam.denominator
        if any(sqfree_form(q_ * Cline[i] - p_ * Aline[i], D) for i in (1, 2)):
            cnt['cand_sqfree'] += 1
            if VERIFY:
                g2 = full_gain(A, Cl, p_, q_)
                cnt['verified_cand'] += 1
                if g2 != -n0:
                    mism.append([A, Cl, [p_, q_], -n0, g2])
            continue
        g, det = full_gain(A, Cl, p_, q_, detail=True)
        cnt['special'] += 1
        hist_spec[g] += 1
        specials.append([g, A, Cl, [p_, q_], det['L'], det['F'], det['N0'], det['interior']])
rec = dict(D=D, orbit=orb, c0=c0, c1=c1, seed=SEED, A=A, counts=dict(cnt), hist_generic=dict(hist_gen),
           hist_pm=dict(hist_pm), hist_special=dict(hist_spec), specials=specials, top=top, mismatches=mism,
           cpu=time.process_time() - t_start)
with open(OUT, 'a') as f:
    f.write(json.dumps(rec, default=str) + '\n')
