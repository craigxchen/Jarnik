# Round-19 abc-gain search, one chunk: D, orbit index o, C-range [c0, c1).  Appends one JSON line to OUT.
# For every pair (A, C) (A = orbit rep, C = multiset of D edges disjoint from A, in itertools order):
#   generic gain  = -N0(A,C)                              (all lam outside the special set; proof in search.md §2)
#   lam = +1, -1 : exact gain; fast path (B squarefree via a line, S by exact initial-form comparison,
#                  nu_T(B) exact by substitution), fallback = exact.exact_gain (referee formula)
#   lam != 0, +-1: candidates = rational roots of Delta_1(lam) = Disc_t(c(t) - lam a(t)) on line 1 (complete
#                  superset, search.md Thm 2.2); each candidate: B squarefree on line 2 or 3 => gain = -N0;
#                  else exact.exact_gain.
import sys, json, itertools, collections, time, random
import numpy as np
from flint import fmpz_poly, fmpz_mpoly_ctx
from fractions import Fraction
from common import edges, Ts, orbit_reps, red_init, N0
import exact

D, o, c0, c1, OUT = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
import os
VERIFY = (len(sys.argv) > 6 and sys.argv[6] == 'verify') or os.environ.get('VERIFY') == '1'   # also run exact_gain on every lam=+-1 pencil
reps = orbit_reps(D, '.')
A = reps[o]; Aset = set(A)
others = [e for e in edges if e not in Aset]
Cs = list(itertools.combinations_with_replacement(others, D))[c0:c1]

# visible T list: 98 Ts then 21 pairs
VT = Ts + [frozenset(e) for e in edges]
M = np.array([[1 if (e[0] in T and e[1] in T) else 0 for T in VT] for e in edges], dtype=np.int64)
eidx = {e: k for k, e in enumerate(edges)}
def cvec(ms): return M[[eidx[e] for e in ms]].sum(axis=0)
# hashes of reduced initial factors along each T in Ts (prefilter only; equality verified exactly)
rng = random.Random(12345); W = {}
def wfac(f):
    if f not in W: W[f] = rng.getrandbits(62)
    return W[f]
H = np.zeros((21, 98), dtype=np.uint64)
for k, e in enumerate(edges):
    for t, T in enumerate(Ts):
        H[k, t] = wfac(red_init([e], T)[1][0])
def hvec(ms): return H[[eidx[e] for e in ms]].sum(axis=0, dtype=np.uint64)

# lines x = p + t q (fixed seed)
import os
lr = random.Random(int(os.environ.get('LINESEED', '2026')))
LINES = [([lr.randint(-60, 60) for _ in range(7)], [lr.randint(-60, 60) for _ in range(7)]) for _ in range(3)]
LP = [{e: fmpz_poly([p[e[0]] - p[e[1]], q[e[0]] - q[e[1]]]) for e in edges} for (p, q) in LINES]
def restr(ms, k):
    r = fmpz_poly([1])
    for e in ms: r *= LP[k][e]
    return r
def sqfree_binary(b, D):
    """True iff the binary form of formal degree D with dehomogenisation b is squarefree on P^1."""
    if b.degree() < D - 1: return False          # double root at infinity (or b == 0)
    return b.gcd(b.derivative()).degree() == 0
ctx2 = fmpz_mpoly_ctx.get(('t', 'L'), 'lex')

aA = cvec(A); hA = hvec(A)
rA = [restr(A, k) for k in range(3)]
redA = {T: red_init(A, T) for T in Ts}

hist_pm = collections.Counter(); hist_gen = collections.Counter(); hist_sp = collections.Counter()
recs_top = []; recs_sp = []; stats = collections.Counter(); t0 = time.time()
for ci, C in enumerate(Cs):
    aC = cvec(C); n0 = int((aA != aC).sum())
    hist_gen[-n0] += 1
    hC = hvec(C)
    cand = [Ts[t] for t in np.nonzero(hA == hC)[0]]
    Sset = {1: [], -1: []}
    for T in cand:
        sa, fa = redA[T]; sc, fc = red_init(C, T)
        if fa == fc: Sset[sc * sa].append(T)       # in_T(C) = (sc*sa) in_T(A)
        else: stats['hash_false_pos'] += 1
    rC = [restr(C, k) for k in range(3)]
    for s in (1, -1):
        g = None
        for k in range(3):
            if sqfree_binary(rC[k] - s * rA[k], D):
                Bm = exact.poly(C) - s * exact.poly(A) if Sset[s] else None
                g = -n0 + sum(exact.nuT(Bm, T) - 1 for T in Sset[s])
                stats['fast'] += 1; break
        if g is None:
            g = exact.exact_gain(A, C, s, 1); stats['fallback'] += 1
        if VERIFY:
            ge = exact.exact_gain(A, C, s, 1)
            if ge != g: stats['VERIFY_MISMATCH'] += 1; recs_top.append(['MISMATCH', A, C, s, g, ge])
            else: stats['verified'] += 1
        g = int(g); hist_pm[g] += 1
        if g > -40: recs_top.append([g, A, C, s])
    # lam != 0, +-1
    F = ctx2.from_dict({**{(i, 0): int(cf) for i, cf in enumerate(rC[0].coeffs()) if cf != 0}})
    F = F - ctx2.from_dict({(i, 1): int(cf) for i, cf in enumerate(rA[0].coeffs()) if cf != 0})
    if F.degrees()[0] != D: stats['deg_drop_line1'] += 1   # formal degree D must be the t-degree
    Dl = F.discriminant('t')
    if Dl.is_zero():
        stats['DISC_ZERO_line1'] += 1; recs_top.append(['DISC_ZERO', A, C]); continue
    dp = fmpz_poly([0] * (Dl.degrees()[1] + 1))
    coeffs = [0] * (Dl.degrees()[1] + 1)
    for mon, cf in zip(Dl.monoms(), Dl.coeffs()): coeffs[mon[1]] += int(cf)
    dp = fmpz_poly(coeffs)
    lams = set()
    for fac, m in dp.factor()[1]:
        if fac.degree() == 1:
            v, u = int(fac.coeffs()[0]), int(fac.coeffs()[1])
            lam = Fraction(-v, u)
            if lam not in (0, 1, -1): lams.add(lam)
    for lam in lams:
        stats['cand_lam'] += 1; p_, q_ = lam.numerator, lam.denominator
        sq = any(sqfree_binary(q_ * rC[k] - p_ * rA[k], D) for k in (1, 2))
        if sq:
            g = -n0; stats['cand_sqfree'] += 1          # B_lam squarefree: gain = -N0 (search.md Prop 2.1)
            if VERIFY:
                ge = exact.exact_gain(A, C, p_, q_)
                if ge != g: stats['VERIFY_MISMATCH_cand'] += 1; recs_top.append(['MISMATCH', A, C, [p_, q_], g, ge])
                else: stats['verified_cand'] += 1
        else:
            g = exact.exact_gain(A, C, p_, q_); stats['cand_exact'] += 1
            g = int(g); recs_sp.append([g, A, C, [p_, q_]]); hist_sp[g] += 1
        if g > -40: recs_top.append([g, A, C, [p_, q_]])
out = dict(D=D, orbit=o, c0=c0, c1=c1, n=len(Cs), A=A, hist_pm=dict(hist_pm), hist_generic=dict(hist_gen),
           hist_special=dict(hist_sp), top=recs_top[:2000], special=recs_sp, stats=dict(stats),
           time=round(time.time() - t0, 1))
with open(OUT, 'a') as fh:
    fh.write(json.dumps(out, default=lambda x: str(x)) + '\n')
