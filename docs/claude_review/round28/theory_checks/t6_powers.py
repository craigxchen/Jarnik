"""T6: which obstructions are powers of linear forms?  For a box A and degree e (mod p): Sigma = sum_{l in L_pair}
M_e(A_l), MA = M_e(A).  phi^e (phi in {-1,0,1}^d, or {-2..2}^d with --wide) is an OBSTRUCTION iff phi^e perp Sigma and
phi^e not perp MA.  Reports obstruction phi's, their ranges on A and A_l, and dim span of their classes vs delta(e).
Usage: t6_powers.py spec e [wide]   spec = syl8:L | sdc7:L | rand:n,seed (k11 convention)"""
import sys, time, random
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/theory_checks')
from tlib import *
from math import factorial
spec, e = sys.argv[1], int(sys.argv[2]); wide = len(sys.argv) > 3
p = P1
def k11_box(seed, n):
    rng = random.Random(seed); d = n - 1
    r = rng.randint(n - 1, 2 * n); pats = set()
    while len(pats) < r:
        S = tuple(x for x in range(n) if rng.random() < 0.5)
        if 0 < len(S) < n: pats.add(canon(n, S))
    L = rng.choice([1, 2])
    return n, [(S, L) for S in sorted(pats)]
kind, arg = spec.split(':')
if kind in ('syl8', 'sdc7'):
    n, pats = syl(3) if kind == 'syl8' else sdc(7); slabs = [(S, int(arg)) for S in pats]
else:
    nn, sd = (int(x) for x in arg.split(',')); n, slabs = k11_box(sd, nn)
d = n - 1
A = build_box(n, slabs); Aset = set(A)
Lp = Lpair_half(d)
MA = Me_modp(A, e, p, d)
HA = 0 if MA is None else MA.rank()
Sg = shift_sum_modp(A, Lp, e, p, d)
dS = 0 if Sg is None else Sg.ncols()
print(spec, "e=%d |A|=%d H_A(e)=%d dim Sigma=%d delta=%d" % (e, len(A), HA, dS, HA - dS), flush=True)
mons = smonos(d, e)
Ex = np.array(mons, dtype=np.int64)
mult = np.array([factorial(e) // int(np.prod([factorial(int(a)) for a in al])) for al in mons], dtype=object)
vals = (-2, -1, 0, 1, 2) if wide else (-1, 0, 1)
phis = [ph for ph in itertools.product(vals, repeat=d) if any(ph) and next(v for v in ph if v) > 0]
rows = []
for ph in phis:
    row = []
    for i, al in enumerate(mons):
        v = int(mult[i])
        for x in range(d):
            if al[x]: v *= int(ph[x]) ** int(al[x])
        row.append(v % p)
    rows.append(row)
PH = flint.nmod_mat(len(rows), len(mons), [v for r in rows for v in r], p)
PS = None if Sg is None else PH * Sg
PM = PH * MA
obs = []
for i, ph in enumerate(phis):
    inSig = True if PS is None else all(int(PS[i, j]) == 0 for j in range(PS.ncols()))
    if not inSig: continue
    if any(int(PM[i, j]) != 0 for j in range(PM.ncols())):
        obs.append((ph, i))
print("number of phi with phi^e an obstruction:", len(obs), flush=True)
if obs:
    # dimension of span of obstruction classes modulo J(A)_e: rank of their pairings with MA
    P_ = flint.nmod_mat(len(obs), MA.ncols(), [int(PM[i, j]) for ph, i in obs for j in range(MA.ncols())], p)
    print("dim span of power obstructions mod J(A)_e:", P_.rank(), " (delta =", HA - dS, ")")
    for ph, c in obs[:12]:
        rA = sorted(set(sum(a * b for a, b in zip(ph, k)) for k in A))
        rl = []
        for l in Lp:
            Al = shiftset(A, l, Aset)
            vs = set(sum(a * b for a, b in zip(ph, k)) for k in Al)
            rl.append(len(vs))
        print("  phi=%s range on A: %d values [%d,%d]; max #values on A_l (l in L_pair): %d" % (ph, len(rA), rA[0], rA[-1], max(rl)))
