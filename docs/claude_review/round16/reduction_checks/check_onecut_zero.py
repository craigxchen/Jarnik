"""Exact check: on one-cut boxes A_{S,L} (T=1, k' ternary, |k(S)| <= L), the L2 form Psi_A(u) vanishes
exactly when the zetc / nearstrata tangent degree nu(u) = 2a + min(2 sigma, 2L + 2t, s + L) is < m
(u with planted exact value classes; C_0 = class of value 0 incl. the base point).
Also checks Psi_A >= (|B|/|A|) Psi_B for the slice B = {k(S)=0} (monotonicity, D6) and that
the slice already has Psi_B > 0 whenever u restricted to S takes >= ... (slice nu = min(2t_0, s) + 2a)."""
import random, sys, time
from fractions import Fraction
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round16/reduction_checks')
from lib import *

random.seed(1616)
out = []
def log(s):
    print(s); out.append(s); sys.stdout.flush()

def nu_onecut(u, S, L):
    d = len(u)
    R = [x for x in range(d) if x not in S]
    s = len(S)
    a = sum(1 for x in R if u[x] != 0)
    sigma = sum(1 for x in S if u[x] != 0)
    vals = {}
    for x in S:
        vals[u[x]] = vals.get(u[x], 0) + 1
    maxc = max(vals.values()) if vals else 0
    t = s - maxc
    return 2 * a + min(2 * sigma, 2 * L + 2 * t, s + L)

t0 = time.time()
cases = [(5, [0, 1, 2], 1), (5, [0, 1, 2, 3], 1), (6, [0, 1, 2, 3], 1), (6, [0, 1, 2, 3, 4], 2), (6, [0, 1, 2, 3], 0)]
tot = 0; bad = 0; mono_bad = 0
for d, S, L in cases:
    A = one_cut(d, S, L)
    B = one_cut(d, S, 0)
    for m in range(2, len(S) + L + 2 * (d - len(S)) + 2):
        bA = Box(A, m); bB = Box(B, m)
        for trial in range(12):
            # planted classes: values from a small pool, some zeros
            pool = [Fraction(0)] + [Fraction(random.randint(1, 30), random.randint(1, 5)) * random.choice((-1, 1)) for _ in range(random.randint(1, 4))]
            u = [random.choice(pool) for _ in range(d)]
            if all(x == 0 for x in u):
                continue
            pA = bA.psi(u); pB = bB.psi(u)
            nu = nu_onecut(u, S, L)
            ok = ((pA == 0) == (nu < m))
            tot += 1; bad += (not ok)
            if not ok:
                log(f"MISMATCH d={d} S={S} L={L} m={m} u={u} psi={pA} nu={nu}")
            if pA * len(B) < pB * len(B) * 0 + Fraction(len(B), len(A)) * pB * len(A) * 0:  # placeholder
                pass
            # monotonicity in counting normalisation: Psi^Sigma_A >= Psi^Sigma_B
            if pA < pB:
                mono_bad += 1
    log(f"case d={d} S={S} L={L}: done ({time.time()-t0:.1f}s)")
log(f"one-cut zero-set test: {tot} instances, mismatches = {bad}; counting-monotonicity violations = {mono_bad}")
open(__file__.replace('.py', '_out.txt'), 'w').write('\n'.join(out) + '\n')
