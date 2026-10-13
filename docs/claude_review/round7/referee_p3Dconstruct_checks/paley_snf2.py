# Greedy: largest row subset R of Paley 0/1 matrix P with column lattice = Z^R (index 1), i.e.
# rows R surjective over Z (hence mod every m_q/2).  Referee check of Consequence 4 caveat.
import random, sys
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round7/referee_p3Dconstruct_checks')
from paley_snf import paley, det_divisor_gcd, rank_mod
random.seed(7)
for q0 in [43, 47, 59, 67, 71, 79, 83, 103, 107]:
    H = paley(q0); M = q0 + 1
    P = [[(1 + H[x][a])//2 for a in range(1, M)] for x in range(M)]
    best = []
    for trial in range(6):
        order = list(range(M)); random.shuffle(order)
        R = []
        for x in order:
            if det_divisor_gcd([P[y] for y in R + [x]]) == 1:
                R.append(x)
        if len(R) > len(best): best = R[:]
    print(f"q0={q0} (mod 8 = {q0%8}) M={M}: largest Z-surjective row subset found |R|={len(best)} ({len(best)/M:.3f} M)")
