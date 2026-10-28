# Referee check r2b: structured (near-extremal) supports for (5.1) and reg < 2 K1.
import sys; sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round25/referee_intra_checks')
from rlib import *
from itertools import product
from r2_regB import cuts
def zs_box(M, T):
    return [k + (-sum(k),) for k in product(range(-T, T+1), repeat=M-1) if abs(sum(k)) <= T]
def l1ball(M, d):
    return [k + (-sum(k),) for k in product(range(-2*d, 2*d+1), repeat=M-1) if sum(map(abs, k))+abs(sum(k)) <= 2*d]
for name, M, S in [("1*R_5", 5, l1ball(5, 1)), ("2*R_5", 5, l1ball(5, 2)), ("cube{-1,0,1}^5 slice", 5, zs_box(5, 1)),
                   ("1*R_6", 6, l1ball(6, 1)), ("cube{-1,0,1}^6 slice", 6, zs_box(6, 1)), ("1*R_7", 7, l1ball(7, 1))]:
    r = reg(S); cs = cuts(M); B = 2**(2-M)*sum(width(S, A) for A in cs); K1 = max(sum(map(abs, k)) for k in S)
    print(f"{name:24s} |S|={len(S):4d} reg={r} bound(5.1)={B:.4f} ok={r <= B+1e-12}  K1={K1} reg<2K1={r < 2*K1}", flush=True)
