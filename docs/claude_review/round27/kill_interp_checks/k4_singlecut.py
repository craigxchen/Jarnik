"""K4: single-cut boxes A = {k' in {-1,0,1}^d : |k(S)| <= L} (round-11 convention, base = last point, S a side
not containing the base).  Where does (SHIFT_e) w.r.t. L_pair fail, compared with the first slab degree 2L+1 and
with rho(A) (= top degree of H_A)?  mod p (two primes)."""
import sys, time
sys.path.insert(0, '.')
from klib import build_box, Lpair
from klib2 import shift_dims
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
cases = [(5, 4, 1), (6, 5, 1), (6, 3, 1), (7, 6, 1), (7, 4, 1), (6, 5, 2), (7, 6, 2)]
if len(sys.argv) > 1:
    cases = [tuple(int(x) for x in a.split(',')) for a in sys.argv[1:]]
for (d, s, L) in cases:
    n = d + 1
    S = tuple(range(s))
    A = build_box(n, [(S, L)])
    t0 = time.time()
    res = []; H = []
    for e in range(0, 2 * d + 2):
        hs, sm = shift_dims(A, Lpair(d), e, 2147483629, d)
        H.append(hs)
        res.append((e, hs, sm))
        if hs == 0:
            break
    fail = [e for e, hs, sm in res if hs != sm]
    pr("d=%d s=%d L=%d |A|=%d H_A=%s first SHIFT failure e=%s (2L+1=%d, rho=%d) [(e,H,sum)]=%s (%.0fs)" % (
        d, s, L, len(A), H, fail[:1], 2 * L + 1, len([h for h in H if h]) - 1, res, time.time() - t0))
open("k4_singlecut_out.txt", "a").write("\n".join(out) + "\n")
