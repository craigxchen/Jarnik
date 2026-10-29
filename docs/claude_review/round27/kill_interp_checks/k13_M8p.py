"""K13: M8' = syl8 with L = 3 (|A| = 1877, rho = 11; record: zeta <= 10, ZE^all(11) holds, round 11 Thm 5.2(a)).
Criterion reach w.r.t. L_pair (mod p, standard monomials), stopping at the first failure."""
import sys, time
sys.path.insert(0, '.')
from klib import build_box, Lpair
from klib3 import shift_dims
exec(open('k6_zono.py').read().split("model = sys.argv[1]")[0])
A = build_box(8, [(S, 3) for S in syl8()])
print("|A| =", len(A), flush=True)
t0 = time.time()
for e in range(0, 12):
    hs, sm = shift_dims(A, Lpair(7), e, 2147483629, 7)
    print("M8' e=%d H=%d dim sum=%d SHIFT=%s (%.0fs)" % (e, hs, sm, hs == sm, time.time() - t0), flush=True)
    if hs != sm: break
