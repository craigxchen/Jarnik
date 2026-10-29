"""K10: (SHIFT_e) with the rich shift set L_cube = {-1,0,1}^d \ 0 (restricted to |A_l| > e), for given model and e.
mod p.  Usage: k10_rich.py syl8|sdc7 e"""
import sys, itertools, time
sys.path.insert(0, '.')
from klib import build_box, shiftset
from klib3 import shift_dims
exec(open('k6_zono.py').read().split("model = sys.argv[1]")[0])
name = sys.argv[1]; e = int(sys.argv[2])
A = build_box(8, [(S, 2) for S in (syl8() if name == 'syl8' else sdc7())])
L = [l for l in itertools.product((-1, 0, 1), repeat=7) if any(l) and len(shiftset(A, l)) > e]
L.sort(key=lambda l: -len(shiftset(A, l)))
t0 = time.time()
hs, sm = shift_dims(A, L, e, 2147483629, 7, verbose=True)
s = "%s e=%d with %d cube shifts: H_A(e)=%d dim sum=%d SHIFT_e=%s (%.0fs)" % (name, e, len(L), hs, sm, hs == sm, time.time() - t0)
print(s, flush=True); open("k10_%s_e%d_out.txt" % (name, e), "w").write(s + "\n")
