"""K9: syl8 / sdc7 (L = 2) at e = 6 with ALL shifts l in {-1,0,1}^7 \ 0 (2186) and with all l in {-2..2}^7 having
|A_l| > 6: is the failure of SHIFT_6 an artefact of L_pair?  mod p."""
import sys, itertools, time
sys.path.insert(0, '.')
from klib import build_box, shiftset
from klib3 import shift_dims
exec(open('k6_zono.py').read().split("model = sys.argv[1]")[0])
name = sys.argv[1]; which = sys.argv[2]
A = build_box(8, [(S, 2) for S in (syl8() if name == 'syl8' else sdc7())])
rng = (-1, 0, 1) if which == 'cube' else (-2, -1, 0, 1, 2)
L = [l for l in itertools.product(rng, repeat=7) if any(l) and len(shiftset(A, l)) > 6]
L.sort(key=lambda l: -len(shiftset(A, l)))
t0 = time.time()
hs, sm = shift_dims(A, L, 6, 2147483629, 7, verbose=True)
print("%s e=6 with %d shifts (%s): H_A(6)=%d dim sum=%d SHIFT_6=%s (%.0fs)" % (name, len(L), which, hs, sm, hs == sm, time.time() - t0), flush=True)
open("k9_%s_%s_out.txt" % (name, which), "w").write("%s e=6 with %d shifts (%s): H_A(6)=%d dim sum=%d SHIFT_6=%s\n" % (name, len(L), which, hs, sm, hs == sm))
