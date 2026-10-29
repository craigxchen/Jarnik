"""K12: syl8 (L=2) at e = 6 with the 98 shifts {+-e_x, +-e_x +- e_y} (mod p); and which cube shifts are needed:
greedy test of adding weight-3 shifts.  Own code."""
import sys, itertools, time
sys.path.insert(0, '.')
from klib import build_box, shiftset
from klib3 import shift_dims
exec(open('k6_zono.py').read().split("model = sys.argv[1]")[0])
A = build_box(8, [(S, 2) for S in syl8()])
d = 7
def unit(x, s): v = [0]*d; v[x] = s; return v
L98 = set()
for x in range(d):
    for s in (1, -1): L98.add(tuple(unit(x, s)))
for x in range(d):
    for y in range(x+1, d):
        for s in (1, -1):
            for t in (1, -1):
                v = [0]*d; v[x] = s; v[y] = t; L98.add(tuple(v))
L98 = sorted(L98)
hs, sm = shift_dims(A, L98, 6, 2147483629, d)
print("syl8 e=6 with %d shifts {+-e_x, +-e_x+-e_y}: H=%d dim sum=%d SHIFT=%s" % (len(L98), hs, sm, hs == sm), flush=True)
L3 = sorted(set(l for l in itertools.product((-1,0,1), repeat=d) if sum(abs(x) for x in l) == 3))
hs, sm = shift_dims(A, L98 + L3, 6, 2147483629, d)
print("syl8 e=6 with the 98 plus all %d weight-3 shifts: H=%d dim sum=%d SHIFT=%s" % (len(L3), hs, sm, hs == sm), flush=True)
