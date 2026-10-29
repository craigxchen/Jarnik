"""K14: M8' (syl8, L = 3) confirmation: (SHIFT_e) w.r.t. L_pair at e = 8, 9, 10, 11 with a second prime; exact check
of the order-10 zonotope; syl8 L = 1 reach / zeta / rho."""
import sys, itertools, time
sys.path.insert(0, '.')
from klib import build_box, Lpair, monos_lt, mval
from klib3 import shift_dims
exec(open('k6_zono.py').read().split("model = sys.argv[1]")[0])
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
A = build_box(8, [(S, 3) for S in syl8()])
for e in (8, 9, 10, 11):
    hs, sm = shift_dims(A, Lpair(7), e, 1000000007, 7)
    pr("M8' p=1000000007 e=%d H=%d dim sum=%d SHIFT=%s" % (e, hs, sm, hs == sm))
k0 = (-1, 0, 0, -1, 0, -1, 0)
gens = [(0,0,0,0,0,0,1), (0,0,0,0,0,1,-1), (0,0,0,0,0,1,0), (0,0,0,0,1,0,0), (0,0,0,1,-1,0,0), (0,0,0,1,0,0,0),
        (0,0,1,0,0,0,0), (0,1,0,0,0,0,0), (1,-1,-1,0,0,0,0), (1,0,0,0,0,0,0)]
Aset = set(A)
c = {}
for t in itertools.product((0, 1), repeat=10):
    p = tuple(k0[x] + sum(t[i] * gens[i][x] for i in range(10)) for x in range(7))
    assert p in Aset
    c[p] = c.get(p, 0) + (-1) ** (10 - sum(t))
c = {k: v for k, v in c.items() if v}
ok10 = all(sum(v * mval(k, a) for k, v in c.items()) == 0 for a in monos_lt(7, 10))
ok11 = all(sum(v * mval(k, a) for k, v in c.items()) == 0 for a in monos_lt(7, 11))
pr("M8' order-10 zonotope: all 1024 points in A: True; c in V(A,10): %s; c in V(A,11): %s; |c|_inf=%d" % (ok10, ok11, max(abs(v) for v in c.values())))
A1 = build_box(8, [(S, 1) for S in syl8()])
H1 = []; reach = None
for e in range(0, 10):
    hs, sm = shift_dims(A1, Lpair(7), e, 2147483629, 7)
    H1.append(hs)
    if hs == 0: break
    if reach is None and hs != sm: reach = e
pr("syl8 L=1: |A|=%d H=%s reach=%s" % (len(A1), H1, reach))
open("k14_M8p_confirm_out.txt", "w").write("\n".join(out) + "\n")
