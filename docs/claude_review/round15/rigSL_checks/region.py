# region.py NAME K depth : boxes of [0,pi/4]^r surviving the label-free constraints m_(v,s) >= 1 (B_box(e) >= 0)
import sys, math
from sl import *
name, K, depth = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
S = setup(name, K); r = S.r
def ok(lo, hi):
    return min(budgets_box(S, lo, hi)) >= 0
boxes = [([0.0] * r, [math.pi / 4] * r)]
for d in range(depth):
    nb = []
    for lo, hi in boxes:
        j = max(range(r), key=lambda k: hi[k] - lo[k]); mid = (lo[j] + hi[j]) / 2
        for a, b in ((lo[j], mid), (mid, hi[j])):
            l2 = list(lo); h2 = list(hi); l2[j] = a; h2[j] = b
            if ok(l2, h2): nb.append((l2, h2))
    boxes = nb
    vol = sum(math.prod(h - l for l, h in zip(lo, hi)) for lo, hi in boxes)
    if d % 4 == 3 or d == depth - 1:
        print(d + 1, len(boxes), f"vol={vol:.3e}", "hull:", [ (round(min(b[0][j] for b in boxes),4), round(max(b[1][j] for b in boxes),4)) for j in range(r)], flush=True)
    if not boxes: break
g = [abs(math.atan2(b, a)) % (math.pi/2) for (a, b) in S.pis]; g = [min(x, math.pi/2 - x) for x in g]
print("genuine angles folded to [0,pi/4]:", [round(x, 4) for x in g])
