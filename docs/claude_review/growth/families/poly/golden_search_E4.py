import itertools, sys
from golden import template_best
best = {}
def upd(res, desc):
    for k, (v, n) in res.items():
        if k not in best or v < best[k][0] - 1e-9:
            best[k] = (v, desc, n)
patterns = []
S = int(sys.argv[1]) if len(sys.argv) > 1 else 10
for s2, s3, s4 in itertools.combinations(range(1, S+1), 3):
    patterns.append(([0, s2, s3, s4], [1, 1, 1, 1]))
for s2, s3 in itertools.combinations(range(1, S+1), 2):
    for w in ([2,1,1], [1,2,1], [1,1,2]):
        patterns.append(([0, s2, s3], w))
for s2 in range(1, S+1):
    for w in ([2,2], [3,1], [1,3]):
        patterns.append(([0, s2], w))
patterns.append(([0], [4]))
for shifts, widths in patterns:
    box = list(itertools.product(*[range(w+1) for w in widths]))
    for par in (0, 1):
        rows = [a for a in box if sum(a) % 2 == par]
        if len(rows) < 2: continue
        res = template_best(shifts, widths, rows, range(60, 120))
        upd(res, "shifts=%s widths=%s parity=%d rows=%d" % (shifts, widths, par, len(rows)))
for k in sorted(best):
    print("k=%d C=%.6f  %s  n=%d" % (k, best[k][0], best[k][1], best[k][2]))
