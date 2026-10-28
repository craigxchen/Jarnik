# list the chunks (D o c0 c1) not yet present in the output jsonl (resumability)
import sys, json, itertools, os
from math import comb
import common
D, CH, OUT = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
done = set()
if os.path.exists(OUT):
    for l in open(OUT):
        try: r = json.loads(l); done.add((r['orbit'], r['c0']))
        except Exception: pass
reps = common.orbit_reps(D, os.path.dirname(os.path.abspath(__file__)))
for o, A in enumerate(reps):
    n = comb(21 - len(set(A)) + D - 1, D)
    for c0 in range(0, n, CH):
        if (o, c0) not in done: print(D, o, c0, c0 + CH)
