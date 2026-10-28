import json, sys, math, os
D, CH, OUT = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
W='/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round19/referee_search_checks'
reps = json.load(open(W+'/rreps_D%d.json' % D))
done = set()
if os.path.exists(OUT):
    for l in open(OUT):
        try:
            r = json.loads(l); done.add((r['orbit'], r['c0']))
        except Exception: pass
for o, A in enumerate(reps):
    s = len(set(tuple(e) for e in A))
    n = math.comb(21 - s + D - 1, D)
    for c0 in range(0, n, CH):
        if (o, c0) not in done:
            print(D, o, c0, min(c0 + CH, n))
