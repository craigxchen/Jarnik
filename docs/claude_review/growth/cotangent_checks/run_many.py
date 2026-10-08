import sys, math
from search_pruned import run
Amax = int(sys.argv[1]); m = int(sys.argv[2]); Cmax = float(sys.argv[3])
Ls = [int(x) for x in sys.argv[4].split(',')]
for L in Ls:
    f = run(L, Amax, m, Cmax, verbose=False)
    big = [r for r in f if r[1] > 10**6]
    best = f[0] if f else None
    bestbig = big[0] if big else None
    print('L=%d count=%d best=%s bestN>1e6=%s' % (L, len(f), best, bestbig), flush=True)
