from faithful_search import *
import time
# Fano plane lines, d=7
fano = [[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]]
for name, d, cols in [('fano', 7, fano), ('fano-compl', 7, [[x for x in range(7) if x not in l] for l in fano])]:
    for L in (1, 2):
        t = time.time()
        r = analyze(d, cols, [L]*len(cols))
        print(name, 'L=%d' % L, r, '%.1fs' % (time.time()-t), flush=True)
