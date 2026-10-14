from zono import *
from tclib import box_points, hilbert
import time
# M8 = Sylvester H8 L=2 (fano complement)
fano = [[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]]
cols = [[x for x in range(7) if x not in l] for l in fano]
for name, d, cols, Ls in [('hexagon', 2, [[0,1]], [1]), ('onecut s5L1', 5, [[0,1,2,3,4]], [1]), ('twocut 1,1,1', 3, [[0,1],[1,2]], [1,1]), ('M8', 7, cols, [2]*7)]:
    pts = box_points(d, list(zip(cols, Ls)))
    t = time.time()
    z, edges = zeta_lower(pts, cols, Ls, d)
    print(name, '|A|=%d zeta_lower=%d' % (len(pts), z), edges, '%.1fs' % (time.time()-t), flush=True)
