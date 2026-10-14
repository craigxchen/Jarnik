# interpolation degree of double slices Lambda_{j,l} = {k in {-1,0,1}^(P+Q+R): s_P+s_Q=j, s_Q+s_R=l}; compare with d-|j|-|l|
import itertools
from tclib import hilbert, cube_monomials
bad = 0
for d in range(2, 8):
    mons = cube_monomials(d)
    for a in range(0, d+1):
        for b in range(0, d+1-a):
            c = d-a-b
            if a < c: continue
            P = range(0,a); Q = range(a,a+b); R = range(a+b,d)
            for j in range(0, a+b+1):
                for l in range(-(b+c), b+c+1):
                    pts = [k for k in itertools.product((-1,0,1), repeat=d) if sum(k[x] for x in P)+sum(k[x] for x in Q)==j and sum(k[x] for x in Q)+sum(k[x] for x in R)==l]
                    if not pts: continue
                    h = hilbert(pts, mons=mons); deg = max(h)
                    bound = d - abs(j) - abs(l)
                    if deg > bound:
                        bad += 1
                        if bad < 40: print('a=%d b=%d c=%d j=%d l=%d deg=%d bound=%d' % (a,b,c,j,l,deg,bound))
print('violations', bad)
