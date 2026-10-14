# rho(A) for two-cut blocks A = {k in {-1,0,1}^(a+b+c): |s_P+s_Q|<=L1, |s_Q+s_R|<=L2}  (mod-p Hilbert function)
import itertools, sys
from tclib import *
out = []
for d in range(2, 8):
    mons = cube_monomials(d)
    for a in range(0, d+1):
        for b in range(0, d+1-a):
            c = d - a - b
            if a + b == 0 or b + c == 0: continue
            if a < c: continue  # symmetry P<->R
            for L1 in range(0, a+b):
                for L2 in range(0, b+c):
                    if (c, b, a, L2, L1) < (a, b, c, L1, L2) and a == c: continue
                    S1 = list(range(0, a+b)); S2 = list(range(a, a+b+c))
                    pts = box_points(d, [(S1, L1), (S2, L2)])
                    h = hilbert(pts, mons=mons)
                    rho = max(h)
                    line = 'a=%d b=%d c=%d L1=%d L2=%d |A|=%d rho=%d  d+L1+L2=%d  H=%s' % (a,b,c,L1,L2,len(pts),rho,d+L1+L2,[h.get(i,0) for i in range(rho+1)])
                    print(line, flush=True)
