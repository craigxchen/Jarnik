import sys, itertools
sys.path.insert(0,'.')
from localgcd import *
from cube_search import roots_points
from fractions import Fraction as Fr
pts=[]
for line in open(sys.argv[1]):
    h,*ys=line.split(); pts.append([Fr(v) for v in ys])
res=[]; done=set()
for y in pts:
    key=tuple(sorted(y)); mkey=tuple(sorted(-v for v in y))
    if key in done or mkey in done: continue
    done.add(key)
    for plus in itertools.combinations(range(4),2):
        minus=tuple(k for k in range(4) if k not in plus)
        yy=[y[plus[0]],y[minus[0]],y[minus[1]],y[plus[1]]]
        pr=roots_points(yy,[1,-1,-1,1])
        if pr is None: continue
        out=C_inf(pr,8)
        res.append((out[0],yy,out[1],[(l,f) for (l,f,_) in out[2] if f!=1.0]))
        print(f"{out[0]:.6g} gam={out[1]:.5g} y={' '.join(map(str,yy))} loc={res[-1][3]}",flush=True)
        break
res.sort(key=lambda r:r[0])
print("BEST")
for r in res[:10]: print(f"{r[0]:.8g}",' '.join(map(str,r[1])),r[3])
