# Search rational points y of the cube surface (sum y=0, prod y=1), all epsilon assignments,
# compute exact asymptotic template constant C_inf of the 8-point cluster (and best k-subclusters).
import sys, itertools, math, time
sys.path.insert(0,'/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round30/bridge_checks')
from localgcd import *
from fractions import Fraction as Fr

def roots_points(y,eps):
    r=[G(Fr(e)/v, v) for e,v in zip(eps,y)]
    roots=[]
    for x in r: roots+=[x,-x]
    if len(set(roots+[x.conj() for x in roots]))<16: return None
    sets=[[], list(range(8)), [0,1,2,3],[4,5,6,7],[0,1,4,5],[2,3,6,7],[0,2,4,6],[1,3,5,7]]
    return [[roots[k].conj() if k in F else roots[k] for k in range(8)] for F in sets]

def rat_points(H):
    sols=set()
    fr=sorted(set(Fr(p,q) for q in range(1,H+1) for p in range(-H,H+1) if p!=0))
    for a in fr:
        for b in fr:
            if b<a: continue
            # a b c^2 + a b (a+b) c + 1 = 0
            A=a*b; B=a*b*(a+b); Cc=Fr(1)
            disc=B*B-4*A*Cc
            if disc<0: continue
            n,d=disc.numerator,disc.denominator
            sn,sd=math.isqrt(n),math.isqrt(d)
            if sn*sn!=n or sd*sd!=d: continue
            s=Fr(sn,sd)
            for c in ((-B+s)/(2*A),(-B-s)/(2*A)):
                y=tuple(sorted((a,b,c,-(a+b+c))))
                if 0 in y: continue
                sols.add(y)
    return sorted(sols,key=lambda y:max(max(abs(v.numerator),v.denominator) for v in y))

if __name__=="__main__":
    H=int(sys.argv[1]) if len(sys.argv)>1 else 6
    t=time.time(); sols=rat_points(H); print("rational points (unordered):",len(sols),time.time()-t,flush=True)
    results=[]
    for y in sols:
        assert sum(y)==0 and y[0]*y[1]*y[2]*y[3]==1
        seen=set()
        for plus in itertools.combinations(range(4),2):
            minus=[k for k in range(4) if k not in plus]
            yy=[y[plus[0]],y[minus[0]],y[minus[1]],y[plus[1]]]
            key=frozenset(plus)
            if frozenset(minus) in seen: continue
            seen.add(key)
            pr=roots_points(yy,[1,-1,-1,1])
            if pr is None: continue
            out=C_inf(pr,8)
            if out is None: print("not a cluster?",yy); continue
            results.append((out[0],yy,out[1],out[2]))
            print(f"C_inf={out[0]:.6g}  |gamma|max={out[1]:.6g}  y={[str(v) for v in yy]}",flush=True)
    results.sort(key=lambda r:r[0])
    print("BEST:")
    for r in results[:8]: print(f"{r[0]:.8g}",[str(v) for v in r[1]],[(l,f) for (l,f,_) in r[3] if f!=1.0])
