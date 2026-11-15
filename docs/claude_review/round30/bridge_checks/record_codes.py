# Exponent codes of the record's integer clusters: cuts, weights, rank, comparison with the cube code.
import sys, math, itertools
sys.path.insert(0,'/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round30/bridge_checks')
from glib import gi_divmod, gi_gcd
from localgcd import split_prime
import sympy as sp

def vpi_int(x,pi):
    v=0
    while True:
        q,r=gi_divmod(x,pi)
        if r!=(0,0): return v
        x=q; v+=1

def code(points):
    pts=[tuple(p) for p in points]
    n=pts[0][0]**2+pts[0][1]**2
    assert all(x*x+y*y==n for x,y in pts) and len(set(pts))==len(pts)
    g=(0,0)
    for p in pts: g=gi_gcd(g,p) if g!=(0,0) else p
    assert g[0]**2+g[1]**2==1, ("not primitive",g)
    layers=[]   # (p, column tuple) one per threshold layer
    for p,e in sp.factorint(n).items():
        assert p%4==1
        a,b=split_prime(p); pi=(a,b)
        ex=[vpi_int(z,pi) for z in pts]
        for tau in range(1,e+1):
            layers.append((p,tuple(1 if x>=tau else 0 for x in ex)))
    return n,pts,layers

def angle_C(pts,n):
    R=math.sqrt(n); ang=sorted(math.atan2(y,x) for x,y in pts)
    span=ang[-1]-ang[0]
    return R*span/math.sqrt(R)

def analyse(name,points):
    n,pts,layers=code(points)
    M=len(pts)
    print(f"== {name}: n={n}={sp.factorint(n)}, M={M}, C={angle_C(pts,n):.6f}")
    cuts={}
    for p,col in layers:
        key=min(col,tuple(1-c for c in col))
        cuts[key]=cuts.get(key,0)+math.log(p)
    for k,w in cuts.items(): print("   cut",k," weight log =",round(w,3))
    vecs=[[c[j] for _,c in layers] for j in range(M)]
    rk=sp.Matrix([[vecs[j][k]-vecs[0][k] for k in range(len(layers))] for j in range(1,M)]).rank()
    print("   rank of difference lattice:",rk," #layers:",len(layers)," #distinct cuts:",len(cuts))
    # complementary-pair test: is there a perfect matching separated by every cut?
    ok=[]
    for perm in itertools.permutations(range(M)):
        pass
    sep=[[all(c[i]!=c[j] for c in cuts) for j in range(M)] for i in range(M)]
    def matchings(rem):
        if not rem: yield []; return
        i=rem[0]
        for j in rem[1:]:
            if sep[i][j]:
                for m in matchings([k for k in rem if k not in (i,j)]): yield [(i,j)]+m
    ms=list(matchings(list(range(M))))
    print("   perfect matchings separated by every cut:",ms[:3], "count",len(ms))
    return cuts,rk,ms

def axis_points(a,y0,depths):
    P=[]
    for d in depths:
        y2=y0*y0+2*d*a-d*d; y=math.isqrt(y2); assert y*y==y2
        P+= [(a-d,y),(a-d,-y)]
    return P

if __name__=="__main__":
    analyse("golden member n=537606725 (a=23186,y0=127)",axis_points(23186,127,(0,1,3,4)))
    analyse("golden member a=7465178,y0=2279",axis_points(7465178,2279,(0,1,3,4)))
    analyse("record 8-point C=5.674 (a=5743586,y0=527)",axis_points(5743586,527,(0,1,3,4)))
    analyse("8-point C=5.754 (a=590749322,y0=12809)",axis_points(590749322,12809,(0,1,3,4)))
    ten=[(24791,23712),(24745,23760),(24636,23873),(24567,23944),(24260,24255)]
    ten=ten+[(y,x) for x,y in ten]
    analyse("10-point C=8.2393",ten)
