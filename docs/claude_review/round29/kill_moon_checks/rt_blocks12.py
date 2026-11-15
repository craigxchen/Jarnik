"""kill_moon: Sperner blocks at D=12, c=0 (s=6, m=5): four disjoint triples in H with equal
(sum Im a^r)_{r<=5}. Numerical (mpmath findroot, 50 digits) -> evidence only; then verify the 6 points' top
coefficients (degrees 7..11) agree."""
import mpmath as mp, random, itertools
mp.mp.dps=50
random.seed(1)
def F(blk): return [sum(mp.im(mp.mpc(x,y)**r) for x,y in blk) for r in range(1,6)]
base=[(0,1),(1,mp.mpf('1.7')),(mp.mpf('-0.6'),mp.mpf('2.3'))]
tau=F(base)
blocks=[base]
tries=0
while len(blocks)<4 and tries<400:
    tries+=1
    x1=mp.mpf(random.uniform(-2.5,2.5))
    st=[random.uniform(0.3,3), random.uniform(-2,2),random.uniform(0.3,3),random.uniform(-2,2),random.uniform(0.3,3)]
    def eqs(y1,x2,y2,x3,y3):
        return [a-b for a,b in zip(F([(x1,y1),(x2,y2),(x3,y3)]),tau)]
    try:
        sol=mp.findroot(eqs,st,tol=mp.mpf(10)**-40,maxsteps=200)
    except Exception: continue
    y1,x2,y2,x3,y3=[sol[i] for i in range(5)]
    blk=[(x1,y1),(x2,y2),(x3,y3)]
    if min(y for _,y in blk)<0.05: continue
    allpts=[mp.mpc(x,y) for b in blocks+[blk] for x,y in b]
    if min(abs(p-q) for p,q in itertools.combinations(allpts,2))<1e-3: continue
    if max(abs(v) for v in eqs(*[sol[i] for i in range(5)]))>mp.mpf(10)**-35: continue
    blocks.append(blk)
print("found blocks:",len(blocks),"after",tries,"tries")
for b in blocks: print([ (mp.nstr(x,10),mp.nstr(y,10)) for x,y in b])
def pr(rs):
    c=[mp.mpc(1)]
    for r in rs:
        c=[mp.mpc(0)]+c
        for k in range(len(c)-1): c[k]-=r*c[k+1]
    return c
if len(blocks)==4:
    pts=[]
    for S in itertools.combinations(range(4),2):
        rs=[]
        for j,b in enumerate(blocks):
            for x,y in b: rs.append(mp.mpc(x,y) if j in S else mp.mpc(x,-y))
        pts.append(pr(rs))
    worst=max(abs(z[k]-w[k]) for z,w in itertools.combinations(pts,2) for k in range(7,12))
    low=min(max(abs(z[k]-w[k]) for k in range(0,7)) for z,w in itertools.combinations(pts,2))
    print("max |coef diff| in degrees 7..11:",mp.nstr(worst,5),"; min over pairs of max diff in degrees<=6:",mp.nstr(low,5))
