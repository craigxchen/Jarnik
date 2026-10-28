# (1) Picard ranks (Lemma 2.3); (2) Prop 5.3 identity and per-section comparison; (3) width part of pullback Lemma R1.
import itertools, random
random.seed(7)
rays=[frozenset(A) for r in range(1,5) for A in itertools.combinations(range(5),r)]
orb={frozenset([A,frozenset(range(5))-A]) for A in rays}
print("(1) rays",len(rays),"Galois orbits",len(orb),"-> rank Pic twisted =",len(orb)-0,", split =",len(rays)-4)
def rS():
    pts=[s for s in itertools.product(range(-2,3),repeat=4) if sum(abs(a) for a in s)+abs(sum(s))<=4]
    return [s+(-sum(s),) for s in random.sample(pts,random.randint(2,6))]
def ip(s,v): return sum(a*b for a,b in zip(s,v))
def rsp(Q,v): return -min(ip(s,v) for s in Q)
def wid(Q,v): return max(ip(s,v) for s in Q)-min(ip(s,v) for s in Q)
cnt=0; stronger=0
for _ in range(20000):
    Q=rS(); v=tuple(random.randint(-2,2) for _ in range(5))
    assert 2*wid(Q,v)/2 == (rsp(Q,v)+rsp(Q,tuple(-a for a in v)))
    cnt+=1; stronger+= (wid(Q,v)/2 < rsp(Q,v))
print("(2) identity r_tw = (r_sp(v)+r_sp(-v))/2 holds on",cnt,"random (Q,v); twisted RHS strictly smaller than split RHS at v in",stronger,"cases")
Q=[(0,)*5,(1,-1,0,0,0)]; v=(0,1,0,0,0)
print("    example f=x1/x2-1, v=1_{2}: r_sp(v)=",rsp(Q,v)," r_tw=",wid(Q,v)/2)
# (3) phi(y1,y2)=(y1, y1*y2): phi^*(s,s') = (s+s', s'); phi_*(v1,v2)=(v1, v1+v2)
bad=0
for _ in range(5000):
    QF=[(s,t) for s,t in zip(rS(),rS())]
    QG=list({tuple(a+b for a,b in zip(s,t))+t for s,t in QF})   # support of pullback (before cancellation) = phi^*(supp F)
    QGc=random.sample(QG,max(1,len(QG)-1))                    # cancellation can only remove points
    v1=tuple(random.randint(-2,2) for _ in range(5)); v2=tuple(random.randint(-2,2) for _ in range(5))
    pv=v1+tuple(a+b for a,b in zip(v1,v2))
    wF=wid([s+t for s,t in QF],pv); wG=wid(QG,v1+v2); wGc=wid(QGc,v1+v2)
    bad+= not (wF==wG and wGc<=wG)
print("(3) width_{Q_F}(phi_* v) = width_{phi^*Q_F}(v) >= width_{Q_G}(v): failures",bad,"of 5000")
