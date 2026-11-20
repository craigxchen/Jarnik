# Judge check for Part VIII: square-block identity (VIII.3), divisor columns (VIII.5),
# moduli-point equality (escape hypothesis), and Theorem D admissibility via (3.2).
import random, math, itertools, sys
from fractions import Fraction
from math import gcd
random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 1)
PR=[2,3,5,7,11,13,17,19,23]
def vp(n,p):
    n=abs(n); 
    if n==0: return 10**9
    c=0
    while n%p==0: n//=p; c+=1
    return c
def det(P,Q): return P[0]*Q[1]-P[1]*Q[0]
def crossratio(a,b,c,d):  # of points of P^1 given as Fractions or columns
    return Fraction((a-c)*(b-d),(a-d)*(b-c))
def colpt(P): return Fraction(P[1],P[0])  # [X:Y] -> Y/X
def lT(g,T,idx):
    best=None
    Tc=[i for i in idx if i not in T]
    for a,b in itertools.combinations(T,2):
        for c,d in itertools.combinations(Tc,2):
            v=min(g[a][b]+g[c][d]-g[a][c]-g[b][d], g[a][b]+g[c][d]-g[a][d]-g[b][c])
            best=v if best is None else min(best,v)
    return max(0,best)
def chain_ok(sets):
    # sets: list of frozensets (with multiplicity irrelevant); check partition into two chains with disjoint unions
    S=list(set(sets))
    for mask in range(1<<len(S)):
        A=[S[i] for i in range(len(S)) if mask>>i&1]; B=[S[i] for i in range(len(S)) if not mask>>i&1]
        def ischain(C): return all(x<=y or y<=x for x,y in itertools.combinations(C,2))
        if ischain(A) and ischain(B):
            uA=frozenset().union(*A) if A else frozenset(); uB=frozenset().union(*B) if B else frozenset()
            if not (uA & uB): return True
    return False
stats=dict(tuples=0,pairs=0,cr=0,lT=0,lTdiff=0)
def check_data(P,n,marks,label):
    # P: 8 columns; n: dict T->n_T ; marks: 8 rationals
    idx=list(range(8))
    assert P[0]==(1,0)
    for i in range(8): assert gcd(P[i][0],P[i][1])==1, (label,'prim')
    t={}
    for i,j in itertools.combinations(idx,2):
        blk=1
        for T,v in n.items():
            if i in T and j in T: blk*=v
        d=det(P[i],P[j]); assert d%blk==0, (label,'div'); t[(i,j)]=d//blk; assert t[(i,j)]!=0
        if i==0: assert t[(0,j)]==P[j][1]
        stats['pairs']+=1
    # moduli point: all cross-ratios of columns equal those of marks (common Mobius)
    pts=[colpt(Pi) for Pi in P]
    for q4 in itertools.combinations(idx,4):
        a,b,c,d=q4
        # column 0 is (1,0) -> point 0 under Y/X; fine (finite)
        assert crossratio(*[pts[k] for k in q4])==crossratio(*[marks[k] for k in q4]),(label,'cr')
    stats['cr']+=1
    # Theorem D (3.2): |l_T(p)-v_p(n_T)| <= 2 max v_p(t_ij), at every prime of the data
    primes=set()
    for v in list(n.values())+[abs(x) for x in t.values()]:
        for p in PR+[29,31,37,41,43,47]:
            if v%p==0: primes.add(p)
    HS=sum(math.log(abs(x)) for x in t.values())
    for p in primes:
        g=[[0]*8 for _ in range(8)]
        for i,j in itertools.combinations(idx,2):
            g[i][j]=g[j][i]=vp(det(P[i],P[j]),p)
        mt=max(vp(x,p) for x in t.values())
        for r in range(2,7):
            for T in itertools.combinations(range(1,8),r):
                l=lT(g,T,idx); nv=vp(n.get(frozenset(T),1),p)
                assert abs(l-nv)<=2*mt,(label,'3.2',p,T,l,nv,mt)
                stats['lT']+=1
                if l!=nv: stats['lTdiff']+=1
    stats['tuples']+=1
# ---------- VIII.2 square blocks ----------
def rand_two_class():
    ps=random.sample(PR[:7],random.randint(1,4)); q=1; res={}
    for p in ps:
        e=random.randint(1,4 if p<5 else 3); pe=p**e; q*=pe
        mode=random.random()
        r1=random.randrange(pe)
        if mode<0.2: r2=r1
        elif mode<0.5 and e>1:
            s=random.randint(1,e-1); r2=(r1+p**s*random.randrange(1,p**(e-s)))%pe
        else: r2=random.randrange(pe)
        res[pe]=(r1,r2)
    return q,res
def crt(rs):
    x,m=0,1
    for pe,r in rs:
        # solve x = r mod pe
        k=((r-x)*pow(m,-1,pe))%pe; x+=m*k; m*=pe
    return x
nsq=0
for inst in range(120):
    q,res=rand_two_class()
    X=set()
    for _ in range(40):
        rs=[(pe,random.choice(res[pe])) for pe in res]
        X.add(crt(rs)+q*random.randint(-3,3))
    X=sorted(X)
    if len(X)<8: continue
    xs=X[0]; g=q
    for x in X: g=gcd(g,x-xs)
    Q=q//g; Y=[(x-xs)//g for x in X]
    fac={}
    QQ=Q
    for p in PR:
        f=vp(QQ,p) if QQ>1 else 0
        if f and f<10**8: fac[p]=f
    # two classes with unit difference at each p^f || Q
    for p,f in fac.items():
        cl=set(y%(p**f) for y in Y); assert len(cl)==2; a,b=cl; assert (a-b)%p!=0
    for _ in range(6):
        y=random.sample(Y,8)
        d=[gcd(Q,yi-y[0]) for yi in y]; d[0]=Q
        P=[(Q//d[i],(y[i]-y[0])//d[i]) for i in range(8)]
        n={}
        for p,f in fac.items():
            T=frozenset(i for i in range(1,8) if (y[i]-y[0])%(p**f)!=0)
            if T: n[T]=n.get(T,1)*p**(2*f)
        # cofactor identity t = (y_j-y_i)/h
        for i,j in itertools.combinations(range(8),2):
            h=gcd(Q,y[j]-y[i]); D2=Q//h  # e^{D/2}
            blk=1
            for T,v in n.items():
                if i in T and j in T: blk*=v
            assert det(P[i],P[j])==((y[j]-y[i])//h)*blk
        check_data(P,n,[Fraction(v) for v in y],'sq'); nsq+=1
print('square-block tuples',nsq,stats)
# ---------- VIII.3 divisors ----------
ndv=0
for inst in range(200):
    q=random.choice([3,4,5,7,8,9,11,13])
    ps=random.sample([p for p in PR if q%p] ,random.randint(2,4))
    if random.random()<0.3:
        ps.append(next(p for p in PR if q%p==0))
    nn=1
    for p in ps: nn*=p**random.randint(1,4)
    divs=[dd for dd in range(1,nn+1) if nn%dd==0]
    classes={}
    for dd in divs:
        if gcd(dd,q)==1: classes.setdefault(dd%q,[]).append(dd)
    big=[c for c in classes.values() if len(c)>=8]
    if not big: continue
    D=random.choice(big)
    for _ in range(4):
        dd=random.sample(D,8); a=dd[0]
        gg=[gcd(a,x) for x in dd]
        P=[(a//gg[i],(dd[i]-a)//(q*gg[i])) for i in range(8)]
        assert all((dd[i]-a)%(q*gg[i])==0 for i in range(8))
        n={}
        for p in set(ps):
            e=vp(nn,p); b=[vp(x,p) for x in dd]
            for s in range(1,e+1):
                T=frozenset(i for i in range(1,8) if (s<=b[i])!=(s<=b[0]))
                if T: n[T]=n.get(T,1)*p
        sets=[]
        for p in set(ps):
            sp=[T for T,v in n.items() if v%p==0]
            assert chain_ok(sp)
        check_data(P,n,[Fraction(x) for x in dd],'div'); ndv+=1
print('divisor tuples',ndv,stats)
