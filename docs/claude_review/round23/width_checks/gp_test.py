# (W) test on lattice points of signed Minkowski sums  P = sum_I y_I Delta_I - sum_J z_J Delta_J  (I,J subsets of [5]),
# translated into M (sum 0 after subtracting a fixed integer vector). reg mod p = certified upper bound for reg over Q,
# so "reg_p <= W" certifies (W) for that set; "reg_p > W" would need an exact rational certificate.
import itertools, random, sys
from reglib import reg
pairs=list(itertools.combinations(range(5),2))
subsets=[I for r in range(2,6) for I in itertools.combinations(range(5),r)]
def mink(A,B): return {tuple(a+b for a,b in zip(x,y)) for x in A for y in B}
def simplex(I,sgn):
    return {tuple(sgn*(1 if j==i else 0) for j in range(5)) for i in I}
random.seed(int(sys.argv[2])); N=int(sys.argv[1]); worst=(0,None)
for t in range(N):
    terms=[(random.choice(subsets),random.choice([1,-1])) for _ in range(random.randint(2,6))]
    P={(0,)*5}
    for I,s in terms: P=mink(P,simplex(I,s))
    # lattice points of conv(P) in the affine lattice: P is already the set of lattice points? Not in general;
    # use the vertex set's lattice hull: enumerate integer points of the box satisfying all 1_A-inequalities of P.
    hA={A:max(sum(x[i] for i in A) for x in P) for r in range(1,5) for A in itertools.combinations(range(5),r)}
    lA={A:min(sum(x[i] for i in A) for x in P) for A in hA}
    s0=sum(next(iter(P)))
    lo=[lA[(i,)] for i in range(5)]; hi=[hA[(i,)] for i in range(5)]
    S=[]
    for x in itertools.product(*[range(lo[i],hi[i]+1) for i in range(4)]):
        y=x+(s0-sum(x),)
        if all(lA[A]<=sum(y[i] for i in A)<=hA[A] for A in hA): S.append(y)
    if len(S)>900 or len(S)<2: continue
    W=sum(hA[A]-lA[A] for A in pairs)/6
    r=reg([y[:4] for y in S])
    if r/W>worst[0]: worst=(r/W,terms)
    print(t,len(S),'reg_p',r,'W',round(W,3),'ratio',round(r/W,3),flush=True)
print('worst',worst)
