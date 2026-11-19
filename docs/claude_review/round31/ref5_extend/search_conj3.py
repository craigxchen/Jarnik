# Char 3, "twisted" design: Y = T real (T in F_3[t]), factors T + w with w in F_9[t] non-real; conj(T+w) = T + wbar.
# Points: transversals of pairs {w, wbar} of a conj-stable set W (full transversals, all have norm prod_W (T+w)).
# Group by top coefficients above floor(D/2) after monic normalisation. Random W of npairs pairs from F_9[t]_{<e}.
import sys, random, itertools
from search_trans3 import MUL, ADD, INV, pmul
def conj(c): a,b=c%3,c//3; return a+3*((-b)%3)
random.seed(int(sys.argv[1])); trials=int(sys.argv[2]); best={}
allw=[w for w in itertools.product(range(9),repeat=2) if any(c//3 for c in w)]  # F_9[t]_{<2}, non-real
reps=[];seen=set()
for w in allw:
    if w in seen: continue
    wb=tuple(conj(c) for c in w); seen|={w,wb}; reps.append((w,wb))
for tr in range(trials):
    npairs=random.randint(4,12); P=random.sample(reps,npairs)
    L=random.randint(2,14); T=[0]*L+[1]
    if random.random()<0.3: T=[random.randrange(3) for _ in range(L)]+[1]
    def fac(w):
        f=list(T)+[0]*max(0,len(w)-len(T))
        for k,c in enumerate(w): f[k]=ADD[f[k]][c]
        while len(f)>1 and f[-1]==0: f.pop()
        return f
    F=[(fac(a),fac(b)) for a,b in P]
    groups={}
    def dfs(k,prod):
        if k==npairs:
            lc=INV[prod[-1]]; p=tuple(MUL[lc][c] for c in prod); D=len(p)-1; K=D//2
            groups.setdefault(p[K+1:],set()).add(p); return
        dfs(k+1,pmul(prod,F[k][0])); dfs(k+1,pmul(prod,F[k][1]))
    dfs(0,[1])
    for key,ps in groups.items():
        m=len(ps)
        if m>=3:
            D=len(next(iter(ps)))-1
            if m not in best or D<best[m][0]: best[m]=(D,T,P)
for m in sorted(best): print(m,best[m][0],"T=",best[m][1])
