# Independent check of theory.md Lemma 2.2(1) and identity (1.1) (referee code).
import random, itertools, math
from sympy import factorint, Matrix
from flint import fmpz_mat
random.seed(12345)
n=7
def contacts(x):
    # d_T for visible T (2<=|T|<=6), via residue classes mod p^t
    diffs=[x[i]-x[j] for i in range(n) for j in range(i+1,n)]
    primes=set()
    for d in diffs: primes|=set(factorint(abs(d)).keys())
    d={}
    classes_all=[]  # (p,t,class)
    for p in primes:
        t=1
        while True:
            cl={}
            for i in range(n): cl.setdefault(x[i]%p**t,[]).append(i)
            cls=[frozenset(c) for c in cl.values() if len(c)>=2]
            if not cls: break
            for c in cls:
                if 2<=len(c)<=6: d[c]=d.get(c,1)*p
                if len(c)==7: raise Exception('not primitive')
                classes_all.append((p,t,c))
            t+=1
    return d,classes_all
def index_lattice(x,classes_all):
    rows=[];mods=[]
    for (p,t,c) in classes_all:
        if len(c)<3: continue
        c=sorted(c); i0=c[0]
        for i in c[1:]:
            r=[0]*n; r[i]=1; r[i0]=-1; rows.append(r); mods.append(p**t)
    if not rows: return 1
    r=len(rows)
    # image of Z^n in prod Z/m_k ; |image| = prod m / [Z^r : L], L = colspan(M) + diag(m)
    big=[[rows[k][j] for j in range(n)]+[mods[k] if kk==k else 0 for kk in range(r)] for k in range(r)]
    Mt=fmpz_mat(big)  # r x (n+r), columns generate L
    H=Mt.transpose().hnf()  # rows generate L
    det=1
    Hl=H.tolist()
    # HNF rows: first r nonzero rows form upper triangular basis
    nz=[row for row in Hl if any(v!=0 for v in row)]
    assert len(nz)==r
    for k in range(r): det*=int(nz[k][k])
    prodm=1
    for m in mods: prodm*=m
    assert prodm%abs(det)==0
    return prodm//abs(det)
ok=0;tot=0;bad=[]
for trial in range(80):
    while True:
        mode=trial%3
        if mode==0: x=[random.randint(0,4000) for _ in range(n)]
        elif mode==1:
            base=random.randint(0,50); x=[base+ (2**random.randint(0,4))*(3**random.randint(0,3))*random.randint(0,30) for _ in range(n)]
        else:
            x=[random.randint(0,3)*360+random.randint(0,2)*60+random.randint(0,1)*7*11 for _ in range(n)]
        if len(set(x))<n: continue
        g=0
        for i in range(n):
            for j in range(i+1,n): g=math.gcd(g,x[i]-x[j])
        x=[(v-x[0])//g for v in x]
        break
    d,cla=contacts(x)
    # (1.1)
    good=True
    for i in range(n):
        for j in range(i+1,n):
            pr=1
            for T,v in d.items():
                if i in T and j in T: pr*=v
            if pr!=abs(x[i]-x[j]): good=False
    idx=index_lattice(x,cla)
    pred=1
    for T,v in d.items():
        if len(T)>=3: pred*=v**(len(T)-1)
    tot+=1
    if good and idx==pred: ok+=1
    else: bad.append((x,idx,pred,good))
print('configs',tot,'pass',ok); print('bad',bad[:3])
