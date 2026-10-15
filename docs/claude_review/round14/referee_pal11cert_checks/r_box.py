# Referee: independent construction of pal11 box, patterns, widths of test directions.
import itertools, numpy as np, pickle, sys
q=11; M=12
def leg(a):
    a%=q
    return 0 if a==0 else (1 if pow(a,(q-1)//2,q)==1 else -1)
# Paley-I matrix (as in round-5 construct / fplib.paley_I): rows/cols indexed inf,0..10
idx=['inf']+list(range(q))
H=np.zeros((M,M),int)
for a,x in enumerate(idx):
    for b,y in enumerate(idx):
        if x=='inf': H[a,b]=1
        elif y=='inf': H[a,b]=-1
        else: H[a,b]=(1 if x==y else 0)+leg(y-x)
assert (H@H.T==12*np.eye(12)).all()
Hn=H.copy(); Hn[1:,:]*=-1
def canon(S):
    S=frozenset(S)
    if 0 in S: S=frozenset(range(M))-S
    return tuple(sorted(S))
pats=set()
for a in range(1,M):
    S=[x for x in range(M) if Hn[x,a]==1]
    if 0<len(S)<M: pats.add(canon(S))
for x in range(M):
    c=Hn[:,x if x>0 else 1].copy()   # lab[0]=1, lab[x]=x
    c[x]*=-1
    S=[y for y in range(M) if c[y]==1]
    if 0<len(S)<M: pats.add(canon(S))
pats=sorted(pats)
print('npat',len(pats),'sizes',sorted(len(P) for P in pats))
core=[P for P in pats if len(P)==6]; twin=[P for P in pats if len(P)==5]
print('core',len(core),'twin',len(twin))
# A
d=11
allk=np.array(list(itertools.product((-1,0,1),repeat=d)),dtype=np.int8)
ok=np.ones(len(allk),bool)
for P in pats:
    s=allk[:,[x-1 for x in P]].sum(1)
    ok&=np.abs(s)<=2
A=allk[ok]
print('|A|',len(A))
k0=-A.sum(1).astype(int); print('base range',k0.min(),k0.max())
for P in pats:
    s=A[:,[x-1 for x in P]].sum(1); assert s.min()==-2 and s.max()==2
# Hadamard columns on 12 rows: all-ones and beta_P=2*1_P-1 (P core, base row 0 not in P)
B=[np.ones(12,int)]+[np.array([1 if r in P else -1 for r in range(12)]) for P in core]
B=np.array(B)
print('Hadamard check', (B@B.T==12*np.eye(12)).all())
# directions: {-1,0,1}^11 up to sign, plus {-2..2}^11 support<=4 containing a +-2, primitive, up to sign
dirs=[]
for v in itertools.product((-1,0,1),repeat=d):
    nz=[t for t in v if t]
    if nz and nz[0]>0: dirs.append(v)
from math import gcd
for s in range(1,5):
    for supp in itertools.combinations(range(d),s):
        for vals in itertools.product((-2,-1,1,2),repeat=s):
            if vals[0]<0: continue
            if 2 not in [abs(t) for t in vals]: continue
            g=0
            for t in vals: g=gcd(g,abs(t))
            if g!=1: continue
            v=[0]*d
            for i,t in zip(supp,vals): v[i]=t
            dirs.append(tuple(v))
D=np.array(dirs,dtype=np.int8); print('ndirs',len(D))
Af=A.astype(np.float32)
tau=np.zeros(len(D),int)
for i in range(0,len(D),4000):
    V=D[i:i+4000].astype(np.float32)@Af.T
    V=np.rint(V).astype(int)+40
    for j in range(V.shape[0]):
        tau[i+j]=np.count_nonzero(np.bincount(V[j],minlength=81))-1
sel=tau<8
print('n tau<8',sel.sum(), 'dist',np.unique(tau[sel],return_counts=True))
pickle.dump(dict(pats=pats,core=core,twin=twin,A=A,D=D[sel],tau=tau[sel],H=H,B=B),open('r_box.pkl','wb'))
