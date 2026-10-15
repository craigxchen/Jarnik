# Referee: pal7 at L (n=8) prep: pool support<=4 entries +-1, dirs {-2..2}^7 with tau<m.
import itertools, numpy as np, sys
L=int(sys.argv[1]); m=int(sys.argv[2]); pre=sys.argv[3]
q=7; M=8
def leg(a):
    a%=q; return 0 if a==0 else (1 if pow(a,(q-1)//2,q)==1 else -1)
idx=['inf']+list(range(q)); H=np.zeros((M,M),int)
for a,x in enumerate(idx):
    for b,y in enumerate(idx):
        H[a,b]=1 if x=='inf' else (-1 if y=='inf' else (1 if x==y else 0)+leg(y-x))
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
    c=Hn[:,x if x>0 else 1].copy(); c[x]*=-1
    S=[y for y in range(M) if c[y]==1]
    if 0<len(S)<M: pats.add(canon(S))
allk=np.array(list(itertools.product((-1,0,1),repeat=7)))
ok=np.ones(len(allk),bool)
for P in pats: ok&=np.abs(allk[:,[x-1 for x in P]].sum(1))<=L
A=allk[ok]; print('pal7 L',L,'|A|',len(A),'tau(1)',len(set(A.sum(1)))-1)
D=np.array([v for v in itertools.product(range(-2,3),repeat=7) if any(v) and next(x for x in v if x)>0])
V=A@D.T; tau=np.array([len(set(V[:,j]))-1 for j in range(V.shape[1])]); keep=tau<m; D=D[keep]; tau=tau[keep]
print('dirs',len(D))
pool=[]
for s in range(2,5):
    for S in itertools.combinations(range(M),s):
        for vals in itertools.product((-1,1),repeat=s):
            if sum(vals)!=0 or vals[0]<0: continue
            z=[0]*M
            for x,v in zip(S,vals): z[x]=v
            pool.append(z)
with open(pre+'_dirs.txt','w') as f:
    f.write('%d\n'%len(D))
    for v,t in zip(D,tau): f.write(' '.join(map(str,v))+' %d\n'%t)
with open(pre+'_pool.txt','w') as f:
    f.write('%d\n'%len(pool))
    for z in pool: f.write(' '.join(map(str,z))+'\n')
np.save(pre+'_A.npy',A)
