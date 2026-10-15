# Referee: M8 = syl8 at L=2 (n=8): pool (support<=4, entries +-1, balanced, primitive, sign-normalised), test dirs {-2..2}^7 with tau<m.
import itertools, numpy as np, sys
L=int(sys.argv[1]); m=int(sys.argv[2]); pre=sys.argv[3]
M=8; pats=set()
for a in range(1,M):
    S=frozenset(x for x in range(M) if bin(a&x).count('1')%2==1)
    if 0 in S: S=frozenset(range(M))-S
    pats.add(tuple(sorted(S)))
allk=np.array(list(itertools.product((-1,0,1),repeat=7)))
ok=np.ones(len(allk),bool)
for P in pats: ok&=np.abs(allk[:,[x-1 for x in P]].sum(1))<=L
A=allk[ok]; print('|A|',len(A),'tau(1)',len(set(A.sum(1)))-1)
D=[v for v in itertools.product(range(-2,3),repeat=7) if any(v) and next(x for x in v if x)>0]
D=np.array(D); V=A@D.T; tau=np.array([len(set(V[:,j]))-1 for j in range(V.shape[1])])
keep=tau<m; D=D[keep]; tau=tau[keep]; print('dirs tau<m',len(D))
pool=[]
for s in range(2,5):
    for S in itertools.combinations(range(M),s):
        for vals in itertools.product((-1,1),repeat=s):
            if sum(vals)!=0 or vals[0]<0: continue
            z=[0]*M
            for x,v in zip(S,vals): z[x]=v
            pool.append(z)
print('pool',len(pool))
with open(pre+'_dirs.txt','w') as f:
    f.write('%d\n'%len(D))
    for v,t in zip(D,tau): f.write(' '.join(map(str,v))+' %d\n'%t)
with open(pre+'_pool.txt','w') as f:
    f.write('%d\n'%len(pool))
    for z in pool: f.write(' '.join(map(str,z))+'\n')
np.save(pre+'_A.npy',A)
