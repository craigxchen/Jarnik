# Referee: verify the claimed length-10 zonotope in M8' = pal7 at L=3 (d=7), and tau_A(1) for the models.
import itertools, numpy as np
def build(q,L):
    M=q+1
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
    d=q
    allk=np.array(list(itertools.product((-1,0,1),repeat=d)))
    ok=np.ones(len(allk),bool)
    for P in pats: ok&=np.abs(allk[:,[x-1 for x in P]].sum(1))<=L
    return allk[ok],pats
A,pats=build(7,3); print("M8' |A|",len(A),'npats',len(pats))
Aset=set(map(tuple,A))
print("tau_A(1) M8'",len(set(A.sum(1)))-1)
e=lambda i: np.eye(7,dtype=int)[i-1]
Z=[-e(1),-e(1),-e(2),-e(2),-e(3),-e(4),-e(5),-e(6),e(3)-e(4),e(5)-e(6)]
sums=set()
for S in itertools.product((0,1),repeat=10):
    sums.add(tuple(sum(s*z for s,z in zip(S,Z))))
sums=[np.array(s) for s in sums]
good=[k0 for k0 in A if all(tuple(k0+s) in Aset for s in sums)]
print('translates k0 with k0+Z subset A:',len(good), good[:3])
A2,_=build(7,2); print('pal7 L2 |A|',len(A2),'tau(1)',len(set(A2.sum(1)))-1)
print('--- 0-indexed reading e_i = coordinate i (0..6)')
e0=lambda i: np.eye(7,dtype=int)[i]
Z=[-e0(1),-e0(1),-e0(2),-e0(2),-e0(3),-e0(4),-e0(5),-e0(6),e0(3)-e0(4),e0(5)-e0(6)]
sums=set(tuple(sum(s*z for s,z in zip(S,Z))) for S in itertools.product((0,1),repeat=10))
good=[k0 for k0 in A if all(tuple(k0+np.array(s)) in Aset for s in sums)]
print('translates:',len(good), [tuple(g) for g in good[:3]])
