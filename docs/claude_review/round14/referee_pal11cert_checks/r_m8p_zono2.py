# Referee: M8' as syl8 (Sylvester n=8) at L=3, M8 = syl8 at L=2; check the claimed zonotope and tau_A(1).
import itertools, numpy as np
def syl(L):
    M=8; pats=set()
    for a in range(1,M):
        S=frozenset(x for x in range(M) if bin(a&x).count('1')%2==1)
        if 0 in S: S=frozenset(range(M))-S
        pats.add(tuple(sorted(S)))
    allk=np.array(list(itertools.product((-1,0,1),repeat=7)))
    ok=np.ones(len(allk),bool)
    for P in pats: ok&=np.abs(allk[:,[x-1 for x in P]].sum(1))<=L
    return allk[ok]
for L in (2,3):
    A=syl(L); print('syl8 L',L,'|A|',len(A),'tau(1)',len(set(A.sum(1)))-1)
A=syl(3); Aset=set(map(tuple,A))
e=lambda i: np.eye(7,dtype=int)[i-1]
Z=[-e(1),-e(1),-e(2),-e(2),-e(3),-e(4),-e(5),-e(6),e(3)-e(4),e(5)-e(6)]
sums=set(tuple(sum(s*z for s,z in zip(S,Z))) for S in itertools.product((0,1),repeat=10))
good=[tuple(k0) for k0 in A if all(tuple(k0+np.array(s)) in Aset for s in sums)]
print('M8\' zonotope translates:',len(good),good[:3])
