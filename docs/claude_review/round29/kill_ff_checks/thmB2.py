# Theorem B: check distinctness of columns AFTER the sign normalisation, for (a) B = basis-first choice, (b) B containing a triple
import itertools, numpy as np
rng=np.random.default_rng(1)
def test(k,c,Bchoice):
    D=2**k; L=D//2-c-1
    us=list(itertools.product([0,1],repeat=k)); nz=[u for u in us if any(u)]
    basis=[tuple(1 if i==j else 0 for i in range(k)) for j in range(k)]
    if Bchoice=='basis': B=basis+[u for u in nz if u not in basis and sum(u)>=3][:L-k] + [u for u in nz if u not in basis and sum(u)==2][:max(0,L-k-len([u for u in nz if sum(u)>=3]))]
    else:
        e12=tuple((basis[0][i]+basis[1][i])%2 for i in range(k))
        B=basis+[e12]+[u for u in nz if u not in basis and u!=e12][:L-k-1]
    B=B[:L]; assert len(B)==L and len(set(B))==L
    W=np.array([[(-1)**(np.dot(b,u)%2) for u in us] for b in B],dtype=float)
    # odd-hyperplane obstruction: x with b.x=1 for all b
    obstr=[x for x in nz if all(np.dot(b,x)%2==1 for b in B)]
    g=rng.standard_normal(L)@W
    sgn=np.sign(g); Wn=W*sgn
    cols={tuple(np.round(Wn[:,j],9)) for j in range(D)}
    return D,L,c,Bchoice,len(obstr),len(cols)==D
for k,c in [(3,0),(4,0),(4,1),(4,2),(4,3),(5,0),(5,5),(5,10)]:
    for ch in ['basis','triple']:
        if ch=='triple' and 2**(k-1)-c-1<k+1: print(k,c,'triple impossible (L<k+1)'); continue
        print(test(k,c,ch))
