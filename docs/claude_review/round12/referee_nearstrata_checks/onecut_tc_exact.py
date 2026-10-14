# Exact check (small cases) of Lemma 2.5 / Thm 2.5: strata of types I, II, III lie in TC_V(A_{S,L}, m),
# and random rational points off Sigma are not in TC_V; also kappa lower bound via the explicit gadget is not
# re-used here (independent): we only test membership.
import itertools, random
from fractions import Fraction as Fr
import flint
random.seed(1)
def setup(d, S, L, m):
    A = [k for k in itertools.product((-1,0,1), repeat=d) if abs(sum(k[x] for x in S)) <= L]
    mons = [e for e in itertools.product((0,1,2), repeat=d) if sum(e) < m]
    M = flint.fmpz_mat(len(mons), len(A))
    for i,e in enumerate(mons):
        for jx,k in enumerate(A):
            v = 1
            for x in range(d):
                if e[x]: v *= k[x]**e[x]
            M[i,jx] = v
    X, nul = M.nullspace()
    basis = [[int(X[r,c]) for r in range(len(A))] for c in range(nul)]
    return A, basis
def Phi(A, c, u, m):
    return sum(ci * sum(k[x]*u[x] for x in range(len(u)))**m for ci,k in zip(c,A) if ci != 0)
def in_TC(A, basis, u, m):
    return all(Phi(A, c, u, m) == 0 for c in basis)
def rq(): return Fr(random.randint(-30,30), random.randint(1,7))
def check(d, S, L, m, trials=6):
    A, basis = setup(d, S, L, m)
    R = [x for x in range(d) if x not in S]; s = len(S)
    out = {'dimV': len(basis), 'U': len(A)}
    res = {'I':[], 'II':[], 'III':[], 'off':[]}
    for X in itertools.chain.from_iterable(itertools.combinations(range(d), r) for r in range(d+1)):
        if 2*len(X) < m:
            for _ in range(2):
                u = [rq() if x in X else Fr(0) for x in range(d)]; res['I'].append(in_TC(A,basis,u,m))
        if 2*L + 2*len(X) < m:
            for _ in range(2):
                a = rq(); u = [(rq() if x in X else Fr(0)) + (a if x in S else 0) for x in range(d)]
                res['II'].append(in_TC(A,basis,u,m))
    for Y in itertools.chain.from_iterable(itertools.combinations(R, r) for r in range(len(R)+1)):
        if s + L + 2*len(Y) < m:
            for _ in range(2):
                u = [rq() if (x in S or x in Y) else Fr(0) for x in range(d)]; res['III'].append(in_TC(A,basis,u,m))
    for _ in range(trials):
        u = [rq() for x in range(d)]; res['off'].append(in_TC(A,basis,u,m))
    return out, {k:(sum(v), len(v)) for k,v in res.items()}
for (d,S,L,m) in [(5,(0,1,2),1,4),(5,(0,1,2),1,5),(6,(0,1,2,3),1,5),(6,(0,1,2,3),1,6),(6,(0,1,2,3),2,6),(6,(0,1,2),1,6),(7,(0,1,2,3),1,6)]:
    print((d,S,L,m), check(d,S,L,m), flush=True)
