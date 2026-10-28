# Lattice points of m*P, P = anticanonical polytope of the permutohedral (Losev-Manin) variety Lbar_k
# M = {u in Z^k : sum u = 0}, P = {u : sum_{i in A} u_i >= -1 for all proper nonempty A}
import itertools, sys
from math import comb
def points(k,m):
    pts=[]
    rng=range(-m*(k-1), m*(k-1)+1)
    # enumerate u_1..u_{k-1}, u_k = -sum
    subsets=[A for r in range(1,k) for A in itertools.combinations(range(k),r)]
    def rec(pref):
        if len(pref)==k-1:
            u=pref+[-sum(pref)]
            if all(sum(u[i] for i in A)>=-m for A in subsets): pts.append(tuple(u))
            return
        for x in range(-m, m*(k-1)+1):
            p=pref+[x]
            # partial pruning: singletons >= -m already by range; pairs within prefix
            ok=True
            for r in range(1,len(p)+1):
                for A in itertools.combinations(range(len(p)),r):
                    if len(A)<k and sum(p[i] for i in A) < -m: ok=False;break
                    # upper: complement condition sum_{A} u <= m (since sum_{A^c} >= -m)
                    if len(A)<k and sum(p[i] for i in A) > m: ok=False;break
                if not ok: break
            if ok: rec(p)
    rec([])
    return pts
if __name__=='__main__':
    k=int(sys.argv[1])
    for m in range(1,int(sys.argv[2])+1):
        P=points(k,m); print(k,m,len(P),flush=True)
