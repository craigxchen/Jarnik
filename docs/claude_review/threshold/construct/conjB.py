"""Test Conjecture B:  reg(S) <= MW(S) := 2 * E_{P iid uniform subset of [M]} w_P(S)
on many finite S in a slice of Z^M; also test the 'good pair' condition
   exists i<j :  E[w_P ; P separates i,j] >= E[w_P ; P does not separate i,j].
reg is computed mod two primes (rigorous upper bound; equality of the two primes reported)."""
import sys, random, itertools
from fractions import Fraction as Fr
from clib import reg2, proj, widths, popc, separates
from gens import Qslice, orbit, SUsym

def MW_and_pairs(S, M):
    w = widths(S, M)          # cuts P with bit M-1 clear; P and complement have same width
    # iid uniform P subset of [M]: each nontrivial cut (up to complement) has prob 2/2^M
    tot = sum(w.values())
    MW = Fr(2 * 2 * tot, 2**M)
    best = None
    for i in range(M):
        for j in range(i+1, M):
            sep = sum(v for P, v in w.items() if separates(P, i, j))
            same = tot - sep
            d = sep - same
            if best is None or d > best: best = d
    return MW, best

random.seed(12345)
M = int(sys.argv[1]); trials = int(sys.argv[2]); R = int(sys.argv[3])
worst = None; badpair = 0
for t in range(trials):
    kind = random.choice(["rand", "randconv", "orbitunion"])
    if kind == "rand":
        npts = random.randint(2, 60)
        S = set()
        while len(S) < npts:
            k = [random.randint(-R, R) for _ in range(M-1)]
            S.add(tuple(k) + (-sum(k),))
        S = sorted(S)
    elif kind == "randconv":
        # lattice points of a random cut polytope {sum=0, k(P) <= U_P}
        U = {}
        for r in range(1, M):
            for P in itertools.combinations(range(M), r):
                U[P] = random.randint(0, R)
        S = []
        for k in itertools.product(range(-R, R+1), repeat=M-1):
            kk = k + (-sum(k),)
            if all(sum(kk[j] for j in P) <= U[P] for P in U): S.append(kk)
        if len(S) < 2: continue
    else:
        S = set()
        for _ in range(random.randint(1, 3)):
            v = [random.randint(-R, R) for _ in range(M-1)]
            v.append(-sum(v))
            S |= set(orbit(tuple(v)))
        S = sorted(S)
        if len(S) < 2: continue
    if len(S) > 2500: continue
    r = reg2(proj(S))
    MW, d = MW_and_pairs(S, M)
    if d < 0: badpair += 1
    marg = MW - r
    if worst is None or marg < worst[0]:
        worst = (marg, kind, r, MW, len(S), S if len(S) < 40 else None)
    if marg < 0:
        print("COUNTEREXAMPLE", kind, r, MW, S, flush=True)
print(f"M={M} trials={trials}: min(MW - reg) = {worst[0]} (kind={worst[1]}, reg={worst[2]}, MW={worst[3]}, |S|={worst[4]}) ; sets with no good pair: {badpair}")
if worst[5] is not None: print("  worst set:", worst[5])
