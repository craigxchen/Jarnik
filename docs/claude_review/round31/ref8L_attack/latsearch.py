"""CRT/lattice constructions for the two-class problem.
For B points, choose bipartitions (S, m) with pairwise coprime moduli, form the lattice of
y in Z^{B-1} (y_0 = 0) occupying the prescribed two classes, LLL-reduce, and scan small
combinations for tuples with distinct coordinates.  Score = log q2(y) / log diam(y), with q2 the
exact two-class modulus of the found tuple (which may exceed the product of prescribed moduli).
"""
import sys, random, itertools
from math import log, prod
import sympy
from tc import q2, split_lattice, lll

def best_in_lattice(B, splits, k=None, coeff=(-1, 0, 1)):
    H, Q = split_lattice(B, splits)
    R = lll(H)
    R.sort(key=lambda v: max(abs(x) for x in v))
    k = min(len(R), k or 6)
    best = None
    for c in itertools.product(coeff, repeat=k):
        if not any(c): continue
        v = [sum(ci * R[i][j] for i, ci in enumerate(c)) for j in range(B - 1)]
        y = [0] + v
        if len(set(y)) < B: continue
        L = max(y) - min(y)
        if L < 2: continue
        q = q2(y)
        r = log(q) / log(L)
        if best is None or r > best[0]:
            best = (r, sorted(x - min(y) for x in y), q, Q)
    return best

def rand_primes(cnt, bits, used):
    out = []
    while len(out) < cnt:
        p = sympy.randprime(2 ** (bits - 1), 2 ** bits)
        if p not in used: used.add(p); out.append(p)
    return out

def all_splits(B, minside=1):
    res = []
    for m in range(1, 1 << (B - 1)):
        S = frozenset(i for i in range(1, B) if (m >> (i - 1)) & 1)
        if min(len(S) + 0, B - len(S)) >= minside and len(S) >= 1:
            res.append(S)
    return res

def cube_splits(n):
    B = 2 ** n
    # point index i <-> binary vector of i; split by coordinate j
    return B, [frozenset(i for i in range(1, B) if (i >> j) & 1) for j in range(n)]

def run(B, nsplits, bits, trials, seed, minside=2):
    random.seed(seed)
    pool = all_splits(B, minside)
    best = None
    for t in range(trials):
        Ss = random.sample(pool, min(nsplits, len(pool)))
        used = set()
        ms = rand_primes(len(Ss), bits, used)
        b = best_in_lattice(B, list(zip(Ss, ms)))
        if b and (best is None or b[0] > best[0]):
            best = b + (nsplits, bits)
    return best

if __name__ == "__main__":
    B = int(sys.argv[1]); bits = int(sys.argv[2]); trials = int(sys.argv[3])
    for ns in range(2, int(sys.argv[4]) + 1):
        b = run(B, ns, bits, trials, seed=1000 * B + ns)
        if b:
            print(f"B={B} nsplits={ns} bits={bits} best_ratio={b[0]:.4f} logq={log(b[2]):.1f} y={b[1]}", flush=True)
