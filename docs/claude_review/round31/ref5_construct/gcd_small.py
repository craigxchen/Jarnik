"""gcd of all z_H via valuations at small irreducibles of F_4[t].
A common irreducible factor f of all z_H must divide two distinct factors Y+v, Y+v' (since the intersection of
all H is {0}), hence f | v+v' in V, so deg f <= d-1.  We enumerate all monic irreducibles of F_4[t] of degree <= d-1.
usage: python3 gcd_small.py r
"""
import sys, itertools
sys.path.insert(0, __import__('os').path.dirname(__file__))
from f4arith import *

def par(x): return bin(x).count('1') & 1

def irreducibles(maxdeg):
    irr = []
    for n in range(1, maxdeg + 1):
        for coeffs in itertools.product(range(4), repeat=n):
            f = list(coeffs) + [1]
            ok = True
            for g in irr:
                if 2 * (len(g) - 1) > n: break
                if not lmod(f, g):
                    ok = False; break
            if ok: irr.append(f)
    return irr

def lpow(f, j):
    r = [1]
    for _ in range(j): r = lmul(r, f)
    return r

def run(r):
    m = 2 ** r; d = r + 1; Q = 2 ** d
    L = m * (r - 1); D = (m - 1) * L
    forms = [w for w in range(1, Q) if w & 1]
    H = {w: [v for v in range(Q) if par(v & w) == 0 and v] for w in forms}
    irr = irreducibles(d - 1)
    Vl = [[(w >> i) & 1 for i in range(w.bit_length())] for w in range(1, Q)]
    irr = [f for f in irr if any(not lmod(w, f) for w in Vl)]
    Y = [0] * (L + 1); Y[L] = 1; Y[0] = 2           # t^L + theta  (theta encoded as 2)
    total = 0; detail = []
    for f in irr:
        val = {}
        for v in range(1, Q):
            Yv = Y[:]
            for i in range(v.bit_length()):
                if (v >> i) & 1: Yv[i] ^= 1
            j = 0
            while True:
                if lmod(Yv, lpow(f, j + 1)): break
                j += 1
            val[v] = j
        e = min(sum(val[v] for v in H[w]) for w in forms)
        if e:
            total += e * (len(f) - 1)
            detail.append((f, e, sorted(set(val.values()))))
    print(f"m={m}: deg gcd = {total}; factors (f coeffs low->high, exponent):",
          [(f, e) for f, e, _ in detail])
    return total

if __name__ == '__main__':
    run(int(sys.argv[1]))
