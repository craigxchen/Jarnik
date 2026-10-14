"""Fast exact J(P) generators: nullspace mod a 62-bit prime, rational reconstruction, EXACT verification over Z.
Only verified vectors are used (soundness: every returned form is an exact element of J(P)_i over Q)."""
import itertools, math
import numpy as np
from fractions import Fraction
from flint import nmod_mat, fmpz_mat
from sympy import prevprime
from loclib import bmonos

P62 = prevprime(2 ** 62)

def ratrec(a, m):
    """rational reconstruction of a mod m with |num|,|den| <= sqrt(m/2)"""
    a %= m
    bound = math.isqrt(m // 2)
    r0, r1 = m, a; s0, s1 = 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    return Fraction(r1, s1)

def jproj_gens_fast(P, b, maxdeg, verbose=False):
    P = [tuple(int(z) for z in p) for p in P]
    out = {}; rho = None
    prev_top = []   # exact integer basis (spanning) of J(P)_{i-1}
    fails = 0
    for i in range(1, maxdeg + 1):
        top = bmonos(b, i); low = [e for t in range(i) for e in bmonos(b, t)]
        mons = low + top
        E = [[math.prod(p[z] ** e[z] for z in range(b)) for e in mons] for p in P]
        Em = nmod_mat(len(P), len(mons), [x % P62 for r in E for x in r], P62)
        X, nul = Em.nullspace()
        # y * J_{i-1}
        topidx = {e: a for a, e in enumerate(top)}
        prevm = bmonos(b, i - 1)
        shifted = []
        for v in prev_top:
            for q in range(b):
                w = [0] * len(top)
                for a, e in enumerate(prevm):
                    if v[a]:
                        e2 = list(e); e2[q] += 1; w[topidx[tuple(e2)]] += v[a]
                shifted.append(w)
        cur = [list(r) for r in shifted]
        r = fmpz_mat(cur).rank() if cur else 0
        new = []
        Ez = None
        for c in range(nul):
            vec = [int(X[k, c]) for k in range(len(mons))]
            if not any(vec[len(low):]):
                continue
            fr = [ratrec(x, P62) for x in vec]
            if any(f is None for f in fr):
                fails += 1; continue
            den = 1
            for f in fr: den = den * f.denominator // math.gcd(den, f.denominator)
            iv = [int(f * den) for f in fr]
            # exact verification
            if Ez is None: Ez = np.array(E, dtype=object)
            if any(int(sum(Ez[k, j] * iv[j] for j in range(len(mons)) if iv[j])) != 0 for k in range(len(P))):
                fails += 1; continue
            tv = iv[len(low):]
            r2 = fmpz_mat(cur + [tv]).rank()
            if r2 > r:
                cur.append(tv); new.append(tv); r = r2
        out[i] = [{e: v[a] for a, e in enumerate(top) if v[a]} for v in new]
        prev_top = cur
        if r == len(top) and rho is None:
            rho = i - 1
        if verbose:
            print(f"      J(P) deg {i}: dim {r}/{len(top)}, new {len(new)}, ratrec/verify fails so far {fails}", flush=True)
    return out, rho, fails
