"""Exact rational kernel vectors of an integer matrix via multi-modular nullspace (62-bit primes) + CRT + rational
reconstruction; each returned vector is VERIFIED exactly (E @ v == 0 over Z)."""
import math, numpy as np
from fractions import Fraction
from flint import nmod_mat
from sympy import prevprime

def primes62(k):
    ps = []; p = 2 ** 62
    while len(ps) < k:
        p = prevprime(p); ps.append(p)
    return ps

def ratrec_big(a, m):
    a %= m; bound = math.isqrt(m // 2)
    r0, r1 = m, a; s0, s1 = 0, 1
    while r1 > bound:
        q = r0 // r1; r0, r1 = r1, r0 - q * r1; s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound: return None
    return Fraction(r1, s1)

def verify(E, iv):
    """exact check E @ iv == 0 over Z (int64 fast path only when no overflow is possible)"""
    mx = max(abs(z) for z in iv)
    if mx * int(np.abs(E).max()) * E.shape[1] < 2 ** 62:
        return not (E @ np.array(iv, dtype=np.int64)).any()
    return not any(E.astype(object).dot(np.array(iv, dtype=object)))

def kernel_exact(E, want=None, nprimes=4, select=None):
    """E: int64 2D array.  Returns list of exact integer kernel vectors (as Python int lists).
    select(vec_mod_p0) -> bool decides (cheaply, mod first prime) whether to lift a basis vector."""
    ps = primes62(nprimes)
    Xs = []; piv_sig = None
    for p in ps:
        M = nmod_mat(E.shape[0], E.shape[1], [int(x) % p for x in E.ravel()], p)
        X, nul = M.nullspace()
        Xs.append((X, nul))
    nul = Xs[0][1]
    assert all(x[1] == nul for x in Xs), "nullity differs between primes"
    mod = 1
    for p in ps: mod *= p
    out = []
    Eo = None
    for c in range(nul):
        v0 = [int(Xs[0][0][k, c]) for k in range(E.shape[1])]
        if select is not None and not select(v0, ps[0]): continue
        # CRT
        vec = []
        ok = True
        for k in range(E.shape[1]):
            r = 0; mm = 1
            for (X, _), p in zip(Xs, ps):
                a = int(X[k, c])
                t = ((a - r) * pow(mm, -1, p)) % p
                r += mm * t; mm *= p
            f = ratrec_big(r, mod)
            if f is None: ok = False; break
            vec.append(f)
        if not ok: continue
        den = 1
        for f in vec: den = den * f.denominator // math.gcd(den, f.denominator)
        iv = [int(f * den) for f in vec]
        if not verify(E, iv): continue
        out.append(iv)
        if want is not None and len(out) >= want: break
    return out
