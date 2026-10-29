# Killer's own point count #S(F_{3^n}) for S: w^2 = f(x,y,z), n = 1..8, and Frobenius char poly on NS.
import sys, itertools, time
import numpy as np
import sympy as sp

p = 3
X, Y, Z = sp.symbols('x y z')
f = X**3*Y + 2*X**2*Y*Z - 2*X*Y**3 - 2*X*Y*Z**2 + 2*X*Z**3 + 2*Y**4 + 2*Y**3*Z + Y**2*Z**2 - Y*Z**3 + Z**4
monos = [(m, int(c) % p) for m, c in sp.Poly(f, X, Y, Z).terms() if int(c) % p]

def irreducible_poly(n):
    # find a primitive polynomial of degree n over F_3 (so that t is a generator)
    T = sp.symbols('T')
    for coeffs in itertools.product(range(p), repeat=n):
        P = sp.Poly([1] + list(coeffs), T, modulus=p)
        if coeffs[-1] == 0 or not P.is_irreducible:
            continue
        # check primitivity: order of T is p^n-1
        q = p**n; ok = True
        for r in sp.primefactors(q - 1):
            if sp.Poly(T, T, modulus=p).pow((q - 1)//r).rem(P) if False else sp.rem(sp.Poly(T, T, modulus=p)**((q-1)//r), P) == sp.Poly(1, T, modulus=p):
                ok = False; break
        if ok:
            return [1] + list(coeffs)
    raise RuntimeError

def build_field(n):
    q = p**n
    poly = irreducible_poly(n)  # monic, highest first
    # element encoding: integer sum d_i 3^i, d_i coefficient of T^i
    exp = np.zeros(q - 1, dtype=np.int64)
    cur = [1] + [0]*(n - 1)
    for k in range(q - 1):
        exp[k] = sum(cur[i]*p**i for i in range(n))
        # multiply by T
        top = cur[n - 1]
        cur = [0] + cur[:n - 1]
        # T^n = -(poly[1] T^{n-1} + ... + poly[n])
        for i in range(n):
            cur[i] = (cur[i] - top*poly[n - i]) % p
    log = np.full(q, -1, dtype=np.int64)
    log[exp] = np.arange(q - 1)
    assert (log[1:] >= 0).all(), "not primitive"
    digits = np.array([[(v // p**i) % p for i in range(n)] for v in range(q)], dtype=np.int64)
    pw = p**np.arange(n)
    return q, exp, log, digits, pw

def run(n):
    q, exp, log, digits, pw = build_field(n)
    def add(u, v):
        return (((digits[u] + digits[v]) % p) @ pw)
    def mul(u, v):
        u = np.asarray(u); v = np.asarray(v)
        out = exp[(log[u] + log[v]) % (q - 1)]
        return np.where((u == 0) | (v == 0), 0, out)
    def powe(u, k):
        u = np.asarray(u)
        if k == 0:
            return np.ones_like(u)
        out = exp[(log[u]*k) % (q - 1)]
        return np.where(u == 0, 0, out)
    def const(c):
        return c % p  # c in F_p encodes as digit 0
    chi = np.zeros(q, dtype=np.int64)
    chi[exp] = np.where(np.arange(q - 1) % 2 == 0, 1, -1)
    allx = np.arange(q)
    total = 0
    t0 = time.time()
    # affine chart z = 1: points (x, y, 1)
    for yv in range(q):
        yy = np.full(q, yv)
        val = np.zeros(q, dtype=np.int64)
        for (i, j, k), c in monos:
            term = mul(mul(powe(allx, i), powe(yy, j)), np.full(q, const(c)))
            val = add(val, term)
        total += int(chi[val].sum())
    # line z = 0: (x, 1, 0) for x in F_q and (1, 0, 0)
    val = np.zeros(q, dtype=np.int64); one = np.ones(q, dtype=np.int64)
    for (i, j, k), c in monos:
        if k == 0:
            val = add(val, mul(mul(powe(allx, i), powe(one, j)), np.full(q, const(c))))
    total += int(chi[val].sum())
    c100 = sum(c for (i, j, k), c in monos if (i, j, k) == (4, 0, 0)) % p
    total += int(chi[np.array([const(c100)])][0])
    npts = q*q + q + 1 + total
    tr = sp.Rational(npts - q*q - 1, q)
    print(f"n={n} q={q} #S={npts} trace={tr} ({time.time()-t0:.1f}s)", flush=True)
    return int(tr)

if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    traces = [run(n) for n in range(1, N + 1)]
    print("traces:", traces)
    # Newton identities -> characteristic polynomial (degree 8) from power sums p_1..p_8
    if N >= 8:
        e = [sp.Integer(1)]
        for k in range(1, 9):
            s = sum((-1)**(i - 1)*e[k - i]*traces[i - 1] for i in range(1, k + 1))
            e.append(sp.Rational(s, k))
        Xs = sp.symbols('X')
        cp = sum((-1)**k*e[k]*Xs**(8 - k) for k in range(9))
        print("char poly:", sp.factor(sp.expand(cp)))
