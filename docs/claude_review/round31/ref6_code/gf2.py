"""Own F_2[t] arithmetic (polynomials as Python ints, bit i = coeff of t^i) and factorisation,
plus F_4[t] arithmetic (F_4 = {0,1,2,3}, 2 = theta, 3 = theta+1 = theta^2, theta^2+theta+1=0)."""
import random

def deg(a):
    return a.bit_length() - 1

def mul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        a <<= 1
        b >>= 1
    return r

def divmod2(a, b):
    assert b
    q = 0
    db = deg(b)
    while a and deg(a) >= db:
        s = deg(a) - db
        q ^= 1 << s
        a ^= b << s
    return q, a

def mod(a, b):
    return divmod2(a, b)[1]

def gcd(a, b):
    while b:
        a, b = b, mod(a, b)
    return a

def mulmod(a, b, m):
    return mod(mul(a, b), m)

def powmod(a, e, m):
    r = 1
    a = mod(a, m)
    while e:
        if e & 1:
            r = mulmod(r, a, m)
        a = mulmod(a, a, m)
        e >>= 1
    return r

def deriv(a):
    # d/dt: keep odd-degree coefficients, shifted down
    r = 0
    i = 1
    a >>= 1
    while a:
        if a & 1 and (i % 2 == 1):
            r |= 1 << (i - 1)
        a >>= 1
        i += 1
    return r

def sqrt_poly(a):
    # a must be a square: a = sum c_{2i} t^{2i}
    r = 0
    i = 0
    while a:
        assert not (a & 2) or True
        if a & 1:
            r |= 1 << i
        a >>= 2
        i += 1
    return r

def is_square(a):
    return deriv(a) == 0

def squarefree_decomp(f):
    """Return list of (g, multiplicity) with g squarefree (not nec. irreducible)."""
    out = []
    def rec(f, mult):
        if deg(f) <= 0:
            return
        d = deriv(f)
        if d == 0:
            rec(sqrt_poly(f), 2 * mult)
            return
        g = gcd(f, d)
        w = divmod2(f, g)[0]
        i = 1
        while deg(w) > 0:
            y = gcd(w, g)
            z = divmod2(w, y)[0]
            if deg(z) > 0:
                out.append((z, i * mult))
            i += 1
            w = y
            g = divmod2(g, y)[0]
        if deg(g) > 0:
            # g is a square
            rec(sqrt_poly(g), 2 * mult)
    rec(f, 1)
    return out

def ddf(f):
    """distinct degree factorisation of squarefree f: list of (product, degree)."""
    res = []
    h = 2  # t
    k = 0
    while deg(f) >= 2 * (k + 1):
        k += 1
        h = mulmod(h, h, f)
        g = gcd(f, h ^ 2)
        if deg(g) > 0:
            res.append((g, k))
            f = divmod2(f, g)[0]
            h = mod(h, f)
    if deg(f) > 0:
        res.append((f, deg(f)))
    return res

def edf(f, k, rng):
    if deg(f) == k:
        return [f]
    while True:
        a = rng.getrandbits(deg(f)) | 0
        a = mod(a, f)
        if deg(a) <= 0:
            continue
        # trace map a + a^2 + ... + a^(2^(k-1)) mod f
        tr = 0
        x = a
        for _ in range(k):
            tr ^= x
            x = mulmod(x, x, f)
        g = gcd(f, tr)
        if 0 < deg(g) < deg(f):
            return edf(g, k, rng) + edf(divmod2(f, g)[0], k, rng)

def factor(f, seed=1):
    """Return dict irreducible -> multiplicity for monic f in F_2[t]."""
    rng = random.Random(seed)
    res = {}
    for g, mlt in squarefree_decomp(f):
        for prod, k in ddf(g):
            for p in edf(prod, k, rng):
                res[p] = res.get(p, 0) + mlt
    # verify
    chk = 1
    for p, e in res.items():
        for _ in range(e):
            chk = mul(chk, p)
    assert chk == f, "factorisation check failed"
    return res

def is_irreducible(f):
    fac = factor(f)
    return len(fac) == 1 and list(fac.values())[0] == 1

# ---------------- F_4 and F_4[t] ----------------
F4_MUL = [[0, 0, 0, 0], [0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2]]
F4_INV = [None, 1, 3, 2]
F4_CONJ = [0, 1, 3, 2]  # Frobenius x -> x^2: theta -> theta^2 = theta+1

def p4_trim(a):
    while a and a[-1] == 0:
        a.pop()
    return a

def p4_add(a, b):
    n = max(len(a), len(b))
    r = [(a[i] if i < len(a) else 0) ^ (b[i] if i < len(b) else 0) for i in range(n)]
    return p4_trim(r)

def p4_mul(a, b):
    if not a or not b:
        return []
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            row = F4_MUL[x]
            for j, y in enumerate(b):
                if y:
                    r[i + j] ^= row[y]
    return p4_trim(r)

def p4_divmod(a, b):
    a = list(a)
    b = p4_trim(list(b))
    assert b
    q = [0] * max(0, len(a) - len(b) + 1)
    inv = F4_INV[b[-1]]
    while len(a) >= len(b) and a:
        c = F4_MUL[a[-1]][inv]
        s = len(a) - len(b)
        q[s] = c
        for j, y in enumerate(b):
            a[s + j] ^= F4_MUL[c][y]
        p4_trim(a)
    return p4_trim(q), a

def p4_monic(a):
    inv = F4_INV[a[-1]]
    return [F4_MUL[inv][x] for x in a]

def p4_gcd(a, b):
    a = p4_trim(list(a)); b = p4_trim(list(b))
    while b:
        a, b = b, p4_divmod(a, b)[1]
    return p4_monic(a) if a else a

def p4_conj(a):
    return [F4_CONJ[x] for x in a]

def p4_from_gf2(f):
    r = []
    i = 0
    while f:
        r.append(f & 1)
        f >>= 1
    return r

def p4_to_gf2(a):
    r = 0
    for i, x in enumerate(a):
        assert x in (0, 1), "not in F_2[t]"
        if x:
            r |= 1 << i
    return r

def p4_ord(f, pi):
    """valuation of f at pi (f nonzero)."""
    k = 0
    while True:
        q, r = p4_divmod(f, pi)
        if r:
            return k
        f = q
        k += 1
