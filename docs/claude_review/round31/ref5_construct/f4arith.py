"""Referee's own F_2[t] / F_4[t] arithmetic (independent of Astra's script).

F_2[t]: python int, bit i = coeff of t^i.
F_4 = F_2[theta], theta^2 = theta + 1.  F_4[t] element z = A + B*theta, stored as tuple (A,B).
Also a slow list-based F_4[t] (coeff ints 0..3, bit0 = 1, bit1 = theta) for cross-checks.
"""

def clmul(x, y):
    # schoolbook over the bits of the shorter operand, via 8-bit window table
    if x.bit_length() < y.bit_length():
        x, y = y, x
    tab = [0] * 256
    for i in range(1, 256):
        lo = i & -i
        tab[i] = tab[i ^ lo] ^ (x << (lo.bit_length() - 1))
    r = 0
    sh = 0
    while y:
        r ^= tab[y & 255] << sh
        y >>= 8
        sh += 8
    return r

def deg2(x):
    return x.bit_length() - 1

def divmod2(a, b):
    db = deg2(b)
    q = 0
    while a and deg2(a) >= db:
        s = deg2(a) - db
        q ^= 1 << s
        a ^= b << s
    return q, a

def gcd2(a, b):
    while b:
        a, b = b, divmod2(a, b)[1]
    return a

def mul4(x, y):
    A, B = x; C, D = y
    AC = clmul(A, C); BD = clmul(B, D)
    return (AC ^ BD, clmul(A, D) ^ clmul(B, C) ^ BD)

def add4(x, y):
    return (x[0] ^ y[0], x[1] ^ y[1])

def conj4(x):
    A, B = x
    return (A ^ B, B)          # A + B(theta+1)

def deg4(x):
    return (x[0] | x[1]).bit_length() - 1

def norm4(x):
    A, B = x
    return clmul(A, A) ^ clmul(A, B) ^ clmul(B, B)

# ---- list-based F_4[t] for gcd and cross checks ----
F4MUL = [[0]*4 for _ in range(4)]
def _f4m(a, b):
    # a = a0 + a1 th
    a0, a1 = a & 1, a >> 1; b0, b1 = b & 1, b >> 1
    c0 = (a0 & b0) ^ (a1 & b1)
    c1 = (a0 & b1) ^ (a1 & b0) ^ (a1 & b1)
    return c0 | (c1 << 1)
for a in range(4):
    for b in range(4):
        F4MUL[a][b] = _f4m(a, b)
F4INV = {1: 1, 2: 3, 3: 2}

def to_list(x):
    A, B = x
    n = max(A.bit_length(), B.bit_length())
    return [((A >> i) & 1) | (((B >> i) & 1) << 1) for i in range(n)]

def from_list(l):
    A = B = 0
    for i, c in enumerate(l):
        if c & 1: A |= 1 << i
        if c & 2: B |= 1 << i
    return (A, B)

def lmod(a, b):
    a = a[:]
    while a and a[-1] == 0: a.pop()
    db = len(b) - 1
    inv = F4INV[b[-1]]
    while len(a) - 1 >= db:
        c = F4MUL[a[-1]][inv]
        s = len(a) - 1 - db
        for i, bc in enumerate(b):
            a[s + i] ^= F4MUL[c][bc]
        while a and a[-1] == 0: a.pop()
    return a

def lgcd(a, b):
    a = a[:]; b = b[:]
    while a and a[-1] == 0: a.pop()
    while b and b[-1] == 0: b.pop()
    while b:
        a, b = b, lmod(a, b)
    if a:
        inv = F4INV[a[-1]]
        a = [F4MUL[c][inv] for c in a]
    return a

def lmul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] ^= F4MUL[x][y]
    return r

def ldiv_exact(a, b):
    # a / b, assert exact
    a = a[:]
    while a and a[-1] == 0: a.pop()
    db = len(b) - 1
    inv = F4INV[b[-1]]
    q = [0] * (len(a) - db)
    while len(a) - 1 >= db:
        c = F4MUL[a[-1]][inv]
        s = len(a) - 1 - db
        q[s] = c
        for i, bc in enumerate(b):
            a[s + i] ^= F4MUL[c][bc]
        while a and a[-1] == 0: a.pop()
    assert not a
    return q

# ---- factorisation over F_2 (squarefree decomposition + DDF + EDF) ----
import random

def deriv2(f):
    r = 0; i = 1
    while (f >> i):
        if (f >> i) & 1 and i & 1:
            r |= 1 << (i - 1)
        i += 1
    return r

def mulmod2(a, b, m):
    return divmod2(clmul(a, b), m)[1]

def sqrt2(f):
    # f is a square in F_2[t]: take even coefficients
    r = 0; i = 0
    while (f >> (2 * i)):
        if (f >> (2 * i)) & 1:
            r |= 1 << i
        i += 1
    return r

def sqfree(f):
    # returns list of (g, e), g squarefree, f = prod g^e
    res = []
    def rec(f, mult):
        if deg2(f) <= 0: return
        df = deriv2(f)
        if df == 0:
            rec(sqrt2(f), mult * 2); return
        c = gcd2(f, df)
        w = divmod2(f, c)[0]
        i = 1
        while deg2(w) > 0:
            y = gcd2(w, c)
            fac = divmod2(w, y)[0]
            if deg2(fac) > 0: res.append((fac, i * mult))
            w = y
            c = divmod2(c, y)[0]
            i += 1
        if deg2(c) > 0:
            rec(sqrt2(c), mult * 2)
    rec(f, 1)
    return res

def ddf(f):
    out = []
    h = 2  # t
    i = 0
    while deg2(f) >= 2 * (i + 1):
        i += 1
        h = mulmod2(h, h, f)
        g = gcd2(f, h ^ 2)
        if deg2(g) > 0:
            out.append((g, i))
            f = divmod2(f, g)[0]
            h = divmod2(h, f)[1]
    if deg2(f) > 0:
        out.append((f, deg2(f)))
    return out

def edf(f, d):
    n = deg2(f)
    if n == d: return [f]
    while True:
        a = random.getrandbits(n) | 1
        # trace map a + a^2 + ... + a^(2^(d-1)) mod f
        tr = a; x = a
        for _ in range(d - 1):
            x = mulmod2(x, x, f); tr ^= x
        g = gcd2(f, tr)
        if 0 < deg2(g) < n:
            return edf(g, d) + edf(divmod2(f, g)[0], d)

def factor2(f):
    fac = {}
    for g, e in sqfree(f):
        for h, d in ddf(g):
            for p in edf(h, d):
                fac[p] = fac.get(p, 0) + e
    return fac
