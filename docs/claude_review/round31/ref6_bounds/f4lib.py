# Own arithmetic: F_2[t] polynomials as Python ints (bit k = coeff of t^k);
# F_4[t] polynomials as pairs (A,B) meaning A + theta*B, theta^2 = theta + 1.
# Conjugation theta -> theta+1:  A + (theta+1)B = (A^B) + theta B.

def deg(a):
    return a.bit_length() - 1  # deg 0 -> -1

def clmul(a, b):
    if a.bit_length() < b.bit_length():
        a, b = b, a
    r = 0
    while b:
        low = b & -b
        r ^= a << (low.bit_length() - 1)
        b ^= low
    return r

def divmod2(a, b):
    if b == 0:
        raise ZeroDivisionError
    db = deg(b); q = 0
    while a and deg(a) >= db:
        s = deg(a) - db
        q ^= 1 << s
        a ^= b << s
    return q, a

def mod2(a, b):
    db = deg(b)
    while a and deg(a) >= db:
        a ^= b << (deg(a) - db)
    return a

def gcd2(a, b):
    while b:
        a, b = b, mod2(a, b)
    return a

def powmod2(a, e, p):
    r = 1; a = mod2(a, p)
    while e:
        if e & 1:
            r = mod2(clmul(r, a), p)
        a = mod2(clmul(a, a), p)
        e >>= 1
    return r

def is_irred2(p):
    # Rabin/Ben-Or style: p irreducible iff gcd(t^(2^i)-t, p)=1 for i<=n/2
    n = deg(p)
    if n <= 0:
        return False
    x = 2; tpow = 2
    for i in range(1, n // 2 + 1):
        tpow = mod2(clmul(tpow, tpow), p)
        if gcd2(tpow ^ x, p) != 1:
            return False
    return True

def factor2(a):
    """trial-division factorisation of small a in F_2[t]; returns dict prime->exp"""
    res = {}
    p = 2
    # iterate over candidate polys in increasing order
    while a != 1 and deg(a) > 0:
        if 2 * deg(p) > deg(a):
            res[a] = res.get(a, 0) + 1
            break
        if is_irred2(p):
            while True:
                q, r = divmod2(a, p)
                if r:
                    break
                a = q
                res[p] = res.get(p, 0) + 1
        p += 1
    return res

# ---------- F_4[t] ----------
def d4(x):
    return max(deg(x[0]), deg(x[1]))

def add4(x, y):
    return (x[0] ^ y[0], x[1] ^ y[1])

def mul4(x, y):
    A, B = x; C, D = y
    AC = clmul(A, C); BD = clmul(B, D)
    ADBC = clmul(A ^ B, C ^ D) ^ AC ^ BD  # AD+BC
    return (AC ^ BD, ADBC ^ BD)

def conj4(x):
    return (x[0] ^ x[1], x[1])

def scal(x, c):
    # multiply by c in {1,2,3} meaning 1, theta, theta^2
    A, B = x
    if c == 1:
        return x
    if c == 2:
        return (B, A ^ B)
    if c == 3:
        return (A ^ B, A)
    raise ValueError

def lc4(x):
    n = d4(x)
    if n < 0:
        return 0
    return ((x[0] >> n) & 1) | (((x[1] >> n) & 1) << 1)  # 1 -> 1, 2 -> theta, 3 -> 1+theta=theta^2

INV = {1: 1, 2: 3, 3: 2}
MULT = {(1, 1): 1, (1, 2): 2, (1, 3): 3, (2, 2): 3, (2, 3): 1, (3, 3): 2}
def fmul(a, b):
    if a == 0 or b == 0:
        return 0
    return MULT[(min(a, b), max(a, b))]

def shl(x, s):
    return (x[0] << s, x[1] << s)

def mod4(a, b):
    db = d4(b); ib = INV[lc4(b)]
    while True:
        da = d4(a)
        if da < db:
            return a
        c = fmul(lc4(a), ib)
        a = add4(a, shl(scal(b, c), da - db))

def divmod4(a, b):
    db = d4(b); ib = INV[lc4(b)]
    q = (0, 0)
    while True:
        da = d4(a)
        if da < db:
            return q, a
        c = fmul(lc4(a), ib)
        q = add4(q, shl(scal((1, 0), c), da - db))
        a = add4(a, shl(scal(b, c), da - db))

def monic4(a):
    c = lc4(a)
    return scal(a, INV[c]) if c else a

def gcd4(a, b):
    while d4(b) >= 0:
        a, b = b, mod4(a, b)
    return monic4(a) if d4(a) >= 0 else a

def prod4(lst):
    r = (1, 0)
    for x in lst:
        r = mul4(r, x)
    return r

def pow4(x, e):
    r = (1, 0)
    while e:
        if e & 1:
            r = mul4(r, x)
        x = mul4(x, x)
        e >>= 1
    return r

def is_real(x):
    return x[1] == 0

def depth_ratio(f, g):
    """ord_inf(f/g - 1) for f,g of the same degree, f != g"""
    return d4(g) - d4(add4(f, g))
