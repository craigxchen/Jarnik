"""Own F_2[t] and F_4[t] arithmetic (referee ref5_reconcile).
F_2[t]: Python int, bit k = coefficient of t^k.
F_4[t]: pair (A,B) meaning A + B*theta, theta^2 = theta + 1, A,B in F_2[t].
Conjugation: theta -> theta+1, i.e. (A,B) -> (A^B, B).
"""

def deg2(a):
    return a.bit_length() - 1  # deg 0 -> -1

def mul2(a, b):
    if a == 0 or b == 0:
        return 0
    if a.bit_length() < b.bit_length():
        a, b = b, a
    # 4-bit window on b
    tab = [0] * 16
    for k in range(1, 16):
        r = 0
        for j in range(4):
            if (k >> j) & 1:
                r ^= a << j
        tab[k] = r
    res = 0
    sh = 0
    while b:
        nib = b & 15
        if nib:
            res ^= tab[nib] << sh
        b >>= 4
        sh += 4
    return res

def divmod2(a, b):
    assert b
    db = deg2(b)
    q = 0
    while a and deg2(a) >= db:
        s = deg2(a) - db
        q ^= 1 << s
        a ^= b << s
    return q, a

def gcd2(a, b):
    while b:
        # a mod b
        db = deg2(b)
        while a and deg2(a) >= db:
            a ^= b << (deg2(a) - db)
        a, b = b, a
    return a

def deriv2(a):
    # formal derivative over F_2: keep odd-degree terms shifted down
    odd = a & int('10' * ((a.bit_length() + 1) // 2 + 1), 2)  # bits at odd positions
    return odd >> 1

# ---------- F_4 ----------
# scalars: 0,1,2=theta,3=theta+1 encoded as (bitA, bitB)
def smul4(c, P):
    A, B = P
    if c == 0:
        return (0, 0)
    if c == 1:
        return P
    if c == 2:  # theta*(A+B theta) = A theta + B(theta+1) = B + (A+B) theta
        return (B, A ^ B)
    # theta+1 = theta^2: theta*(theta*P)
    A1, B1 = B, A ^ B
    return (B1, A1 ^ B1)

INV4 = {1: 1, 2: 3, 3: 2}

def deg4(P):
    return max(deg2(P[0]), deg2(P[1]))

def lc4(P):
    d = deg4(P)
    return ((P[0] >> d) & 1) | (((P[1] >> d) & 1) << 1)

def add4(P, Q):
    return (P[0] ^ Q[0], P[1] ^ Q[1])

def mul4(P, Q):
    A, B = P
    C, D = Q
    AC = mul2(A, C)
    BD = mul2(B, D)
    X = mul2(A ^ B, C ^ D)
    return (AC ^ BD, X ^ AC)

def conj4(P):
    return (P[0] ^ P[1], P[1])

def norm4(P):
    A, B = P
    return mul2(A, A) ^ mul2(A, B) ^ mul2(B, B)

def shift4(P, s):
    return (P[0] << s, P[1] << s)

def iszero4(P):
    return P[0] == 0 and P[1] == 0

def divmod4(P, Q):
    dq = deg4(Q)
    inv = INV4[lc4(Q)]
    quo = (0, 0)
    R = P
    while not iszero4(R) and deg4(R) >= dq:
        s = deg4(R) - dq
        c = smul4(inv, (lc4(R) & 1, lc4(R) >> 1)) if False else None
        # coefficient: lc(R)/lc(Q)
        cl = lc4(R)
        coef = mulscal(cl, inv)
        T = shift4(smul4(coef, Q), s)
        R = add4(R, T)
        quo = add4(quo, ((coef & 1) << s, ((coef >> 1) & 1) << s))
    return quo, R

def mulscal(a, b):
    # F_4 scalar product, encoding 0,1,2(theta),3(theta+1)
    if a == 0 or b == 0:
        return 0
    if a == 1:
        return b
    if b == 1:
        return a
    # theta^i representation: 2 -> theta^1, 3 -> theta^2
    e = {2: 1, 3: 2}
    k = (e[a] + e[b]) % 3
    return {0: 1, 1: 2, 2: 3}[k]

def monic4(P):
    if iszero4(P):
        return P
    return smul4(INV4[lc4(P)], P)

def gcd4(P, Q):
    while not iszero4(Q):
        _, R = divmod4(P, Q)
        P, Q = Q, R
    return monic4(P)

def from2(a):
    return (a, 0)

def poly4_from_coeffs(coeffs):
    """coeffs: list of scalars 0..3, index = degree"""
    A = B = 0
    for k, c in enumerate(coeffs):
        if c & 1:
            A |= 1 << k
        if c & 2:
            B |= 1 << k
    return (A, B)
