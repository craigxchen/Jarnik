"""Independent F_4[t] / F_2[t] arithmetic (referee ref5_lower).

F_4 element: int 0..3, bit0 = constant part, bit1 = theta part, theta^2 = theta + 1.
F_4[t] polynomial: tuple of F_4 ints, index = power of t, no trailing zeros (zero poly = ()).
F_2[t] elements of V are ints (bit i = coeff of t^i); they embed with theta-part 0.
Conjugation: theta -> theta + 1, i.e. a + b theta -> (a+b) + b theta.
"""

# multiplication table built from the rule, not typed by hand
def _mul_f4(x, y):
    a, b = x & 1, x >> 1
    c, d = y & 1, y >> 1
    # (a + b th)(c + d th) = ac + (ad+bc) th + bd (th+1)
    r0 = (a & c) ^ (b & d)
    r1 = (a & d) ^ (b & c) ^ (b & d)
    return r0 | (r1 << 1)

MUL = [[_mul_f4(x, y) for y in range(4)] for x in range(4)]
INV = [None] + [next(y for y in range(1, 4) if MUL[x][y] == 1) for x in range(1, 4)]
CONJ = [((x & 1) ^ (x >> 1)) | (x & 2) for x in range(4)]
assert MUL[2][2] == 3  # theta^2 = theta + 1
assert all(MUL[x][CONJ[x]] in (0, 1) for x in range(4))  # norms in F_2


def trim(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return tuple(p)


def deg(p):
    return len(p) - 1 if p else -10**9


def add(p, q):
    n = max(len(p), len(q))
    r = [0] * n
    for i, x in enumerate(p):
        r[i] = x
    for i, x in enumerate(q):
        r[i] ^= x
    return trim(r)


def mul(p, q):
    if not p or not q:
        return ()
    r = [0] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        if x:
            row = MUL[x]
            for j, y in enumerate(q):
                if y:
                    r[i + j] ^= row[y]
    return trim(r)


def conj(p):
    return tuple(CONJ[x] for x in p)


def divmod_(p, q):
    p = list(p)
    dq = deg(q)
    il = INV[q[-1]]
    quo = [0] * max(0, len(p) - dq)
    while len(p) - 1 >= dq and p:
        sh = len(p) - 1 - dq
        c = MUL[p[-1]][il]
        quo[sh] = c
        row = MUL[c]
        for j, y in enumerate(q):
            if y:
                p[sh + j] ^= row[y]
        while p and p[-1] == 0:
            p.pop()
    return trim(quo), tuple(p)


def gcd(p, q):
    while q:
        p, q = q, divmod_(p, q)[1]
    if not p:
        return p
    c = INV[p[-1]]
    return trim([MUL[c][x] for x in p])


def from_f2(v):
    """int bitmask in F_2[t] -> F_4[t] poly"""
    r = []
    while v:
        r.append(v & 1)
        v >>= 1
    return tuple(r)


def f2deg(v):
    return v.bit_length() - 1  # deg 0 = -1 (only used for nonzero)


def linfac(T, v):
    """Y + v with Y = T + theta, T, v in F_2[t] as ints"""
    p = list(from_f2(T ^ v)) or [0]
    p[0] ^= 2
    return trim(p)


def prod(polys):
    r = (1,)
    for p in polys:
        r = mul(r, p)
    return r


def f2mul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        a <<= 1
        b >>= 1
    return r


def parity(x):
    return bin(x).count("1") & 1


def forms_and_sides(d):
    """forms l with l(1)=1: l(v) = parity(v & w), w odd. sides: (w, side bit) -> frozenset."""
    Q = 1 << d
    out = []
    for w in range(1, Q, 2):
        ker = frozenset(v for v in range(Q) if parity(v & w) == 0)
        out.append((w, ker, frozenset(range(Q)) - ker))
    return out
