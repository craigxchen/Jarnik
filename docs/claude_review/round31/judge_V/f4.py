# F_2[t]: python ints as bit vectors. F_4[t]: pairs (A,B) = A + B*theta, theta^2 = theta + 1.
def mul2(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        a <<= 1; b >>= 1
    return r
def deg2(a): return a.bit_length() - 1   # deg 0 poly = -1
def divmod2(a, b):
    q = 0; db = deg2(b)
    while a and deg2(a) >= db:
        s = deg2(a) - db; q ^= 1 << s; a ^= b << s
    return q, a
def gcd2(a, b):
    while b: a, b = b, divmod2(a, b)[1]
    return a
def der2(a):  # formal derivative over F_2
    r = 0; i = 1
    while (a >> i):
        if i & 1 and (a >> i) & 1: r |= 1 << (i - 1)
        i += 1
    return r
def mul4(x, y):
    A, B = x; C, D = y
    BD = mul2(B, D)
    return (mul2(A, C) ^ BD, mul2(A, D) ^ mul2(B, C) ^ BD)
def add4(x, y): return (x[0] ^ y[0], x[1] ^ y[1])
def conj4(x): return (x[0] ^ x[1], x[1])
def deg4(x): return max(deg2(x[0]), deg2(x[1]))
def norm4(x):
    p = mul4(x, conj4(x)); assert p[1] == 0; return p[0]
def der4(x): return (der2(x[0]), der2(x[1]))
ONE = (1, 0); THETA = (0, 1)
