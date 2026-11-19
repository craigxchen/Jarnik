# Own finite-field and polynomial arithmetic in characteristic 2.
import itertools, random

class GF2k:
    """F_{2^k} with polynomial basis, elements are ints < 2^k."""
    POLYS = {4: 0b10011, 5: 0b100101, 6: 0b1000011, 8: 0b100011011, 16: (1 << 16) | 0b101101,
             20: (1 << 20) | 0b1001}
    def __init__(self, k):
        self.k = k; self.mod = self.POLYS[k]; self.q = 1 << k
        # check irreducibility crudely by x^(2^k) == x and no smaller subfield for proper divisors
        assert self.frob_pow(2, k) == 2
        for d in range(1, k):
            if k % d == 0:
                assert self.frob_pow(2, d) != 2, "modulus not irreducible"
    def mul(self, a, b):
        r = 0
        while b:
            if b & 1: r ^= a
            b >>= 1; a <<= 1
            if a >> self.k: a ^= self.mod
        return r
    def pow(self, a, e):
        r = 1
        while e:
            if e & 1: r = self.mul(r, a)
            a = self.mul(a, a); e >>= 1
        return r
    def frob_pow(self, a, j):
        for _ in range(j): a = self.mul(a, a)
        return a
    def inv(self, a):
        assert a
        return self.pow(a, self.q - 2)
    def rand(self): return random.randrange(self.q)

def rank_gf(F, M):
    """rank of matrix (list of rows) over field F (GF2k)."""
    M = [row[:] for row in M]; r = 0
    ncol = len(M[0]) if M else 0
    for c in range(ncol):
        piv = None
        for i in range(r, len(M)):
            if M[i][c]: piv = i; break
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        iv = F.inv(M[r][c])
        M[r] = [F.mul(iv, x) for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [x ^ F.mul(f, y) for x, y in zip(M[i], M[r])]
        r += 1
    return r

# ---- multivariate polynomials over F_2: frozenset-like dict {exp tuple: 1} represented as set
def mp_add(a, b): return a ^ b
def mp_mul(a, b):
    r = set()
    for e in a:
        for f in b:
            m = tuple(x + y for x, y in zip(e, f))
            if m in r: r.remove(m)
            else: r.add(m)
    return r
def mp_var(i, n):
    e = [0] * n; e[i] = 1; return {tuple(e)}
def mp_const1(n): return {tuple([0] * n)}
def mp_diff(a, i):
    r = set()
    for e in a:
        if e[i] % 2 == 1:
            f = list(e); f[i] -= 1; f = tuple(f)
            if f in r: r.remove(f)
            else: r.add(f)
    return r
def mp_pow2(a, j):
    # Frobenius: (sum m)^(2^j) = sum m^(2^j) in char 2
    return {tuple(x << j for x in e) for e in a}
def mp_eval(F, a, pt):
    s = 0
    for e in a:
        t = 1
        for x, ex in zip(pt, e):
            if ex: t = F.mul(t, F.pow(x, ex))
        s ^= t
    return s

# ---- F_4 = {0,1,w,w+1} encoded 0..3 (bit0 = 1, bit1 = w), w^2 = w+1
def f4mul(a, b):
    a0, a1 = a & 1, a >> 1; b0, b1 = b & 1, b >> 1
    # (a0 + a1 w)(b0 + b1 w) = a0b0 + (a0b1+a1b0) w + a1b1 (w+1)
    c0 = (a0 & b0) ^ (a1 & b1); c1 = (a0 & b1) ^ (a1 & b0) ^ (a1 & b1)
    return c0 | (c1 << 1)
F4INV = {1: 1, 2: 3, 3: 2}
def f4conj(a):  # w -> w^2 = w+1
    a0, a1 = a & 1, a >> 1
    return (a0 ^ a1) | (a1 << 1)

def p_trim(p):
    while p and p[-1] == 0: p.pop()
    return p
def p_add(a, b):
    n = max(len(a), len(b)); r = [0] * n
    for i, x in enumerate(a): r[i] ^= x
    for i, x in enumerate(b): r[i] ^= x
    return p_trim(r)
def p_mul(a, b):
    if not a or not b: return []
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y: r[i + j] ^= f4mul(x, y)
    return p_trim(r)
def p_divmod(a, b):
    a = a[:]; q = [0] * max(0, len(a) - len(b) + 1)
    ib = F4INV[b[-1]]
    while len(a) >= len(b) and a:
        c = f4mul(a[-1], ib); s = len(a) - len(b); q[s] = c
        for j, y in enumerate(b):
            a[s + j] ^= f4mul(c, y)
        p_trim(a)
    return p_trim(q), a
def p_gcd(a, b):
    while b:
        _, r = p_divmod(a, b); a, b = b, r
    if a:
        iv = F4INV[a[-1]]; a = [f4mul(iv, x) for x in a]
    return a
def p_conj(a): return [f4conj(x) for x in a]
def deg(a): return len(a) - 1
def f2poly_to_f4(bits):  # int bitmask over F_2 -> list
    r = []
    while bits: r.append(bits & 1); bits >>= 1
    return r
