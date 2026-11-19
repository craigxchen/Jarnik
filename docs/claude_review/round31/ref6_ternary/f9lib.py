"""Own arithmetic: F_3, F_9 = F_3[i]/(i^2+1), F_9[t] as pairs of int64 numpy arrays (low degree first),
F_q (q=3^r) as lists of ints via a log table. No external finite-field library."""
import numpy as np

P = 3

def trim(a):
    nz = np.nonzero(a)[0]
    return a[:nz[-1]+1] if len(nz) else a[:0]

class F9P:
    """polynomial re + i*im over F_3, arrays reduced mod 3"""
    __slots__ = ('re', 'im')
    def __init__(self, re, im=None):
        re = np.asarray(re, dtype=np.int64) % 3
        im = np.zeros(len(re), np.int64) if im is None else np.asarray(im, dtype=np.int64) % 3
        n = max(len(re), len(im))
        r = np.zeros(n, np.int64); r[:len(re)] = re
        s = np.zeros(n, np.int64); s[:len(im)] = im
        nz = np.nonzero(r | s)[0]
        k = nz[-1]+1 if len(nz) else 0
        self.re, self.im = r[:k], s[:k]
    def deg(self):
        return len(self.re) - 1
    def iszero(self):
        return len(self.re) == 0
    def __add__(s, o):
        n = max(len(s.re), len(o.re))
        r = np.zeros(n, np.int64); r[:len(s.re)] += s.re; r[:len(o.re)] += o.re
        m = np.zeros(n, np.int64); m[:len(s.im)] += s.im; m[:len(o.im)] += o.im
        return F9P(r, m)
    def __neg__(s):
        return F9P(-s.re, -s.im)
    def __sub__(s, o):
        return s + (-o)
    def __mul__(s, o):
        if isinstance(o, tuple):  # scalar (a,b) = a+bi
            a, b = o
            return F9P(a*s.re - b*s.im, a*s.im + b*s.re)
        if s.iszero() or o.iszero():
            return F9P([])
        ac = conv(s.re, o.re); bd = conv(s.im, o.im)
        ad = conv(s.re, o.im); bc = conv(s.im, o.re)
        return F9P(ac - bd, ad + bc)
    def conj(s):
        return F9P(s.re, -s.im)
    def lead(s):
        return (int(s.re[-1]), int(s.im[-1]))
    def __eq__(s, o):
        return np.array_equal(s.re, o.re) and np.array_equal(s.im, o.im)
    def shift(s, k):
        return F9P(np.concatenate([np.zeros(k, np.int64), s.re]), np.concatenate([np.zeros(k, np.int64), s.im]))
    def is_f3(s):
        return not np.any(s.im)

def conv(a, b):
    """exact integer convolution of small (|x|<=2 after centering) vectors, result mod 3"""
    if len(a) == 0 or len(b) == 0:
        return np.zeros(0, np.int64)
    if min(len(a), len(b)) < 600:
        return np.convolve(a % 3, b % 3) % 3
    ac = ((a + 1) % 3) - 1; bc = ((b + 1) % 3) - 1   # centered in {-1,0,1}
    n = len(a) + len(b) - 1
    N = 1 << (n - 1).bit_length()
    fa = np.fft.rfft(ac.astype(np.float64), N); fb = np.fft.rfft(bc.astype(np.float64), N)
    c = np.fft.irfft(fa * fb, N)[:n]
    ci = np.rint(c)
    err = np.max(np.abs(c - ci)) if n else 0
    assert err < 0.2, err
    return ci.astype(np.int64) % 3

# F_9 scalars as pairs
def f9mul(x, y):
    return ((x[0]*y[0] - x[1]*y[1]) % 3, (x[0]*y[1] + x[1]*y[0]) % 3)
def f9inv(x):
    for a in range(3):
        for b in range(3):
            if f9mul(x, (a, b)) == (1, 0):
                return (a, b)
    raise ZeroDivisionError

def pw(a, e):
    r = F9P([1]); b = a
    while e:
        if e & 1: r = r * b
        e >>= 1
        if e: b = b * b
    return r

def divmod9(a, b):
    """Euclidean division in F_9[t]"""
    ar, ai = a.re.copy(), a.im.copy()
    db = b.deg(); inv = f9inv(b.lead())
    br, bi = b.re, b.im
    qn = max(len(ar) - db, 0)
    qr = np.zeros(qn, np.int64); qi = np.zeros(qn, np.int64)
    for k in range(len(ar) - 1, db - 1, -1):
        c = (int(ar[k]) % 3, int(ai[k]) % 3)
        if c == (0, 0): continue
        f = f9mul(c, inv)
        qr[k-db], qi[k-db] = f
        s = k - db
        ar[s:k+1] -= f[0]*br - f[1]*bi
        ai[s:k+1] -= f[0]*bi + f[1]*br
        ar[s:k+1] %= 3; ai[s:k+1] %= 3
    return F9P(qr, qi), F9P(ar[:db], ai[:db])

def gcd9(a, b):
    while not b.iszero():
        _, r = divmod9(a, b)
        a, b = b, r
    if a.iszero(): return a
    return a * f9inv(a.lead())

class Fq:
    """F_{3^r} via a primitive-element log table; elements are ints 0..q-1 (base-3 digit vectors)"""
    def __init__(self, r):
        self.r = r; q = self.q = 3**r
        # find irreducible primitive modulus x^r + ... by brute force
        import itertools
        for coeffs in itertools.product(range(3), repeat=r):
            if coeffs[0] == 0: continue
            # multiply-by-x map on digit vectors
            def mulx(v):
                d = [0] + v[:-1]; top = v[-1]
                return [(d[j] - top*coeffs[j]) % 3 for j in range(r)]
            v = [1] + [0]*(r-1); seen = []
            for _ in range(q-1):
                seen.append(sum(c*3**j for j, c in enumerate(v)))
                v = mulx(v)
            if len(set(seen)) == q-1:
                break
        self.exp = seen; self.log = {e: k for k, e in enumerate(seen)}
        def add(x, y):
            return sum(((x//3**j + y//3**j) % 3)*3**j for j in range(r))
        self.addt = [[add(x, y) for y in range(q)] for x in range(q)]
        self.neg = [sum(((-(x//3**j)) % 3)*3**j for j in range(r)) for x in range(q)]
        self.chi = [0] + [1 if self.log[x] % 2 == 0 else -1 for x in range(1, q)]
    def smul(self, c, x):   # c in F_3
        c %= 3
        return 0 if c == 0 else (x if c == 1 else self.neg[x])
    def mul(self, x, y):
        if x == 0 or y == 0: return 0
        return self.exp[(self.log[x] + self.log[y]) % (self.q - 1)]
    def frob(self, x, e):   # x^(3^e)
        if x == 0: return 0
        return self.exp[(self.log[x]*3**e) % (self.q - 1)]
