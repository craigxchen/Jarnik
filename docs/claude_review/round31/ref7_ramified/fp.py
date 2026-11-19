# Own F_p[t] arithmetic (numpy int64 arrays, low degree first) and the quadratic order F_p[t,eta], eta^2 = Delta.
import numpy as np

P = 3

def trim(a):
    a = np.asarray(a, dtype=np.int64) % P
    nz = np.nonzero(a)[0]
    if len(nz) == 0:
        return np.zeros(0, dtype=np.int64)
    return a[:nz[-1] + 1].copy()

def deg(a):
    return len(a) - 1  # -1 for zero

def add(a, b):
    n = max(len(a), len(b))
    c = np.zeros(n, dtype=np.int64)
    c[:len(a)] += a
    c[:len(b)] += b
    return trim(c)

def neg(a):
    return trim(-a)

def sub(a, b):
    return add(a, neg(b))

def scal(c, a):
    return trim(c * a)

def mul(a, b):
    if len(a) == 0 or len(b) == 0:
        return np.zeros(0, dtype=np.int64)
    if min(len(a), len(b)) < 64:
        return trim(np.convolve(a, b))
    # exact via FFT in float with splitting not needed: entries < 3, sums < 4*len
    n = len(a) + len(b) - 1
    N = 1 << (n - 1).bit_length()
    fa = np.fft.rfft(a.astype(np.float64), N)
    fb = np.fft.rfft(b.astype(np.float64), N)
    c = np.rint(np.fft.irfft(fa * fb, N)[:n]).astype(np.int64)
    return trim(c)

def mono(k, c=1):
    a = np.zeros(k + 1, dtype=np.int64); a[k] = c
    return trim(a)

def const(c):
    return trim(np.array([c]))

def powp(a, e):
    r = const(1)
    b = a
    while e:
        if e & 1: r = mul(r, b)
        b = mul(b, b); e >>= 1
    return r

def inv(x):
    return pow(int(x), P - 2, P)

def divmod_(a, b):
    a = a.copy(); db = deg(b); ib = inv(b[-1])
    if deg(a) < db:
        return np.zeros(0, dtype=np.int64), trim(a)
    q = np.zeros(deg(a) - db + 1, dtype=np.int64)
    for i in range(deg(a) - db, -1, -1):
        c = (a[i + db] * ib) % P
        if c:
            q[i] = c
            a[i:i + db + 1] = (a[i:i + db + 1] - c * b) % P
    return trim(q), trim(a)

def gcd(a, b):
    a = trim(a); b = trim(b)
    while len(b):
        _, r = divmod_(a, b)
        a, b = b, r
    if len(a):
        a = scal(inv(a[-1]), a)
    return a

def deriv(a):
    if len(a) <= 1:
        return np.zeros(0, dtype=np.int64)
    return trim(a[1:] * np.arange(1, len(a)))

def squarefree(a):
    return deg(gcd(a, deriv(a))) == 0

def max_exponent(a):
    """largest multiplicity of an irreducible factor (square-free decomposition, char p aware)."""
    # Yun-like in char p: use repeated gcd with derivative and p-th roots.
    a = trim(a)
    if deg(a) <= 0:
        return 0
    best = 0
    # factor a = prod a_i^i ; do generic algorithm
    def sqf(f):
        # returns dict multiplicity->poly
        res = {}
        i = 1
        c = gcd(f, deriv(f))
        w, _ = divmod_(f, c)
        while deg(w) > 0:
            y = gcd(w, c)
            z, _ = divmod_(w, y)
            if deg(z) > 0: res[i] = z
            i += 1
            w = y
            c, _ = divmod_(c, y)
        if deg(c) > 0:
            # c is a p-th power
            croot = trim(c[::P])
            sub_ = sqf(croot)
            for k, v in sub_.items():
                res[k * P] = v
        return res
    d = sqf(a)
    return max(d.keys()), d

# ---------- quadratic order elements: (A,B) = A + B eta ----------
class Order:
    def __init__(self, Delta):
        self.D = trim(Delta)
    def mul(self, x, y):
        A, B = x; C, E = y
        return (add(mul(A, C), mul(self.D, mul(B, E))), add(mul(A, E), mul(B, C)))
    def add(self, x, y):
        return (add(x[0], y[0]), add(x[1], y[1]))
    def sub(self, x, y):
        return (sub(x[0], y[0]), sub(x[1], y[1]))
    def conj(self, x):
        return (x[0], neg(x[1]))
    def norm(self, x):
        return sub(mul(x[0], x[0]), mul(self.D, mul(x[1], x[1])))
    def pole2(self, x):
        """twice the pole order at infinity: max(2 deg A, 2 deg B + deg Delta); -inf for 0."""
        A, B = x
        vals = []
        if len(A): vals.append(2 * deg(A))
        if len(B): vals.append(2 * deg(B) + deg(self.D))
        return max(vals) if vals else -10**9
    def pow(self, x, e):
        r = (const(1), np.zeros(0, dtype=np.int64)); b = x
        while e:
            if e & 1: r = self.mul(r, b)
            b = self.mul(b, b); e >>= 1
        return r
