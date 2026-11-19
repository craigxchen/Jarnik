"""Independent GF(2^d) and GF(2^d)[t] arithmetic (log/antilog tables, dense coefficient lists).

Polynomials: lists of field elements, index = degree, trailing zeros stripped.
Field elements: ints 0..q-1 (bit vectors w.r.t. a primitive polynomial found here by search).
"""


def find_primitive(d):
    """Smallest primitive polynomial of degree d over F_2 (as int with bit d set)."""
    q = 1 << d
    for poly in range(q | 1, 2 * q, 2):
        # generate powers of x mod poly, check order q-1
        x, seen = 1, set()
        ok = True
        for k in range(q - 1):
            if x in seen:
                ok = False
                break
            seen.add(x)
            x <<= 1
            if x & q:
                x ^= poly
        if ok and x == 1 and len(seen) == q - 1:
            return poly
    raise ValueError


class GF:
    def __init__(self, d):
        self.d, self.q = d, 1 << d
        self.poly = find_primitive(d)
        q = self.q
        self.exp = [0] * (2 * q)
        self.log = [None] * q
        x = 1
        for k in range(q - 1):
            self.exp[k] = x
            self.log[x] = k
            x <<= 1
            if x & q:
                x ^= self.poly
        for k in range(q - 1, 2 * q):
            self.exp[k] = self.exp[k - (q - 1)]

    def mul(self, a, b):
        if a == 0 or b == 0:
            return 0
        return self.exp[self.log[a] + self.log[b]]

    def inv(self, a):
        assert a
        return self.exp[(self.q - 1 - self.log[a]) % (self.q - 1)]

    def pw(self, a, n):
        if a == 0:
            return 0 if n > 0 else 1
        return self.exp[(self.log[a] * n) % (self.q - 1)]

    def tr(self, a):
        s, y = 0, a
        for _ in range(self.d):
            s ^= y
            y = self.mul(y, y)
        assert s in (0, 1)
        return s

    # ---- polynomials ----
    @staticmethod
    def strip(p):
        p = list(p)
        while p and p[-1] == 0:
            p.pop()
        return p

    @staticmethod
    def deg(p):
        return len(p) - 1  # -1 for zero polynomial

    def add(self, p, r):
        n = max(len(p), len(r))
        out = [0] * n
        for i, c in enumerate(p):
            out[i] ^= c
        for i, c in enumerate(r):
            out[i] ^= c
        return self.strip(out)

    def pmul(self, p, r):
        if not p or not r:
            return []
        out = [0] * (len(p) + len(r) - 1)
        for i, c in enumerate(p):
            if c:
                lc = self.log[c]
                for j, e in enumerate(r):
                    if e:
                        out[i + j] ^= self.exp[lc + self.log[e]]
        return self.strip(out)

    def scal(self, c, p):
        return self.strip([self.mul(c, e) for e in p])

    def pdivmod(self, p, r):
        p = list(p)
        r = self.strip(r)
        assert r
        dr = len(r) - 1
        ilc = self.inv(r[-1])
        quo = [0] * max(0, len(p) - dr)
        for k in range(len(p) - 1, dr - 1, -1):
            c = p[k]
            if c:
                f = self.mul(c, ilc)
                quo[k - dr] = f
                for j, e in enumerate(r):
                    if e:
                        p[k - dr + j] ^= self.mul(f, e)
        return self.strip(quo), self.strip(p[:dr])

    def pgcd(self, p, r):
        p, r = self.strip(p), self.strip(r)
        while r:
            p, r = r, self.pdivmod(p, r)[1]
        if p:
            p = self.scal(self.inv(p[-1]), p)
        return p

    def peval(self, p, x):
        s = 0
        for c in reversed(p):
            s = self.mul(s, x) ^ c
        return s

    def from_roots(self, roots):
        p = [1]
        for x in roots:
            p = self.pmul(p, [x, 1])  # t - x = t + x in char 2
        return p
