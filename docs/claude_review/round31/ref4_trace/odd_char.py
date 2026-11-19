"""Odd-characteristic analogue of the trace construction: F_p-hyperplanes of F_q, q = p^d.

Marks: hyperplanes H_a = ker(x -> Tr(a x)), a in F_q^*/F_p^*  ((q-1)/(p-1) marks).
Primes: t - x, x in F_q^*.  Row factor n_a = prod_{x in H_a, x != 0} (t-x), degree q/p - 1.
Column P_a = (f_a, 1), f_a = (t^(q/p) - L_a)/t with L_a the monic subspace polynomial of H_a.
We check: f_a(x) = x^(q/p-1) on H_a\\0, f_a - f_b = const * gcd(n_a,n_b), and report
pair separation and degree budget, which FAIL for p odd: |H_a xor H_b| = 2q(p-1)/p^2 <= q/2 - (q-1)/2...
"""
import sys
from fractions import Fraction
from itertools import combinations


def irreducible(p, d):
    # find monic irreducible of degree d over F_p by brute force (coeff lists low->high)
    from itertools import product
    for tail in product(range(p), repeat=d):
        f = list(tail) + [1]
        if f[0] == 0:
            continue
        # check no roots in extensions: use Rabin-free brute force: x^(p^d) = x mod f and no smaller
        if is_irred(f, p, d):
            return f
    raise ValueError


def polmulmod(a, b, f, p):
    d = len(f) - 1
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] = (out[i + j] + x * y) % p
    for k in range(len(out) - 1, d - 1, -1):
        c = out[k]
        if c:
            for j in range(d + 1):
                out[k - d + j] = (out[k - d + j] - c * f[j]) % p
    out = out[:d] + [0] * (d - len(out[:d]))
    return out


def is_irred(f, p, d):
    # f irreducible iff x^(p^d) = x mod f and gcd(x^(p^(d/r)) - x, f) = 1 for prime r | d.
    # Simpler for small sizes: count roots via exhaustive: build field candidates; use order test.
    q = p ** d
    # brute force: f has no factor of degree <= d/2
    from itertools import product
    def pdivides(g, f):
        # does monic g divide f over F_p
        r = f[:]
        dg = len(g) - 1
        for k in range(len(r) - 1, dg - 1, -1):
            c = r[k]
            if c:
                for j in range(dg + 1):
                    r[k - dg + j] = (r[k - dg + j] - c * g[j]) % p
        return not any(r[:dg])
    for e in range(1, d // 2 + 1):
        for tail in product(range(p), repeat=e):
            if pdivides(list(tail) + [1], f):
                return False
    return True


class GFq:
    def __init__(self, p, d):
        self.p, self.d, self.q = p, d, p ** d
        f = irreducible(p, d)
        q = self.q
        # elements as tuples -> int index base p
        def toint(v):
            return sum(c * p ** i for i, c in enumerate(v))
        def tovec(n):
            return [(n // p ** i) % p for i in range(d)]
        self.toint, self.tovec = toint, tovec
        self.addt = [[toint([(x + y) % p for x, y in zip(tovec(a), tovec(b))]) for b in range(q)] for a in range(q)]
        self.negt = [toint([(-x) % p for x in tovec(a)]) for a in range(q)]
        # find primitive element
        for g in range(1, q):
            seen, x = [], toint([1] + [0] * (d - 1))
            xs = set()
            for k in range(q - 1):
                if x in xs:
                    break
                xs.add(x); seen.append(x)
                x = toint(polmulmod(tovec(x), tovec(g), f, p))
            if len(xs) == q - 1:
                break
        self.exp = seen + seen
        self.log = {v: k for k, v in enumerate(seen)}

    def mul(self, a, b):
        if a == 0 or b == 0:
            return 0
        return self.exp[self.log[a] + self.log[b]]

    def pw(self, a, n):
        if a == 0:
            return 0
        return self.exp[(self.log[a] * n) % (self.q - 1)]

    def add(self, a, b):
        return self.addt[a][b]

    def tr(self, a):
        s, y = 0, a
        for _ in range(self.d):
            s = self.add(s, y)
            y = self.pw(y, self.p)
        return s

    # polys
    def strip(self, P):
        P = list(P)
        while P and P[-1] == 0:
            P.pop()
        return P

    def padd(self, A, B):
        n = max(len(A), len(B))
        A = A + [0] * (n - len(A)); B = B + [0] * (n - len(B))
        return self.strip([self.add(x, y) for x, y in zip(A, B)])

    def pneg(self, A):
        return [self.negt[x] for x in A]

    def pmul(self, A, B):
        out = [0] * (len(A) + len(B) - 1)
        for i, x in enumerate(A):
            for j, y in enumerate(B):
                out[i + j] = self.add(out[i + j], self.mul(x, y))
        return self.strip(out)

    def from_roots(self, R):
        P = [1]
        for x in R:
            P = self.pmul(P, [self.negt[x], 1])
        return P

    def peval(self, P, x):
        s = 0
        for c in reversed(P):
            s = self.add(self.mul(s, x), c)
        return s

    def pdivmod(self, A, B):
        A = list(A); dB = len(B) - 1
        inv = self.pw(B[-1], self.q - 2)
        Q = [0] * max(0, len(A) - dB)
        for k in range(len(A) - 1, dB - 1, -1):
            c = A[k]
            if c:
                f = self.mul(c, inv)
                Q[k - dB] = f
                for j, e in enumerate(B):
                    A[k - dB + j] = self.add(A[k - dB + j], self.negt[self.mul(f, e)])
        return self.strip(Q), self.strip(A[:dB])

    def pgcd(self, A, B):
        A, B = self.strip(A), self.strip(B)
        while B:
            A, B = B, self.pdivmod(A, B)[1]
        return A


def run(p, d):
    F = GFq(p, d)
    q = F.q
    nz = list(range(1, q))
    # representatives of F_q^*/F_p^*: canonical = min over scalar multiples
    Fp = [F.toint([c] + [0] * (d - 1)) for c in range(1, p)]
    reps = sorted({min(F.mul(c, a) for c in Fp) for a in nz})
    assert len(reps) == (q - 1) // (p - 1)
    rows = {}
    for a in reps:
        H = [x for x in nz if F.tr(F.mul(a, x)) == 0]
        assert len(H) == q // p - 1
        LH = F.from_roots([0] + H)               # monic subspace polynomial
        assert all(c == 0 for i, c in enumerate(LH) if i not in [p ** j for j in range(d)])  # p-linearized
        num = F.padd([0] * (q // p) + [1], F.pneg(LH))
        assert num[0] == 0
        fa = F.strip(num[1:])
        assert len(fa) - 1 == q // p ** 2 - 1
        for x in H:
            assert F.peval(fa, x) == F.pw(x, q // p - 1)
        rows[a] = (set(H), fa, F.from_roots(H))
    L = q - 1
    worst_excess, sep_vals = None, set()
    for a, b in combinations(reps, 2):
        Ha, fa, na = rows[a]; Hb, fb, nb = rows[b]
        D = F.padd(fa, F.pneg(fb))
        g = F.pgcd(na, nb)
        Qt, R = F.pdivmod(D, g)
        assert R == [] and len(Qt) == 1, (a, b)
        degD = len(D) - 1
        budget = Fraction(len(na) - 1 + len(nb) - 1, 2) - Fraction(L, 4)
        ex = degD - budget
        worst_excess = ex if worst_excess is None else max(worst_excess, ex)
        sep_vals.add(len(Ha ^ Hb))
    m = len(reps)
    print(f"p={p}, q={q}: m={m} marks, deg D = deg g = {q//p**2-1}, cofactor constant on all pairs;"
          f" separation {sorted(sep_vals)} vs L/2={L/2} -> {'OK' if min(sep_vals) > L/2 else 'FAILS'};"
          f" degree excess over critical budget = {worst_excess} ({float(worst_excess):.3f})")


if __name__ == '__main__':
    for p, d in [(3, 2), (3, 3), (3, 4), (5, 2), (5, 3), (7, 2)]:
        run(p, d)
