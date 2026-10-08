# Golden-unit analogue of the eight-point parity-class template:
# J_r = f_{r+1} + i f_r  (norm f_{2r+1}),  blocks J_{n+s}, s in shifts,
# rows: sign words in the odd (or even) parity class, prefactor aligning leading phases.
import itertools, math, sys
from gpoly import gmul, ggcd, gdiv

def fib(r):
    a, b = 0, 1
    for _ in range(r):
        a, b = b, a + b
    return a

def J(r): return (fib(r+1), fib(r))
def conj(z): return (z[0], -z[1])
F = (1, 2)

def gpow(z, e):
    p = (1, 0)
    for _ in range(e):
        p = gmul(p, z)
    return p

def rows(n, shifts, parity, block=J):
    bl = [block(n+s) for s in shifts]
    K = len(shifts)
    pts = []; words = []
    for w in itertools.product([0, 1], repeat=K):
        p = sum(w)
        if p % 2 != parity: continue
        # leading phase (2p-K)*arg(a), arg a = arg(F)/2 ; compensate with F^(K-p) Fbar^p up to common
        c = gmul(gpow(F, max(0, (K - 2*p)//2) if False else 0), (1, 0))
        # general: need prefactor with phase -(2p-K)*arg(F)/2 = (K-2p)/2 * arg F; K-2p even here
        e = (K - 2*p)//2
        c = gmul(gpow(F, max(e, 0)), gpow(conj(F), max(-e, 0)))
        # equalize norms: multiply by (F Fbar)^(|emax|-|e|) = 5^(...)
        pts.append((c, w))
    emax = max(abs((K - 2*sum(w))//2) for _, w in pts)
    out = []
    for c, w in pts:
        e = abs((K - 2*sum(w))//2)
        z = gmul(c, (5**(emax - e), 0))
        for b, s in zip(bl, w):
            z = gmul(z, b if s else conj(b))
        out.append(z)
    return out, [w for _, w in pts]

def spans(pts):
    g = pts[0]
    for p in pts[1:]: g = ggcd(g, p)
    P = [gdiv(p, g) for p in pts]
    N = P[0][0]**2 + P[0][1]**2
    assert all(p[0]**2 + p[1]**2 == N for p in P)
    p0c = conj(P[0]); angs = []
    for p in P:
        q = gmul(p, p0c)
        q = max([q, (-q[1], q[0]), (-q[0], -q[1]), (q[1], -q[0])], key=lambda t: t[0])
        angs.append(math.atan(q[1]/q[0]))
    r4 = math.isqrt(math.isqrt(N))
    r4f = float(r4) if r4 < 10**300 else None
    return angs, N, g

if __name__ == "__main__":
    shifts = [0, 1, 2, 3]
    for parity in (1, 0):
        for n in [40, 41, 80, 81, 160, 161]:
            pts, words = rows(n, shifts, parity)
            angs, N, g = spans(pts)
            lr = math.log(N)/4  # log sqrt(R)
            span = max(angs) - min(angs)
            C = span*math.exp(lr)
            print("parity", parity, "n", n, "k", len(set(pts)), "C=%.6f" % C, "gcdnorm", g[0]**2 + g[1]**2)
