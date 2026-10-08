# Golden-unit (Fibonacci) fixed-shift templates: blocks J_r = f_{r+1} + i f_r, |J_r|^2 = f_{2r+1},
# limiting phase theta=arg(F)/2, F=1+2i.  A row is an exponent vector a (0<=a_s<=e_s) on shifts s;
# z = c_a * prod_s J_{n+s}^{a_s} conj(J_{n+s})^{e_s-a_s}, with prefactor c_a = F^al Fbar^be,
# al-be = (E-2|a|)/2, al+be fixed (rows of one parity class), so all rows have equal norm and
# equal limiting argument.  Exact evaluation, division by the cluster's Gaussian gcd, and the
# minimal normalized span for every sub-cluster size.
import itertools, math
from gpoly import gmul, ggcd, gdiv

_fib = [0, 1]
def fib(r):
    while len(_fib) <= r: _fib.append(_fib[-1] + _fib[-2])
    return _fib[r]
def J(r): return (fib(r+1), fib(r))
def conj(z): return (z[0], -z[1])
F = (1, 2)
def gpow(z, e):
    p = (1, 0)
    for _ in range(e): p = gmul(p, z)
    return p

def evaluate(n, shifts, widths, rows):
    E = sum(widths)
    assert E % 2 == 0
    es = [(E - 2*sum(a))//2 for a in rows]
    par = {sum(a) % 2 for a in rows}
    assert len(par) == 1
    tot = max(abs(e) for e in es)
    pts = []
    bl = [J(n+s) for s in shifts]
    for a, e in zip(rows, es):
        al, be = max(e, 0), max(-e, 0)
        extra = (tot - abs(e))//2  # multiply by 5^extra to equalize norm
        z = gmul(gmul(gpow(F, al), gpow(conj(F), be)), (5**extra, 0))
        for b, w, x in zip(bl, widths, a):
            z = gmul(z, gmul(gpow(b, x), gpow(conj(b), w - x)))
        pts.append(z)
    return pts

def subspans(pts):
    g = pts[0]
    for p in pts[1:]: g = ggcd(g, p)
    P = [gdiv(p, g) for p in pts]
    N = P[0][0]**2 + P[0][1]**2
    assert all(p[0]**2 + p[1]**2 == N for p in P)
    if len(set(P)) < len(P): return None, N, g
    p0c = conj(P[0]); angs = []
    for p in P:
        q = gmul(p, p0c)
        q = max([q, (-q[1], q[0]), (-q[0], -q[1]), (q[1], -q[0])], key=lambda t: t[0])
        angs.append(math.atan(q[1]/q[0]))
    order = sorted(range(len(angs)), key=lambda i: angs[i])
    a = [angs[i] for i in order]
    lr = math.exp(math.log(N)/4)
    best = {}
    for k in range(2, len(a)+1):
        best[k] = min((a[i+k-1]-a[i])*lr for i in range(len(a)-k+1))
    return best, N, g

def template_best(shifts, widths, rows, nrange):
    res = {}
    for n in nrange:
        pts = evaluate(n, shifts, widths, rows)
        b, N, g = subspans(pts)
        if b is None: continue
        for k, v in b.items():
            if k not in res or v < res[k][0]: res[k] = (v, n)
    return res

if __name__ == "__main__":
    rows = [a for a in itertools.product([0,1], repeat=4) if sum(a) % 2 == 1]
    print(template_best([0,1,2,3], [1,1,1,1], rows, range(60, 90)))
