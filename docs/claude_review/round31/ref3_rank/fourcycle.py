"""Independent referee check of the balanced four-cycle of astra_free_III III.5.
Uses sympy for the symbolic quartics, and own Gaussian-integer arithmetic for points."""
import sys
from math import gcd, lcm, isqrt
from itertools import combinations
from fractions import Fraction
import sympy as sp
from mpmath import mp, mpf, atan2, log as mlog
mp.dps = 80

n = sp.symbols('n', integer=True)
roots = [(0,6),(-135,-36),(-135,-51),(-135,75),(-145,81),(-117,-45),(-108,-30)]
rows = [set(), {0,1,2,4}, {0,1,3,5}, {0,2,3,6}]
H = [n - (a + sp.I*b) for a,b in roots]
Hb = [n - (a - sp.I*b) for a,b in roots]

# 1. edge quartics B_ij
edges = list(combinations(range(4),2))
B = {}
for (i,j) in edges:
    e = 1
    for k in rows[i]-rows[j]: e *= H[k]
    for k in rows[j]-rows[i]: e *= Hb[k]
    e = sp.expand(e)
    re = sp.Poly(sp.re(sp.expand(e.subs(n, sp.Symbol('x', real=True)))).subs(sp.Symbol('x', real=True), n), n)
    im = sp.expand(sp.im(e.subs(n, sp.Symbol('x', real=True))))
    B[(i,j)] = (re, im)
    print((i,j), 'deg', re.degree(), 'lead', re.LC(), 'Im =', im)

bs = [abs(int(B[e][1])) for e in edges]
print('|b|', bs)
assert bs == [19966500,16110900,14871600,3855600,5094900,1239300]
cyc = [(0,1),(0,2),(1,3),(2,3)]
P = lcm(*[abs(int(B[e][1])) for e in cyc])
print('lcm of cycle b', P, P == 86434978500)
n0 = 79658575767
for e in cyc:
    re = int(B[e][0].eval(n0)); b = abs(int(B[e][1]))
    g = gcd(re, b)
    print(e, 'gcd', g, 'b/g', b//g, 'Re/g parity', (re//g) % 2)

# 2. actual points at several n on the progression
def gmul(z,w): return (z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0])
def gconj(z): return (z[0], -z[1])
def gnorm(z): return z[0]*z[0]+z[1]*z[1]
def gdivmod_round(z,w):
    num = gmul(z, gconj(w)); d = gnorm(w)
    q = ((2*num[0]+d)//(2*d), (2*num[1]+d)//(2*d))
    r = (z[0]-gmul(q,w)[0], z[1]-gmul(q,w)[1])
    return q, r
def ggcd(z,w):
    while gnorm(w):
        _, r = gdivmod_round(z,w); z, w = w, r
    return z
def gexact(z,w):
    num = gmul(z, gconj(w)); d = gnorm(w)
    assert num[0]%d == 0 and num[1]%d == 0
    return (num[0]//d, num[1]//d)

# bounded-prime set: primes dividing norms of all root differences (incl. conjugates)
allr = [(a,b) for a,b in roots] + [(a,-b) for a,b in roots]
K = 1
for z,w in combinations(allr,2):
    K = lcm(K, gnorm((z[0]-w[0], z[1]-w[1])))
Kprimes = sorted(sp.factorint(K).keys())
print('bounded primes', Kprimes)

def smooth_part(x, ps):
    s = 1
    for p in ps:
        while x % p == 0: x//=p; s*=p
    return s

for t in [0, 1, 7, 10**6, 10**20]:
    nn = n0 + t*P
    h = [(nn-a, -b) for a,b in roots]
    pts = []
    for S in rows:
        z = (1,0)
        for k in range(7):
            z = gmul(z, h[k] if k in S else gconj(h[k]))
        pts.append(z)
    G = pts[0]
    for z in pts[1:]: G = ggcd(G, z)
    pts = [gexact(z, G) for z in pts]
    N = gnorm(pts[0]); assert all(gnorm(z)==N for z in pts)
    assert len(set(pts)) == 4
    offs = []; halfs = []
    for (i,j) in edges:
        w = gmul(pts[i], gconj(pts[j]))
        # half-angle tangent of arg(z_i/z_j): tan(theta/2) = Im w / (N + Re w)
        f = Fraction(w[1], N + w[0])
        s_, t_ = f.denominator, f.numerator  # A = s + i t primitive, s>0
        assert gcd(s_, t_) == 1 and (s_*s_+t_*t_) % 2 == 1 and N % (s_*s_+t_*t_) == 0
        # check A/conj(A) equals z_i/z_j exactly: z_i * conj(A) * z_j-bar ... cross multiply
        A = (s_, t_)
        lhs = gmul(pts[i], gconj(A)); rhs = gmul(pts[j], A)
        assert lhs == rhs, 'ratio mismatch'
        offs.append(abs(t_))
    # arc: angles relative to point 0
    ang = [atan2(mpf(gmul(z,gconj(pts[0]))[1]), mpf(gmul(z,gconj(pts[0]))[0])) for z in pts]
    width = max(ang)-min(ang)
    Rq = mpf(N)**mpf(0.25)  # sqrt(R), R = sqrt(N)
    arcconst = width*mpf(N)**0.5 / Rq
    logN = mlog(mpf(N))
    # cut weights: column k weight = log N(H_k) minus bounded-prime part
    cw = []
    for k in range(7):
        x = gnorm(h[k]); cw.append(float(mlog(mpf(x // smooth_part(x, Kprimes)))))
    print(f't={t}: digits N={len(str(N))}, N(G)={gnorm(G)}, offsets={offs}, arc const={float(arcconst):.3e},'
          f' 1/128={1/128:.3e}, width*N^(2/7)={float(width*mpf(N)**(mpf(2)/7)):.3f}')
    print('   logN/7=%.2f  residual col weights' % float(logN/7), [round(c,2) for c in cw],
          ' bounded part of logN = %.2f' % (float(logN) - sum(cw)))

# extra checks at t=0 and t=1
for t in [0,1]:
    nn = n0 + t*P
    h = [(nn-a, -b) for a,b in roots]
    pts = []
    for S in rows:
        z = (1,0)
        for k in range(7):
            z = gmul(z, h[k] if k in S else gconj(h[k]))
        pts.append(z)
    G = pts[0]
    for z in pts[1:]: G = ggcd(G, z)
    pts = [gexact(z, G) for z in pts]
    G2 = pts[0]
    for z in pts[1:]: G2 = ggcd(G2, z)
    N = gnorm(pts[0])
    bad = [p for p in sp.primerange(2, 200000) if N % p == 0 and p % 4 != 1]
    print('t',t,'gcd of normalized points has norm', gnorm(G2), 'N mod 4 =', N%4, 'non-split small primes:', bad)
    print('  small prime factors of N:', [(p, sp.multiplicity(p,N)) for p in sp.primerange(2,200000) if N%p==0][:30])
