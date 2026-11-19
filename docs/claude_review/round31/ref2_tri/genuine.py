"""Referee check (ref2_tri): (II.1), its reverse form, the profile condition and the
determinant (Question A, k=3) dictionary on genuine triples of lattice points.
Own code; standard library + sympy.factorint only."""
import sys, random, math
from math import gcd, isqrt
from itertools import combinations
from functools import reduce
from sympy import factorint

def mul(z, w): return (z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0])
def cj(z): return (z[0], -z[1])
def nm(z): return z[0]*z[0]+z[1]*z[1]
def det(p, q): return p[0]*q[1]-p[1]*q[0]

def ggcd(z, w):
    while w != (0, 0):
        n = nm(w); x, y = mul(z, cj(w))
        q = ((2*x+n)//(2*n), (2*y+n)//(2*n))
        t = mul(q, w); z, w = w, (z[0]-t[0], z[1]-t[1])
    return z

def gprime(p):
    if p == 2: return (1, 1)
    # find x with x^2 = -1 mod p
    for g in range(2, p):
        x = pow(g, (p-1)//4, p)
        if x*x % p == p-1: break
    return ggcd((p, 0), (x, 1))

def half(zi, z0, N):
    """primitive A with zi/z0 = A/conj(A), or None if unit class is +-i."""
    G = mul(zi, cj(z0))
    X, Y = N+G[0], G[1]
    if X == 0 and Y == 0: return None
    g = gcd(X, Y); A = (X//g, Y//g)
    if nm(A) % 2 == 0: return None
    return A

def points_on(N):
    f = factorint(N)
    pts = [(1, 0)]
    for p, e in f.items():
        assert p % 4 == 1
        pi = gprime(p); pic = cj(pi)
        new = []
        for z in pts:
            for k in range(e+1):
                t = z
                for _ in range(k): t = mul(t, pi)
                for _ in range(e-k): t = mul(t, pic)
                new.append(t)
        pts = new
    out = set()
    for z in pts:
        for un in [(1,0),(0,1),(-1,0),(0,-1)]:
            out.add(mul(z, un))
    assert all(nm(z) == N for z in out)
    return sorted(out, key=lambda z: math.atan2(z[1], z[0]))

stats = dict(triples=0, skipped_unit=0, ok=0)

def check_triple(z0, z1, z2, N):
    P1 = half(z1, z0, N); P2 = half(z2, z0, N); A12 = half(z2, z1, N)
    if P1 is None or P2 is None or A12 is None:
        stats['skipped_unit'] += 1; return
    stats['triples'] += 1
    d01, d02, d12 = nm(P1), nm(P2), nm(A12)
    a = nm(ggcd(P1, P2)); b = d01//a; c = d02//a
    assert a*b == d01 and a*c == d02
    assert b*c == d12, (b, c, d12)                        # label of third edge
    assert d01*d02 % d12 == 0 and isqrt(d01*d02//d12)**2 == d01*d02//d12 == a*a
    assert gcd(gcd(a, b), c) == 1                          # chain condition, k=3
    assert N % (a*b*c) == 0                                 # profile condition
    s, w = abs(P1[0]), abs(P1[1]); v = abs(P2[1])
    D = det(P1, P2); assert D % a == 0
    u = abs(D)//a
    assert u == abs(A12[1])                                 # t_12 = det/n_12 is the third offset
    assert w > 0 and v > 0 and u > 0
    num_p = a*u*u + b*v*v + 2*s*u*v; num_m = a*u*u + b*v*v - 2*s*u*v
    assert c*w*w in (num_p, num_m)                          # (II.1)
    eps = 1 if c*w*w == num_p else -1
    # reverse parametrisation: coordinates recovered
    assert (a*u + eps*s*v) % w == 0 and (s*u + eps*b*v) % w == 0
    X2 = (a*u + eps*s*v)//w; ReQ = (s*u + eps*b*v)//w
    assert X2*X2 == a*c - v*v and abs(X2) == abs(P2[0])
    assert ReQ*ReQ == b*c - u*u and abs(ReQ) == abs(A12[0])
    # Heron
    Aq, Bq, Cq, U = a*u*u, b*v*v, c*w*w, (u*v*w)**2
    assert (Aq+Bq-Cq)**2 == 4*Aq*Bq - 4*U
    stats['ok'] += 1
    return (a, b, c, u, v, w, s, eps)

def window(pts, N, k, C):
    R = math.isqrt(N)  # R^2 = N radius sqrt(N)
    rad = math.sqrt(N)
    best = None
    L = len(pts)
    for i in range(L):
        sub = [pts[(i+j) % L] for j in range(k)]
        th0 = math.atan2(sub[0][1], sub[0][0]); th1 = math.atan2(sub[-1][1], sub[-1][0])
        arc = ((th1-th0) % (2*math.pi))*rad
        if arc <= C*math.sqrt(rad):
            best = sub; break
    return best

def run_cluster(N, k, C, label):
    pts = points_on(N)
    cl = window(pts, N, k, C)
    assert cl is not None, label
    before = stats['ok']
    examples = []
    for i, j, l in combinations(range(k), 3):
        for (x, y, z) in [(i, j, l), (j, l, i), (l, i, j)]:
            r = check_triple(cl[x], cl[y], cl[z], N)
            if r and len(examples) < 2: examples.append(r)
    print(f'{label}: N={N}, {len(pts)} points, cluster of {k} found; triples checked ok: {stats["ok"]-before}; sample (a,b,c,u,v,w,s,eps)={examples}')

if __name__ == '__main__':
    run_cluster(32988780417125, 8, 5.674, 'record 8-point cluster (families.md)')
    run_cluster(235370525, 7, 5.2700, 'record 7-point cluster (families.md)')
    # Part I parent N=131490079981: build the four-point profile 0,v1,v2,v1+v2
    N = 131490079981
    pts = points_on(N)
    # all triples among points in a short arc window of 3 (and full random sample)
    random.seed(1)
    before = stats['ok']
    tri = None
    for _ in range(20000):
        z0, z1, z2 = random.sample(pts, 3)
        r = check_triple(z0, z1, z2, N)
        if r and sorted(r[:3]) == sorted((1417, 20381, 4553)): tri = r
    print(f'N=131490079981 random triples ok: {stats["ok"]-before}; found Astra triangle as genuine triple: {tri}')
    # random circles, all triples on random N
    for N in [5*13*17*29*37, 5**3*13**2*17, 5**2*13*29*41*53, 65*65*85]:
        pts = points_on(N); before = stats['ok']
        for _ in range(5000):
            check_triple(*random.sample(pts, 3), N)
        print(f'random N={N}: ok {stats["ok"]-before}')
    print(stats)
