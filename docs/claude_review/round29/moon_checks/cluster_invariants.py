"""Idea 17 (machine look at invariants of record clusters).  Exact Gaussian arithmetic.
For each N: all lattice points, the best k-point arc, gcd-normalise, reduced chords c_xy, factor N(c_xy),
and compare collision mass at each prime l with the pigeonhole-forced minimum.
"""
import sys, math, itertools
from sympy import factorint
from collections import defaultdict
sys.path.insert(0, '/home/user/Jarnik/docs/claude_review/growth/families')

def gmul(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def conj(z): return (z[0], -z[1])
def ggcd(a, b):
    while b != (0, 0):
        n = b[0]*b[0]+b[1]*b[1]; t = gmul(a, conj(b))
        q = ((2*t[0]+n)//(2*n), (2*t[1]+n)//(2*n)); r = gmul(b, q)
        a, b = b, (a[0]-r[0], a[1]-r[1])
    return a
def gdiv(a, b):
    n = b[0]*b[0]+b[1]*b[1]; t = gmul(a, conj(b))
    assert t[0] % n == 0 and t[1] % n == 0
    return (t[0]//n, t[1]//n)
def nrm(z): return z[0]*z[0]+z[1]*z[1]
def gprime(p):  # pi with norm p, p = 1 mod 4
    for a in range(1, math.isqrt(p)+1):
        b2 = p - a*a; b = math.isqrt(b2)
        if b*b == b2: return (a, b)

def points(N):
    f = factorint(N); pts = [(1, 0)]
    for p, e in f.items():
        if p == 2:
            for _ in range(e): pts = [gmul(z, (1, 1)) for z in pts]
        elif p % 4 == 3:
            assert e % 2 == 0
            pts = [(z[0]*p**(e//2), z[1]*p**(e//2)) for z in pts]
        else:
            pi = gprime(p); pb = conj(pi); new = []
            for z in pts:
                for a in range(e+1):
                    w = z
                    for _ in range(a): w = gmul(w, pi)
                    for _ in range(e-a): w = gmul(w, pb)
                    new.append(w)
            pts = new
    allp = []
    for z in pts:
        w = z
        for _ in range(4): allp.append(w); w = gmul(w, (0, 1))
    return sorted(set(allp), key=lambda z: math.atan2(z[1], z[0]))

def best_arc(N, k):
    P = points(N); n = len(P); R = math.sqrt(N)
    ang = [math.atan2(z[1], z[0]) for z in P]
    best = None
    for s in range(n):
        e = (s+k-1) % n
        span = (ang[e]-ang[s]) % (2*math.pi)
        if best is None or span < best[0]: best = (span, [P[(s+j) % n] for j in range(k)])
    return best[0]*R/math.sqrt(R), best[1]

def analyse(N, k):
    C, pts = best_arc(N, k)
    g = pts[0]
    for z in pts[1:]: g = ggcd(g, z)
    Z = [gdiv(z, g) for z in pts]; N1 = nrm(Z[0])
    print(f"N={N} k={k} C~{C:.6f} normalised N'={N1} = {factorint(N1)}")
    mass = defaultdict(int); pairs = defaultdict(int)
    for x, y in itertools.combinations(range(k), 2):
        d = (Z[x][0]-Z[y][0], Z[x][1]-Z[y][1]); gg = ggcd(Z[x], Z[y]); c = gdiv(d, gg)
        for l, v in factorint(nrm(c)).items(): mass[int(l)] += int(v); pairs[int(l)] += 1
    tot = sum(int(v)*math.log(int(l)) for l, v in mass.items())
    small = {l: (pairs[l], mass[l]) for l in sorted(mass) if l <= 4*k}
    big = sorted(((pairs[l], l) for l in mass if l > 4*k), reverse=True)[:6]
    print("  sum log N(c) =", round(tot, 2), "; archimedean cap binom(k,2)*log(C^2 R)/... ~", round(math.comb(k,2)*math.log(max(C,1)**2*math.sqrt(N1)), 2))
    print("  small primes l<=4k: (pairs with l|N(c), total exponent):", small)
    print("  largest pair-counts at primes l>4k (count, l):", big)
    # forced pigeonhole minimum at l: residues of points mod l (split: l-1 classes on uv=N; inert: l+1; 2: special)
    forced = {}
    for l in small:
        if l == 2 or N1 % l == 0: continue
        cls = l-1 if l % 4 == 1 else l+1
        q, r = divmod(k, cls); forced[l] = r*math.comb(q+1, 2)+(cls-r)*math.comb(q, 2)
    print("  pigeonhole-forced colliding pairs:", forced)

for N, k in [(537606725, 8), (32988780417125, 8), (348984761607530165, 8), (235370525, 7), (683617125425, 5)]:
    analyse(N, k)
