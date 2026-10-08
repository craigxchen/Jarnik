# Translated Pell eight-point family (primitive_eight_point_translated_pell_family.md):
# lambda=9+4 sqrt5, x_r+y_r sqrt5 = lambda^r, H_r=(x_r+y_r)+2i y_r, F=1+2i,
# D=H_{n-2},A=H_n,B=H_{n+2},C=H_{n+4}; words with odd number p of unbarred factors,
# prefactor F (p=1) or conj F (p=3).
import itertools, math, sys
from gpoly import gmul, ggcd, gdiv, cluster_C

def pell(r):
    x, y = 1, 0
    for _ in range(r):
        x, y = 9*x + 20*y, 4*x + 9*y
    return x, y

def H(r):
    x, y = pell(r)
    return (x + y, 2*y)

def conj(z): return (z[0], -z[1])

def family(n):
    blocks = [H(n-2), H(n), H(n+2), H(n+4)]
    F = (1, 2)
    pts = []
    for signs in itertools.product([0, 1], repeat=4):
        p = sum(signs)
        if p % 2 == 0: continue
        z = F if p == 1 else conj(F)
        for b, s in zip(blocks, signs):
            z = gmul(z, b if s else conj(b))
        pts.append(z)
    return pts

def best_sub(pts):
    # min normalized span over k-subsets of consecutive angles, for each k
    C, N, k, g = cluster_C(pts, True)
    pts2 = [gdiv(p, g) for p in pts]
    p0c = conj(pts2[0])
    angs = []
    for p in pts2:
        q = gmul(p, p0c)
        cands = [q, (-q[1], q[0]), (-q[0], -q[1]), (q[1], -q[0])]
        q = max(cands, key=lambda t: t[0])
        angs.append(math.atan(q[1]/q[0]))
    angs.sort()
    r4 = float(math.isqrt(math.isqrt(N)))
    out = {}
    for kk in range(2, len(angs)+1):
        out[kk] = min(angs[i+kk-1]-angs[i] for i in range(len(angs)-kk+1))*r4
    return out, g, N

for n in [61, 121, 181, 241]:
    pts = family(n)
    out, g, N = best_sub(pts)
    print("n=%d log10 R=%.1f gcd=%s" % (n, math.log10(N)/2, g), " ".join("k%d:%.6f" % (k, v) for k, v in out.items()))
