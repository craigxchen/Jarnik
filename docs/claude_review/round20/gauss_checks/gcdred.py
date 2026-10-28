"""GCD reduction (exact, Python integers). For actual clusters: every k-subset T, g_T = gcd_{Z[i]}(T),
reduced radius R'^2 = N/Norm(g_T), reduced max chord c' = maxchord/|g_T|, exponent alpha = log c'/log R'.
Checks: (1) R' >= sqrt(R) for every subset of >=2 distinct points; (2) average of W_empty/W over ordered k-tuples
>= rho_{M,k} = min_S [C(S,k)+C(M-S,k)]/C(M,k) (normalised cluster, W = log N after removing cluster gcd);
(3) fraction with W_empty >= W/2^k; (4) distribution of alpha.  Also (5) lower bound rho_{M,k} >= 2^{1-k}(1-k(k-1)/(2(M-k))).
EVIDENCE ONLY for (3),(4); (1),(2),(5) are exact checks of proved lemmas."""
import math, sys
from math import comb, isqrt
from itertools import combinations
def gmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gnorm(a): return a[0]*a[0]+a[1]*a[1]
def gdivmod(a,b):
    n = gnorm(b); c = gmul(a,(b[0],-b[1]))
    q = ((2*c[0]+n)//(2*n), (2*c[1]+n)//(2*n))
    r = (a[0]-gmul(q,b)[0], a[1]-gmul(q,b)[1]); return q, r
def ggcd(a,b):
    while b != (0,0):
        _, r = gdivmod(a,b); a, b = b, r
    return a
def gexdiv(a,b):
    q, r = gdivmod(a,b); assert r == (0,0); return q
def points(N):
    P=[]
    for x in range(-isqrt(N), isqrt(N)+1):
        y2 = N-x*x; y = isqrt(y2)
        if y*y == y2:
            P.append((x,y));
            if y: P.append((x,-y))
    P.sort(key=lambda p: math.atan2(p[1],p[0])); return P
def best_window(N, M):
    P = points(N); n = len(P); best=None
    for i in range(n):
        W = [P[(i+j)%n] for j in range(M)]
        d = max(math.dist(a,b) for a in W for b in W)
        if best is None or d < best[0]: best = (d, W)
    return best
def rho(M,k): return min((comb(S,k)+comb(M-S,k))/comb(M,k) for S in range(1,M))
def analyse(N, M, k, label):
    d, Wpts = best_window(N, M)
    g = Wpts[0]
    for z in Wpts[1:]: g = ggcd(g, z)
    Z = [gexdiv(z,g) for z in Wpts]; Nn = N//gnorm(g); R = math.sqrt(Nn); W = math.log(Nn)
    C = d/math.sqrt(R)
    ok1 = True; sumWe = 0.0; cnt = 0; big = 0; alphas = []
    for T in combinations(range(M), k):
        gt = Z[T[0]]
        for i in T[1:]: gt = ggcd(gt, Z[i])
        We = math.log(gnorm(gt)); Rp = math.sqrt(Nn/gnorm(gt))
        if Rp < math.sqrt(R) - 1e-9: ok1 = False
        ch = max(math.dist(Z[a],Z[b]) for a in T for b in T)/math.sqrt(gnorm(gt))
        alphas.append(math.log(ch)/math.log(Rp)); sumWe += We; cnt += 1
        if We >= W/2**k: big += 1
    print(f"{label}: N'={Nn} M={M} C=arc/sqrt(R)~{C:.4f} k={k}: R'>=sqrt(R) all: {ok1}; "
          f"E[W_empty/W]={sumWe/cnt/W:.4f} >= rho={rho(M,k):.4f}: {sumWe/cnt/W >= rho(M,k)-1e-12}; "
          f"frac(W_empty>=W/2^k)={big/cnt:.3f}; alpha min/median/max = {min(alphas):.4f}/{sorted(alphas)[len(alphas)//2]:.4f}/{max(alphas):.4f}")
if __name__ == "__main__":
    for M in (10, 40, 200):
        for k in (3,4,5,6):
            lb = 2**(1-k)*(1-k*(k-1)/(2*(M-k))) if M>2*k else float('nan')
            print(f"rho({M},{k})={rho(M,k):.5f}  bound 2^(1-k)(1-k(k-1)/(2(M-k)))={lb:.5f}  2^(1-k)={2**(1-k):.5f}")
    analyse(1176852625, 10, 5, "referee 10-pt cluster")
    analyse(1176852625, 10, 3, "referee 10-pt cluster")
    analyse(1176852625, 9, 5, "9-pt subcluster")
