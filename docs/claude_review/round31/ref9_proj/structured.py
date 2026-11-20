"""Windows on structured large circles: all points from the factorisation; certificates and projection chain.
Referee's own code."""
import sys; sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref9_proj')
from cert import *
from math import exp
from sympy import factorint, primerange
import mpmath
mpmath.mp.dps = 40

def all_points(N):
    f = factorint(N)
    pts = [(1, 0)]
    for p, e in f.items():
        if p == 2:
            new = [mul(z, powg((1, 1), e)) for z in pts]
        elif p % 4 == 3:
            assert e % 2 == 0
            new = [(z[0]*p**(e//2), z[1]*p**(e//2)) for z in pts]
        else:
            pi_ = gprime(p)
            new = [mul(z, mul(powg(pi_, a), powg(conj(pi_), e-a))) for z in pts for a in range(e+1)]
        pts = new
    out = set()
    for z in pts:
        for u in UNITS:
            out.add(mul(u, z))
    return sorted(out, key=lambda z: float(mpmath.atan2(z[1], z[0])))

def windows(N, C, MMIN):
    pts = all_points(N); k = len(pts)
    R = mpmath.sqrt(N); lim = C/mpmath.sqrt(R)
    ang = [mpmath.atan2(z[1], z[0]) for z in pts]
    ang2 = ang + [a + 2*mpmath.pi for a in ang]
    res = []; prevj = -1; j = 0
    for i in range(k):
        j = max(j, i)
        while j+1 < i+k and ang2[j+1] - ang2[i] <= lim: j += 1
        if j-i+1 >= MMIN and j > prevj:
            res.append([pts[t % k] for t in range(i, j+1)])
        prevj = max(prevj, j)
    return res

def analyse(N, C, MMIN, small=True, label=""):
    W = windows(N, C, MMIN)
    for cl in W:
        g, zs, Np = normalise(cl)
        M = len(zs)
        primes, pis, e, rows = rows_of(zs, Np)
        a, th, wid = angles(zs)
        best = None; bestpairL1 = None; pairT1 = 0
        cs = []
        for i, j in combinations(range(M), 2):
            c = [0]*M; c[i] += 1; c[j] += 1; cs.append(c)
        for i in range(M):
            c = [0]*M; c[i] = 2; cs.append(c)
        if small and M <= 8:
            for c in product(range(-1, 3), repeat=M):
                if sum(c) == 2 and 3 <= sum(map(abs, c)) <= 4: cs.append(list(c))
        for c in cs:
            T, b, u, v = certificate(zs, Np, primes, pis, e, rows, c)
            l1 = sum(map(abs, c))
            xs = projections(zs, b); Xb = xs[-1]-xs[0]
            assert 2*len(xs) >= M
            q = (1+l1)*sqrt(T)/M**1.5
            if best is None or q < best[0]: best = (q, T, l1, Xb, b)
            if l1 == 2 and max(c) == 1:
                L1 = sum(map(abs, v)); bestpairL1 = L1 if bestpairL1 is None else min(bestpairL1, L1)
                if T == 1: pairT1 += 1
        q, T, l1, Xb, b = best
        print(f"{label} N={N} N'={Np} s={len(primes)} M={M} arcC={float(mpmath.sqrt(N)**0.5*wid):.3f} "
              f"min(1+|c|)sqrtT/M^1.5={q:.3f} at T={T} (log T/log N'={log(T)/log(Np) if T>1 else 0:.3f}) |c|_1={l1} X_b={Xb} "
              f"min pair sum|v|={bestpairL1} pairs with T=1: {pairT1}")
    return len(W)

if __name__ == '__main__':
    for N, C, MM, lab in [(235370525, 5.3, 7, "sporadic7"), (32988780417125, 5.68, 8, "record8"),
                          (537606725, 5.92, 8, "golden8"), (1176852625, 8.24, 10, "record10")]:
        analyse(N, C, MM, label=lab)
    # products of the first k primes = 1 mod 4 (squarefree), generic windows
    P1 = [p for p in primerange(5, 200) if p % 4 == 1]
    for k in (6, 7, 8, 9):
        N = 1
        for p in P1[:k]: N *= p
        for C, MM in ((2, 4), (4, 5), (8, 6)):
            nW = analyse(N, C, MM, small=(k <= 8), label=f"k={k},C={C}")
            print(f"# k={k} C={C}: {nW} windows with M>={MM}")
