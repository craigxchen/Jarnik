"""Exact checks for the abc section of conditional.md (Section 8).

(E) Pair identity: for lattice points z != w on x^2+y^2=N with |z-w| small,
    g = gcd(z,w), z = g u, w = g eta conj(u) (eta a unit), u conjugate-primitive,
    and with (s,t) as in Prop. 6.1: s^2+t^2 in {n, 2n}, n = N/|g|^2,
    |u - eta conj(u)|^2 = 4s^2 (eta=+-1) or 2 s^2 (eta=+-i).
    Report log rad(n)/log R (abc predicts >= 1/2 - o(1) for C sqrt(R) pairs)
    and the abc quality q = log c0 / log rad(s^2 t^2 c0).
(F) Ptolemy relation for 4 cyclically ordered points: the three products are
    integer multiples a,b,c of one primitive omega with a+b=c; abc quality.
"""
import sys, random
from math import gcd, log, isqrt, atan2
from itertools import combinations
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from gauss import mul, conj, norm, sub, ggcd, divexact, UNITS, rad, two_squares

def pair_data(z, w):
    g = ggcd(z, w)
    u = divexact(z, g); up = divexact(w, g)
    n = norm(u); assert n == norm(up)
    eta = None
    for e in UNITS:
        if mul(e, conj(u)) == up:
            eta = e
    assert eta is not None, "w/g is not a unit times conj(z/g)"
    A, B = u
    assert gcd(A, B) == 1
    if eta in [(1, 0), (-1, 0)]:
        s, t = (B, A) if eta == (1, 0) else (A, B)
        c0 = n; dd = norm(sub(u, mul(eta, conj(u))))
        assert dd == 4*s*s
    else:
        s, t = (A - B, A + B) if eta == (0, 1) else (A + B, A - B)
        c0 = 2*n; dd = norm(sub(u, mul(eta, conj(u))))
        assert dd == 2*s*s
    assert s*s + t*t == c0 and s != 0 and gcd(s, t) == 1
    return g, n, s, t, c0

def angle_sorted(pts):
    return sorted(pts, key=lambda z: atan2(z[1], z[0]))

def close_pairs(pts, N, C):
    R2 = N  # R^2
    out = []
    P = angle_sorted(pts)
    m = len(P)
    for i in range(m):
        for k in range(1, m):
            j = (i + k) % m
            d2 = norm(sub(P[i], P[j]))
            if d2 * d2 > C**4 * R2:   # |z-w|^2 > C^2 R  <=>  d2^2 > C^4 N
                break
            if i < j or k < m - i:
                out.append((P[i], P[j]))
    return out

SPLIT = [5, 13, 17, 29, 37, 41, 53, 61, 73, 89, 97, 101, 109, 113, 137, 149, 157, 173, 181, 193, 197]

def run_E(C=4, trials=60, seed=7, Nmax=3*10**10):
    rng = random.Random(seed)
    worst = None; maxq = 0; cnt = 0; circles = 0
    for _ in range(trials):
        N = 1
        while True:
            p = rng.choice(SPLIT); e = rng.choice((1, 1, 1, 2))
            if N * p**e > Nmax: break
            N *= p**e
        pts = two_squares(N); circles += 1
        for z, w in close_pairs(pts, N, C):
            g, n, s, t, c0 = pair_data(z, w)
            logR = 0.5*log(N)
            r = rad(n)
            ex = log(r)/logR if r > 1 else 0.0
            q = log(c0)/log(rad(s*s*t*t*c0))
            cnt += 1
            if worst is None or ex < worst[0]:
                worst = (ex, N, z, w, n, r, s, t)
            maxq = max(maxq, q)
    print(f"(E) {circles} random circles N<={Nmax:.0e}, C={C}: {cnt} pairs with |z-w|<=C sqrt(R)")
    print(f"    min log rad(n)/log R = {worst[0]:.3f}  at N={worst[1]}, n={worst[4]}, rad={worst[5]}, s={worst[6]}")
    print(f"    max abc quality over the triples s^2+t^2=c0: {maxq:.3f}")

def run_powerful(C=4):
    print("(E') powerful circles N = 5^a 13^b (rad N = 65): pairs at distance <= C sqrt(R)")
    for a in range(0, 13):
        for b in range(0, 9):
            N = 5**a * 13**b
            if N < 10**4 or N > 10**12: continue
            pts = two_squares(N)
            prs = close_pairs(pts, N, C)
            # remove pairs that are unit rotations / conjugate-trivial? keep all; report count
            if prs:
                ex = [log(rad(pair_data(z, w)[1]) or 1)/(0.5*log(N)) for z, w in prs]
                print(f"    N=5^{a} 13^{b}: {len(prs)} close pairs, max log rad(n)/log R = {max(ex):.3f}")
    print("    (no line printed for a given (a,b) means no close pair)")

def ptolemy(z1, z2, z3, z4):
    P1 = mul(sub(z1, z3), sub(z2, z4))
    P2 = mul(sub(z1, z2), sub(z3, z4))
    P3 = mul(sub(z1, z4), sub(z2, z3))
    assert P1 == (P2[0] + P3[0], P2[1] + P3[1])
    # primitive omega on the common line
    g = 0
    for P in (P1, P2, P3):
        g = gcd(g, gcd(P[0], P[1]))
    # line direction: use P1
    d = gcd(P1[0], P1[1]); om = (P1[0]//d, P1[1]//d)
    def coeff(P):
        if om[0] != 0:
            assert P[0] % om[0] == 0; k = P[0]//om[0]
        else:
            k = P[1]//om[1]
        assert (k*om[0], k*om[1]) == P, "not collinear with omega"
        return k
    c, a, b = coeff(P1), coeff(P2), coeff(P3)
    assert a + b == c
    return a, b, c

def run_F(avals=(20, 50, 100, 300)):
    from itertools import product
    print("(F) Ptolemy relation a+b=c on 4-point subsets of the six-point CG family")
    sig = [s for s in product((1, -1), repeat=4) if sum(s) == 0]
    for a0 in avals:
        pts = []
        for s in sig:
            z = (1, 0)
            for j, sj in enumerate(s, start=1):
                z = mul(z, (a0 + j, sj))
            pts.append(z)
        P = angle_sorted(pts)
        N = norm(P[0])
        qmax = 0; ratios = []
        for S in combinations(range(6), 4):
            a, b, c = ptolemy(*[P[i] for i in S])
            g = gcd(gcd(a, b), c); a, b, c = a//g, b//g, c//g
            if 0 in (a, b, c): continue
            r = rad(abs(a*b*c))
            q = log(max(abs(a), abs(b), abs(c)))/log(r)
            qmax = max(qmax, q)
            ratios.append(log(r)/log(max(abs(a), abs(b), abs(c))))
        print(f"    a={a0}: N={N:.3e}, 15 quadruples: max abc quality {qmax:.3f}, "
              f"min log rad(abc)/log max = {min(ratios):.3f}")

if __name__ == '__main__':
    run_E(); run_powerful(); run_F()
