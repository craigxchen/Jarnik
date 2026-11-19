# Independent test of III.1/III.7 local uniqueness and projective reconstruction
# on explicit genuine clusters (golden 8-point family, record clusters, a 10-point window).
import sys, itertools, math
from math import gcd
from fractions import Fraction
import sympy
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref3_proj')
from gauss import *

def all_points(n):
    fac = sympy.factorint(n)
    opts = [[(1, 0)]]
    for p, e in fac.items():
        if p == 2: opts.append([gpow((1, 1), e)])
        elif p % 4 == 3:
            if e % 2: return []
            opts.append([(p**(e//2), 0)])
        else:
            pi = split_prime(p)
            opts.append([gm(gpow(pi, a), gpow(gc(pi), e-a)) for a in range(e+1)])
    pts = set()
    for combo in itertools.product(*opts):
        z = (1, 0)
        for c in combo: z = gm(z, c)
        for u in [(1, 0), (0, 1), (-1, 0), (0, -1)]:
            pts.add(gm(z, u))
    return sorted(pts, key=lambda z: math.atan2(z[1], z[0]))

def best_window(n, k):
    pts = all_points(n)
    ang = [math.atan2(z[1], z[0]) for z in pts]
    L = len(pts); best = None
    for s in range(L):
        e = (s+k-1) % L
        w = (ang[e]-ang[s]) % (2*math.pi)
        if best is None or w < best[0]: best = (w, [pts[(s+j) % L] for j in range(k)])
    R = math.sqrt(n)
    return best[1], best[0]*math.sqrt(R)

def fib(k):
    a, b = 0, 1
    for _ in range(k): a, b = b, a+b
    return a

def golden(n):
    J = [(fib(n+s+1), fib(n+s)) for s in range(4)]
    F = (1, 2)
    pts = []
    for w in itertools.product([0, 1], repeat=4):
        p = sum(w)
        if p % 2 == 0: continue
        z = F if p == 1 else gc(F)
        for s in range(4):
            z = gm(z, J[s] if w[s] else gc(J[s]))
        pts.append(z)
    pts.sort(key=lambda z: math.atan2(z[1], z[0]))
    nn = gn(pts[0])
    ang = [math.atan2(z[1], z[0]) for z in pts]
    return pts, (max(ang)-min(ang))*nn**0.25

def minor_prod(S, P, pairing):
    (p, q), (r, s) = pairing
    return det(P[p], P[q])*det(P[r], P[s])

def kernel_small(A, B, C, Ux, Uy, Uz):
    """All primitive (x,y,z), nonzero, |x|<=Ux,|y|<=Uy,|z|<=Uz with A x - B y + C z = 0,
    found via modulus B: A x + C z = 0 mod B. Exhaustive over the lattice by reduced basis."""
    sols = set()
    if B == 1:
        return None  # no modulus information
    r = (-C*pow(A, -1, B)) % B
    # lattice {(x,z)}: basis (B,0),(r,1). Weighted Gauss reduction with weights 1/Ux,1/Uz
    def nrm(v): return Fraction(v[0], 1)**2/Fraction(Ux)**2 + Fraction(v[1])**2/Fraction(Uz)**2
    def dot(u, v): return Fraction(u[0]*v[0])/Fraction(Ux)**2 + Fraction(u[1]*v[1])/Fraction(Uz)**2
    u, v = (B, 0), (r, 1)
    if nrm(u) < nrm(v): u, v = v, u
    while True:
        # make v shortest
        if nrm(u) < nrm(v): u, v = v, u
        mu = round(dot(u, v)/nrm(v))
        if mu == 0: break
        u = (u[0]-mu*v[0], u[1]-mu*v[1])
        if nrm(u) >= nrm(v): break
    if nrm(u) < nrm(v): u, v = v, u
    # now v short, u reduced; any lattice vector with weighted norm <= 2 has bounded coeffs
    G = [[nrm(v), dot(u, v)], [dot(u, v), nrm(u)]]
    detG = G[0][0]*G[1][1]-G[0][1]**2
    bound = 2  # weighted norm^2 <= 2 inside the box
    c1 = int(math.isqrt(int(bound*G[1][1]/detG))+2)
    c2 = int(math.isqrt(int(bound*G[0][0]/detG))+2)
    for a1 in range(-c1, c1+1):
        for a2 in range(-c2, c2+1):
            x = a1*v[0]+a2*u[0]; z = a1*v[1]+a2*u[1]
            if x == 0 or z == 0 or abs(x) > Ux or abs(z) > Uz: continue
            if (A*x+C*z) % B: continue
            y = (A*x+C*z)//B
            if y == 0 or abs(y) > Uy: continue
            gg = gcd(gcd(x, y), z)
            s = (x//gg, y//gg, z//gg)
            if s[0] < 0: s = (-s[0], -s[1], -s[2])
            sols.add(s)
    return sols

ORACLE=False
def analyse(name, pts, Cconst, verbose=True):
    S = stars(pts)
    m, P, g, d, t, N = S['m'], S['P'], S['g'], S['d'], S['t'], S['N']
    W = math.log(N)
    # arc width of the points (after gcd removal angles unchanged)
    ang = [math.atan2(z[1], z[0]) for z in S['zs']]
    Phi = max(ang)-min(ang)
    def G(i, j): i, j = min(i, j), max(i, j); return d[i, j]
    def tb(i, j):  # a priori bound on |t_ij| from d_ij and the arc width only
        if ORACLE: return abs(t[min(i,j),max(i,j)])  # best-case: exact sizes known
        return 2*math.sqrt(G(i, j))*math.sin(Phi/2)
    sepmin = min(math.log(d[k])/W for k in d)
    Cconst = Phi*N**0.25
    print(f"== {name} (oracle U={ORACLE}): m={m}, C={Cconst:.6f}, N has {len(str(N))} digits, primes={S['fac']}")
    print(f"   min log d_ij / W = {sepmin:.4f} (pair separation needs > 1/2);",
          f"max |t| = {max(abs(v) for v in t.values())}")
    stats = dict(quartets=0, plucker_ok=0, lcm_ok=0, certified=0, rec_ok=0, rec_bad=0, ambiguous=0)
    rec = {}
    for q in itertools.combinations(range(m), 4):
        i, j, k, l = q
        stats['quartets'] += 1
        X1 = det(P[i], P[j])*det(P[k], P[l]); X2 = det(P[i], P[k])*det(P[j], P[l]); X3 = det(P[i], P[l])*det(P[j], P[k])
        if X1-X2+X3 == 0: stats['plucker_ok'] += 1
        g1, g2, g3 = g[i, j]*g[k, l], g[i, k]*g[j, l], g[i, l]*g[j, k]
        c = gcd(gcd(g1, g2), g3); A, B, C = g1//c, g2//c, g3//c
        D1, D2, D3 = d[i, j]*d[k, l], d[i, k]*d[j, l], d[i, l]*d[j, k]
        L = sympy.ilcm(D1, D2, D3)
        Ad, Bd, Cd = sympy.integer_nthroot(L//D1, 2)[0], sympy.integer_nthroot(L//D2, 2)[0], sympy.integer_nthroot(L//D3, 2)[0]
        if (Ad, Bd, Cd) == (A, B, C) and L % D1 == 0: stats['lcm_ok'] += 1
        x = t[i, j]*t[k, l]; y = t[i, k]*t[j, l]; z = t[i, l]*t[j, k]
        assert A*x-B*y+C*z == 0
        Ux = tb(i, j)*tb(k, l); Uy = tb(i, k)*tb(j, l); Uz = tb(i, l)*tb(j, k)
        assert all(abs(v) <= U*(1+1e-12)+1 for v, U in [(x, Ux), (y, Uy), (z, Uz)])
        Ux, Uy, Uz = int(Ux), int(Uy), int(Uz)
        # certification: some modulus exceeds 2*product of other two bounds
        cert = None
        for (M, U1, U2, perm) in [(B, Ux, Uz, 'B'), (A, Uy, Uz, 'A'), (C, Ux, Uy, 'C')]:
            if M > 2*U1*U2: cert = perm; break
        if cert is None: continue
        stats['certified'] += 1
        # recover from d's only, using modulus B if possible else permute equation
        if cert == 'B':
            sols = kernel_small(Ad, Bd, Cd, Ux, Uy, Uz)
        elif cert == 'A':
            # A x = B y - C z  -> modulus A on (y,z): B y - C z = 0 mod A ; reuse by renaming
            s2 = kernel_small(Bd, Ad, -Cd, Uy, Ux, Uz)  # B y - A x + (-C) z = 0
            sols = set()
            for (yy, xx, zz) in s2:
                v = (xx, yy, zz)
                if v[0] < 0: v = (-v[0], -v[1], -v[2])
                sols.add(v)
        else:
            s2 = kernel_small(Bd, Cd, -Ad, Uy, Uz, Ux)  # B y - C z - A x = 0, modulus C
            sols = set()
            for (yy, zz, xx) in s2:
                v = (xx, yy, zz)
                if v[0] < 0: v = (-v[0], -v[1], -v[2])
                sols.add(v)
        gg = gcd(gcd(x, y), z); act = (x//gg, y//gg, z//gg)
        if act[0] < 0: act = tuple(-u for u in act)
        if len(sols) > 1: stats['ambiguous'] += 1
        if sols == {act}:
            stats['rec_ok'] += 1
            rec[q] = (Ad*act[0], Bd*act[1], Cd*act[2])  # recovered raw products up to scale
        else:
            stats['rec_bad'] += 1
            print("   RECOVERY MISMATCH", q, sols, act)
    print("   ", stats)
    # projective reconstruction: anchors, chart, propagation using only certified quartets
    def raw(q4, pair1, pair2):
        """recovered value proportional to [pair1][pair2] in quartet q4 (sorted), sign-correct"""
        a_, b_, c_, d_ = q4
        X = rec[q4]
        def sgn_pair(p): return (1 if p[0] < p[1] else -1), tuple(sorted(p))
        s1, p1 = sgn_pair(pair1); s2, p2 = sgn_pair(pair2)
        mt = {frozenset([(a_, b_), (c_, d_)]): X[0], frozenset([(a_, c_), (b_, d_)]): X[1],
              frozenset([(a_, d_), (b_, c_)]): X[2]}
        return s1*s2*mt[frozenset([p1, p2])]
    def chart_val(anc, x):
        i, j, k = anc
        q4 = tuple(sorted((i, j, k, x)))
        if q4 not in rec: return None
        num = raw(q4, (x, j), (k, i)); den = raw(q4, (x, i), (k, j))
        return Fraction(num, den)
    def chart_act(anc, x):
        i, j, k = anc
        return Fraction(det(P[x], P[j])*det(P[k], P[i]), det(P[x], P[i])*det(P[k], P[j]))
    best = None
    for anc in itertools.permutations(range(m), 3):
        placed = {}
        for x in range(m):
            if x in anc: continue
            v = chart_val(anc, x)
            if v is not None: placed[x] = v
        if best is None or len(placed) > len(best[1]): best = (anc, placed)
    anc, placed = best
    ok = all(placed[x] == chart_act(anc, x) for x in placed)
    # iterative propagation: place x from any certified quartet {x,y,z,w} with y,z,w placed
    coords = {anc[0]: None, anc[1]: Fraction(0), anc[2]: Fraction(1)}
    coords.update(placed)
    # represent placed marks by actual chart vectors to propagate via cross ratio equality
    changed = True
    while changed:
        changed = False
        for x in range(m):
            if x in coords: continue
            for trip in itertools.combinations(sorted(coords), 3):
                q4 = tuple(sorted(trip+(x,)))
                if q4 in rec:
                    # recovered cross ratio for (trip) chart equals actual (checked) -> x determined
                    v = chart_val(trip, x); 
                    if v != chart_act(trip, x): ok = False
                    coords[x] = 'prop'; changed = True; break
    print(f"    after iterative propagation: {len(coords)} of {m} marks determined")
    print(f"    best anchor triple {anc}: {len(placed)} of {m-3} other marks placed directly; all equal to actual: {ok}")
    return stats, len(placed), ok

if __name__ == '__main__':
    out = []
    import sys as _s
    ORACLE = len(_s.argv) > 1
    for n in [4, 7, 10, 13, 16, 19, 22]:
        pts, Cc = golden(n)
        out.append(analyse(f"golden n={n}", pts, Cc))
    for (nn, k) in [(32988780417125, 8), (235370525, 7)]:
        pts, Cc = best_window(nn, k)
        out.append(analyse(f"record R^2={nn} k={k}", pts, Cc))
    nn = 5**2*13**2*17*29*37*41*53*61
    for k in [8, 10]:
        pts, Cc = best_window(nn, k)
        out.append(analyse(f"window R^2={nn} k={k}", pts, Cc))
