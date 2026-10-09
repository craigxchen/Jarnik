"""Referee check R2: Lemma 1.2 of construct.md ("(R) holds for every nonzero zero-sum c and every unit u
on actual normalised clusters") fails when v(c) = sum c_x a_x = 0, i.e. on multiplicative rectangles.

Exact integer arithmetic.  Example: N = 65 = 5*13, points 1+8i, 8+i, 4+7i, 7+4i (gcd 1, first quadrant),
(1+8i)(8+i) = 65 i = (4+7i)(7+4i).  Then c = (1,1,-1,-1) has v(c) = 0 and Z_+ = Z_-, so
prod rho_x^{c_x} = 1 modulo every odd q^a: D_1(c) = prod of all odd prime powers <= Q0 coprime to N,
while |Z_+ - Z_-| = 0.  (R) would require |sin(0 - psi_1)| = 0 >= D_1/sqrt2.  The construct's own
check (check_axiomR.py) skips v = 0 ('if not any(v): continue').

The second part re-verifies (R) in the integer form |Z_+ - u Z_-|^2 Norm(A_v) >= 2 D_u^2 N^(n/2) for every
character with v != 0 of this 4-point set and of all 4-subsets of an 8-point arc on N = 5*13*17*29,
including the multi-unit cases, and records the arc constants."""
import itertools, math

def gmul(a, b): return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def gconj(a): return (a[0], -a[1])
def gnorm(a): return a[0]*a[0] + a[1]*a[1]
def gpow(a, k):
    r = (1, 0)
    for _ in range(k): r = gmul(r, a)
    return r
def gdiv_exact(a, b):
    n = gnorm(b); w = gmul(a, gconj(b))
    if w[0] % n or w[1] % n: return None
    return (w[0]//n, w[1]//n)

def odd_prime_powers(Q0, N):
    out = []
    for q in range(3, Q0 + 1, 2):
        if all(q % d for d in range(2, int(q**0.5) + 1)) and N % q:
            m = q
            while m <= Q0:
                out.append((q, m)); m *= q
    return out

UNITS = [(1, 0), (0, 1), (-1, 0), (0, -1)]

def D_u(Zp, Zm, u, pps):
    D = 1; best = {}
    for q, m in pps:
        d = (Zp[0] - gmul(u, Zm)[0], Zp[1] - gmul(u, Zm)[1])
        if d[0] % m == 0 and d[1] % m == 0:
            best[q] = max(best.get(q, 1), m)
    for q, m in best.items(): D *= m
    return D

def exponents(z, pis):
    e = []
    for pi in pis:
        k = 0; w = z
        while True:
            w2 = gdiv_exact(w, pi)
            if w2 is None: break
            w = w2; k += 1
        e.append(k)
    return e

def arc_constant(pts, N):
    angs = sorted(math.atan2(y, x) % (2*math.pi) for x, y in pts)
    # smallest arc containing all points
    gaps = [(angs[(i+1) % len(angs)] - angs[i]) % (2*math.pi) for i in range(len(angs))]
    span = 2*math.pi - max(gaps)
    R = math.sqrt(N)
    return R*span/math.sqrt(R)

def check_R(pts, N, primes, pis, Q0, report_v0=True):
    pps = odd_prime_powers(Q0, N)
    A = [exponents(z, pis) for z in pts]
    M = len(pts)
    viol_v0 = []; ok = 0
    for c in itertools.product(range(-2, 3), repeat=M):
        if sum(c) != 0 or not any(c): continue
        n = sum(abs(t) for t in c)
        Zp = (1, 0); Zm = (1, 0)
        for z, t in zip(pts, c):
            if t > 0: Zp = gmul(Zp, gpow(z, t))
            elif t < 0: Zm = gmul(Zm, gpow(z, -t))
        v = [sum(t*a[j] for t, a in zip(c, A)) for j in range(len(primes))]
        normA = 1
        for p, vj in zip(primes, v): normA *= p**abs(vj)
        for u in UNITS:
            D = D_u(Zp, Zm, u, pps)
            uZm = gmul(u, Zm)
            d = (Zp[0] - uZm[0], Zp[1] - uZm[1])
            lhs = gnorm(d) * normA; rhs = 2 * D * D * N**(n // 2)
            if not any(v):
                if lhs < rhs: viol_v0.append((c, u, D, gnorm(d)))
            else:
                assert lhs >= rhs, (c, u)
                ok += 1
    return ok, viol_v0

# ---- the N = 65 rectangle
N = 65; primes = [5, 13]; pis = [(2, 1), (3, 2)]
pts = [(1, 8), (8, 1), (4, 7), (7, 4)]
assert all(x*x + y*y == N for x, y in pts)
assert gmul(pts[0], pts[1]) == gmul(pts[2], pts[3]) == (0, 65)
A = [exponents(z, pis) for z in pts]
print("N=65 points", pts, "exponent rows", A, " arc constant C = %.4f (>= 2 sqrt 2 = 2.8284)" % arc_constant(pts, N))
for Q0 in (8, 60):
    ok, bad = check_R(pts, N, primes, pis, Q0)
    print("  Q0=%d: (R) holds on all %d (character,unit) instances with v(c) != 0;  violations with v(c)=0: %d"
          % (Q0, ok, len(bad)))
    for c, u, D, nd in bad[:3]:
        print("     e.g. c=%s u=%s: Z+ - u Z- has norm %d, D_u(c)=%d, (R) needs 0 >= %.2f" % (c, u, nd, D, D / 2**0.5))
    assert len(bad) > 0

# ---- 8-point arc on N = 5*13*17*29 (all 4-subsets), Q0 = 60
N2 = 5*13*17*29; primes2 = [5, 13, 17, 29]; pis2 = [(2, 1), (3, 2), (4, 1), (5, 2)]
allpts = []
for x in range(-int(N2**0.5) - 1, int(N2**0.5) + 2):
    y2 = N2 - x*x
    if y2 < 0: continue
    y = int(math.isqrt(y2))
    if y*y == y2:
        for yy in {y, -y}: allpts.append((x, yy))
prim = [z for z in allpts if all(exponents(z, [pi])[0] == 0 or exponents(z, [gconj(pi)])[0] == 0 for pi in pis2)]
prim.sort(key=lambda z: math.atan2(z[1], z[0]))
win = prim[:8]
print("N=%d: %d primitive points; testing the first 8 by angle, arc constant %.2f" % (N2, len(prim), arc_constant(win, N2)))
tot_ok = 0; tot_bad = 0
for sub in itertools.combinations(win, 4):
    ok, bad = check_R(list(sub), N2, primes2, pis2, 60)
    tot_ok += ok; tot_bad += len(bad)
print("  all 70 four-subsets, entries in [-2,2]: (R) holds on %d instances with v != 0; v=0 violations: %d"
      % (tot_ok, tot_bad))
print("R2 DONE")
