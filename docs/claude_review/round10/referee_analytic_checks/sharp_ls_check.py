"""Referee check: the sharp large sieve (constant ~ (2K+1) + O(X)) is NOT the prime cross-gap PCG.

analytic.md Remark 3.4.2 and the Sec. 6 table put "the sharp (X+K) large sieve" in the row "pointwise prime
separation ... yes under GPA (Prop. D.1)".  The sharp large sieve only controls the spacing of the points
4*theta(varpi) at scale 1/X (X = height), i.e. it can only enforce
    WCG_c:  |sin(phi_p - arg varpi)| >= c / max(p, Norm varpi)        (max-scale spacing),
whereas PCG_p needs (p Norm varpi)^(-1/2).

Part 1 (exact integers).  For split p and every primitive COMPOSITE beta (no Gaussian prime is parallel to it)
with p < Norm beta <= 4p:
  (a) LCG_p(arg beta; p) is violated (Lemma 5.3 converse: some alpha, Norm alpha <= p, has p Im(alpha conj beta)^2 < Norm beta);
  (b) PCG_p(arg beta; Y) is violated for Y = 10^5 when some prime lies on a line |Im(w conj beta)| < (Norm beta/p)^(1/2), checked;
  (c) WCG_{1/2} holds against every Gaussian prime varpi (all associates) with Norm varpi <= Y:
      4 Im(varpi conj beta)^2 max(p, Norm varpi)^2 >= Norm varpi Norm beta.
Part 2 (floating point ILLUSTRATION).  Optimal large-sieve constant Lambda(X,K) = largest eigenvalue of the Gram
matrix G_rs = sum_{|k|<=K} e^{i k (psi_r - psi_s)}, psi = 4 * angle, over the elements of G (first-quadrant
Gaussian primes, split p) of norm <= X, for: the genuine angles; a graft putting one angle at distance 1/(40 X) from
another genuine angle; a graft putting it at arg(beta) for a composite beta with p < Norm beta <= 4p.
"""
import math, sys
import numpy as np
from sympy import isprime

def gprime(x, y):
    if x == 0 or y == 0:
        m = abs(x + y); return m % 4 == 3 and isprime(m)
    return isprime(x * x + y * y)

def two_squares(p):
    for b in range(1, math.isqrt(p // 2) + 1):
        a = math.isqrt(p - b * b)
        if a * a + b * b == p: return max(a, b), min(a, b)

Y = 10**5
R = math.isqrt(Y)
primes = [(x, y) for x in range(-R, R + 1) for y in range(-R, R + 1) if 0 < x*x + y*y <= Y and gprime(x, y)]
print("Gaussian primes (all associates) of norm <= %d: %d" % (Y, len(primes)))
P = np.array(primes, dtype=np.int64); NP = (P ** 2).sum(1)

def lcg_violated(beta, n):
    u, v = beta; Nb = u*u + v*v; r = math.isqrt(n)
    for x in range(-r, r + 1):
        for y in range(-r, r + 1):
            Na = x*x + y*y
            if Na == 0 or Na > n: continue
            im = v*x - u*y
            if im != 0 and n * im * im < Nb: return (x, y)
    return None

tot = 0; badA = 0; badC = 0; pcg_viol = 0
for p in (101, 1009, 10009):
    cnt = 0
    for u in range(0, math.isqrt(4 * p) + 1):
        for v in range(0, math.isqrt(4 * p) + 1):
            Nb = u*u + v*v
            if not (p < Nb <= 4 * p) or math.gcd(u, v) != 1: continue
            if gprime(u, v): continue           # composite beta only
            cnt += 1
            if lcg_violated((u, v), p) is None: badA += 1
            im = P[:, 1] * u - P[:, 0] * v       # Im(varpi conj beta)
            mx = np.maximum(NP, p)
            # exact: 4 im^2 mx^2 >= Nvarpi * Nb   (python ints to avoid overflow)
            ok = all(4 * int(a) * int(a) * int(b) * int(b) >= int(c) * Nb for a, b, c in zip(im, mx, NP))
            if not ok: badC += 1
            # PCG_p violated?  prime w (Norm != p) with p*Im^2 < Nb
            if np.any((p * im.astype(object) ** 2 < Nb) & (NP != p)): pcg_viol += 1
    tot += cnt
    print("p=%d: composite primitive beta with p < Norm <= 4p: %d" % (p, cnt))
print("Part 1: %d directions; LCG_p not violated: %d (expect 0); WCG_1/2 failures vs primes of norm <= 1e5: %d (expect 0); "
      "PCG_p violated by a prime of norm <= 1e5: %d of %d" % (tot, badA, badC, pcg_viol, tot))

# Part 2: illustration with the optimal large-sieve constant
X = 1500
G = []
for p in range(5, X + 1, 4):
    if isprime(p):
        a, b = two_squares(p)
        G.append((p, math.atan2(b, a))); G.append((p, math.atan2(a, b)))
psi0 = np.array([4 * t for _, t in G])
def lam(psi, K):
    d = psi[:, None] - psi[None, :]
    with np.errstate(invalid='ignore', divide='ignore'):
        D = np.sin((K + 0.5) * d) / np.sin(d / 2)
    D[np.isclose(np.sin(d / 2), 0.0, atol=1e-15)] = 2 * K + 1
    return np.linalg.eigvalsh(D).max()
# graft target: the element of G above p0 = largest split prime <= X
j0 = len(G) - 2; p0, th0 = G[j0]
# near-coincidence graft: put phi at distance 1/(40 X) from the angle of another prime (index 10)
psi_near = psi0.copy(); psi_near[j0] = psi0[10] + 4.0 / (40 * X); psi_near[j0 + 1] = 4 * (math.pi / 2) - psi_near[j0]
# lattice-direction graft: composite beta with p0 < Norm beta <= 4 p0, nearest in angle to th0
best = None
for u in range(1, math.isqrt(4 * p0) + 1):
    for v in range(0, u + 1):
        Nb = u*u + v*v
        if p0 < Nb <= 4 * p0 and math.gcd(u, v) == 1 and not gprime(u, v):
            d = abs(math.atan2(v, u) - th0)
            if best is None or d < best[0]: best = (d, u, v)
_, ub, vb = best
psi_lat = psi0.copy(); psi_lat[j0] = 4 * math.atan2(vb, ub); psi_lat[j0 + 1] = 4 * (math.pi / 2) - psi_lat[j0]
print("Part 2 (ILLUSTRATION, float): |G| = %d elements of norm <= %d; graft at p0 = %d; lattice graft beta = %d+%di (Norm %d)"
      % (len(G), X, p0, ub, vb, ub*ub + vb*vb))
for K in (X, 4 * X, 16 * X, 64 * X):
    l0, l1, l2 = lam(psi0, K), lam(psi_near, K), lam(psi_lat, K)
    print("  K=%6d: Lambda/(2K+1): genuine %.3f, near-coincidence graft %.3f, lattice-direction graft %.3f"
          % (K, l0 / (2*K+1), l1 / (2*K+1), l2 / (2*K+1)))
