"""End-to-end exact check of the factored residue design of round6/sharp.md §3 (Lemma 3.2):
realisability through the exponent matrix with the forced parity, prime powers, level
compatibility, units (k = 0), and the resulting collision pattern.

Profile: Paley q0 = 43 (M = 44), b = 33 copies, canonical assignment (profile.py).
Primes: r = 33*43 actual primes p_j = 1 mod 4 above 10^7 (only their Legendre symbols matter).
Moduli: every odd prime power Q = q^a <= X = 2 M^A (A = 2 here).
Design (Def. 3.1):
  small q < 4M:  level-1 class  sigma_q(x) + 2 (x mod p'(q))  in Z/m_q,  CRT with x mod q^(a-1)
                 for a <= a*(q) = min{a : p'(q) q^(a-1) >= M}; random lifts above a*;
  large q > 4M:  level-1 class  sigma_q(x) + 2 y_x,  y a uniformly random injection into Z/(m_q/2);
                 random lifts at higher levels.
  sigma_q(x) = (A nu(q))_x mod 2,  nu_j(q) = [p_j is a non-residue mod q]  (forced parity).
Realisation: ell = nu + 2 ell',  ell' = sum_x y'_x (-H(x,alpha x)) (e_f(x) - e_s(x))  (twins),
then r_j with r_j/conj(r_j) = g^ell_j and Norm r_j = p_j mod Q (Hilbert 90 + square root).
Checked exactly, at every level a of every q:
  (1) Norm rho_x = N mod q^a for rho_x = prod_j r_j^(a_xj) conj(r_j)^(1 - a_xj);
  (2) rho_x/rho_y = eta mod q^a (eta in {1,i,-1,-i})  <=>  kappa_x - kappa_y = iota*s mod m_(q^a)
      (iota = log_g i), for all pairs;
  (3) the line-1 pair collisions are exactly those predicted: small q, a <= a*:
      sigma equal, p'(q) | x - y, q^(a-1) | x - y;  never otherwise;
  (4) for 300 random characters c per modulus: prod rho^c = eta  <=>  sum c kappa = iota*s;
  (5) a parity-blind design (sigma ignored) is NOT realisable at the moduli where sigma_q is
      nonconstant (no r_j of norm p_j exists for some j): counted.
Prints ALL REALISATION CHECKS PASSED."""
import random, math, sys
import numpy as np
from gres import *
from profile import build
from check_design import assignment

A_EXP = float(sys.argv[1]) if len(sys.argv) > 1 else 2.0
q0, b = 43, 33
H, S, Amat, flip, sib, alpha = build(q0, b)
M, r = Amat.shape
X = int(2 * M ** A_EXP)
rng = random.Random(2026)
pp = assignment(int(math.log2(4 * M)) + 2)

# actual primes = 1 mod 4 above 10^7
plist = []
n = 10 ** 7 + 1
while len(plist) < r:
    if n % 4 == 1 and is_prime(n):
        plist.append(n)
    n += 4
assert min(plist) > X

ok = True
stats = dict(moduli=0, levels=0, nonconst_sigma=0, blind_unrealisable=0, pairs_line1=0)
qs = [q for q in primes_upto(X) if q > 2]
for q in qs:
    Atop = 1
    while q ** (Atop + 1) <= X:
        Atop += 1
    Q = q ** Atop
    mtop = m_of(q, Atop)
    nu = np.array([1 if legendre(p, q) == -1 else 0 for p in plist], dtype=np.int64)
    sigma = (Amat @ nu) % 2
    stats["nonconst_sigma"] += int(sigma.min() != sigma.max())
    small = q < 4 * M
    # ---- design: kappa at every level, as integers mod m_(q^a)
    m1 = m_of(q, 1)
    if small:
        pq = pp[q]
        k1 = [(int(sigma[x]) + 2 * (x % pq)) % m1 for x in range(M)]
        astar = 1
        while pq * q ** (astar - 1) < M:
            astar += 1
    else:
        ys = rng.sample(range(m1 // 2), M)
        k1 = [(int(sigma[x]) + 2 * ys[x]) % m1 for x in range(M)]
        astar = 1
    kap = {1: k1}
    for a in range(2, Atop + 1):
        ma, mprev = m_of(q, a), m_of(q, a - 1)
        if small and a <= astar:
            qa1 = q ** (a - 1)
            # CRT of (k1 mod m1, x mod q^(a-1)) in Z/(m1 q^(a-1))
            kap[a] = []
            for x in range(M):
                t = (x - k1[x]) * pow(m1, -1, qa1) % qa1
                kap[a].append((k1[x] + m1 * t) % ma)
        else:
            kap[a] = [(kap[a - 1][x] + mprev * rng.randrange(q)) % ma for x in range(M)]
        ok &= all(kap[a][x] % mprev == kap[a - 1][x] for x in range(M))
    ktop = kap[Atop]
    # ---- realisation through the exponent matrix
    aNu = (Amat @ nu)
    yprime = []
    for x in range(M):
        diff = (ktop[x] - int(aNu[x])) % mtop
        ok &= diff % 2 == 0
        yprime.append(diff // 2)
    half = mtop // 2
    ellp = np.zeros(r, dtype=object)
    for x in range(M):
        h = int(H[x, alpha[x]])
        ellp[flip[x]] = (ellp[flip[x]] + yprime[x] * (-h)) % half
        ellp[sib[x]] = (ellp[sib[x]] - yprime[x] * (-h)) % half
    ell = [(int(nu[j]) + 2 * int(ellp[j])) % mtop for j in range(r)]
    # check A ell = ktop (mod mtop)
    for x in range(M):
        ok &= sum(int(Amat[x, j]) * ell[j] for j in range(r)) % mtop == ktop[x]
    g = generator_T(q, Atop, seed=q)
    rj = []
    for j in range(r):
        t = gpow(g, ell[j], Q)
        r0 = hilbert90(t, q, Q, rng)
        n0 = norm(r0, Q)
        s = sqrt_unit(plist[j] * pow(n0, -1, Q) % Q, q, Atop)
        ok &= s is not None
        rr = ((r0[0] * s) % Q, (r0[1] * s) % Q)
        ok &= norm(rr, Q) == plist[j] % Q
        rj.append(rr)
    rho = []
    for x in range(M):
        u = (1, 0)
        for j in range(r):
            u = mul(u, rj[j] if Amat[x, j] else conj(rj[j], Q), Q)
        rho.append(u)
    Nmod = 1
    for p in plist:
        Nmod = Nmod * p % Q
    # ---- checks at every level
    for a in range(1, Atop + 1):
        Qa, ma = q ** a, m_of(q, a)
        ga = (g[0] % Qa, g[1] % Qa)
        iota = None
        t = (1, 0)
        for k in range(ma):
            if t == (0, 1):
                iota = k
                break
            t = mul(t, ga, Qa)
        ok &= iota is not None and (4 * iota) % ma == 0
        units = {(1, 0): 0, (0, 1): 1, (Qa - 1, 0): 2, (0, Qa - 1): 3}
        rh = [(u[0] % Qa, u[1] % Qa) for u in rho]
        ok &= all(norm(u, Qa) == Nmod % Qa for u in rh)
        for x in range(M):
            for y in range(x + 1, M):
                quot = mul(rh[x], inv(rh[y], Qa), Qa)
                dk = (kap[a][x] - kap[a][y]) % ma
                pred = {s for s in range(4) if dk == (iota * s) % ma}
                act = {units[quot]} if quot in units else set()
                ok &= pred == act
                line1 = (quot == (1, 0))
                stats["pairs_line1"] += line1
                if small and a <= astar:
                    exp = (sigma[x] == sigma[y]) and ((x - y) % pp[q] == 0) and ((x - y) % q ** (a - 1) == 0)
                else:
                    exp = False
                ok &= line1 == bool(exp)
        for _ in range(300):
            c = [0] * M
            supp = rng.sample(range(M), rng.choice([3, 4, 5, 6, 8]))
            for x in supp[:-1]:
                c[x] = rng.choice([-2, -1, 1, 2])
            c[supp[-1]] = -sum(c)
            u = (1, 0)
            for x in range(M):
                if c[x] > 0:
                    u = mul(u, gpow(rh[x], c[x], Qa), Qa)
                elif c[x] < 0:
                    u = mul(u, gpow(inv(rh[x], Qa), -c[x], Qa), Qa)
            sk = sum(c[x] * kap[a][x] for x in range(M)) % ma
            pred = {s for s in range(4) if sk == (iota * s) % ma}
            act = {units[u]} if u in units else set()
            ok &= pred == act
        stats["levels"] += 1
    # ---- (5) parity-blind design: kappa' = x mod p'(q) style without sigma
    if sigma.min() != sigma.max():
        blind = [(2 * (x % (m1 // 2))) % m1 for x in range(M)]
        # realisable iff blind = A nu mod 2 (+ const): fails since sigma nonconstant
        diffs = {(blind[x] - int(aNu[x])) % 2 for x in range(M)}
        stats["blind_unrealisable"] += int(len(diffs) == 2)
    stats["moduli"] += 1

print("parameters: q0 =", q0, " M =", M, " b =", b, " r =", r, " X = 2 M^A =", X)
print("stats:", stats)
print("ALL REALISATION CHECKS PASSED" if ok else "FAILURE")
