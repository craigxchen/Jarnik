"""Independent finite instance in the spirit of Theorem 2: genuine labels at q <= qstar, a random collision-free
parity-respecting level-1 vector at every Q_7 prime in (qstar, Y]. Exact max over pairs (zeta_a, zeta_b),
Norm <= H, of the Q_7 log-mass on which (c_a, c_b) = u (L^zeta_a, L^zeta_b) mod n for some unit u
(non-unit zeta counted as coherent). Also reports the share inside F_K = {n_q <= A_K}.
usage: ref_instance.py K r first|last Y H seed"""
import sys, math, random
import numpy as np
from rlib import *
K, r, which, Y, H, seed = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
sp = split_primes_upto(200)
ps = sp[:r] if which == "first" else sp[-r:]
qstar = {(2,4,'first'):41,(2,4,'last'):31,(2,6,'first'):47,(2,6,'last'):71,(2,8,'first'):79,(2,8,'last'):89,
         (3,4,'first'):103,(3,4,'last'):137}[(K,r,which)]
AK = K * (2 * (K + 1) ** 2 * r + 2 * K + 1)
rng = random.Random(seed)
V = VK(r, K)
Varr = np.array(V)
def random_free(n, nuv):
    even = n % 2 == 0
    for attempt in range(20000):
        c = []
        okall = True
        for j in range(r):
            cand = [x for x in range(n) if (not even) or x % 2 == nuv[j]]
            rng.shuffle(cand)
            found = False
            for x in cand:
                c2 = c + [x]
                # check all v supported on 0..j with v_j != 0
                sub = Varr[(Varr[:, j] != 0) & np.all(Varr[:, j + 1:] == 0, axis=1)][:, :j + 1]
                if np.all((sub @ np.array(c2)) % n != 0):
                    c = c2; found = True; break
            if not found: okall = False; break
        if okall: return c
    return None
# zeta list: first-quadrant reps (a>0,b>=0), norm <= H
zetas = [(a, b) for a in range(1, int(math.isqrt(H)) + 1) for b in range(0, int(math.isqrt(H)) + 1) if a * a + b * b <= H]
pis = gaussian_primes(ps)
Q7 = [q for q in primes_upto(Y) if q % 8 == 7]
theta_all = sum(math.log(q) for q in Q7)
theta_gen = sum(math.log(q) for q in Q7 if q <= qstar)
FK = [q for q in Q7 if (q + 1) // 4 <= AK]
theta_FK = sum(math.log(q) for q in FK)
freeQ = [q for q in Q7 if q > qstar]
pairs = [(0, 1), (2, 3)] if r >= 4 else [(0, 1)]
best = {pr: 0.0 for pr in pairs}; bestF = {pr: 0.0 for pr in pairs}
Z = len(zetas)
acc = {pr: np.zeros((Z, Z), dtype=np.float64) for pr in pairs}
accF = {pr: np.zeros((Z, Z), dtype=np.float64) for pr in pairs}
for q in freeQ:
    n = (q + 1) // 4
    c = random_free(n, nu(ps, q))
    assert c is not None
    assert all(sum(a * b for a, b in zip(v, c)) % n for v in V)
    T = Torus(q)
    units = [u for u in range(1, n) if math.gcd(u, n) == 1] if n > 1 else [0]
    Lz = []; nonunit = np.zeros(Z, bool)
    for i, zt in enumerate(zetas):
        if (zt[0] ** 2 + zt[1] ** 2) % q == 0: nonunit[i] = True; Lz.append(0)
        else: Lz.append(T.ratio_log(zt) % n)
    Lz = np.array(Lz)
    U = np.array(units)
    prod = (Lz[:, None] * U[None, :]) % n  # Z x phi
    for pr in pairs:
        A = (prod == c[pr[0]]).astype(np.float32); A[nonunit, :] = 1
        B = (prod == c[pr[1]]).astype(np.float32); B[nonunit, :] = 1
        coh = (A @ B.T) > 0
        acc[pr] += coh * math.log(q)
        if q in FK: accF[pr] += coh * math.log(q)
print(f"K={K} r={r} {which} ps={ps} q*={qstar} Y={Y} H={H} (#zeta={Z}) seed={seed} A_K={AK}")
print(f"  genuine share theta_7([1,q*])/theta_7([1,Y]) = {theta_gen/theta_all:.4f}")
for pr in pairs:
    mx = acc[pr].max() / theta_all
    print(f"  columns {pr}: max pair-coherent share of free part = {mx:.4f}; upper bound for coherent share = {mx + theta_gen/theta_all:.4f}")
    if FK[-1] > qstar:
        mxF = accF[pr].max() / theta_FK
        print(f"     restricted to F_K={{n_q<=A_K}} (Q_7 mass {theta_FK:.1f}, forced Q_7 share {theta_gen/theta_FK:.4f}): max pair-coherent share of free part = {mxF:.4f}")
