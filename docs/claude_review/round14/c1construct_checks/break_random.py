# Thm B1: the uniform random fibre (sharp Def 6.3; non-twin labels independent, uniform in their parity
# class, independent over q) in regime (1.4): p_j <= M^h0, X = 3 M^A.  Norm-1 per-unit budget:
#   e_j may collide (unit s) only at prime powers with product <= nu_s^{-1} sqrt(p_j);
#   any level-1 collision at a prime q > sqrt(2 p_j) violates it.
# Exact: P(e_j has no collision at any q in (sqrt(2p_j), X]) =
#   prod_{q=+-3 (8)} (1 - 4/m_q) * prod_{q=+-1 (8), nu_j(q)=0} (1 - 8/m_q)
# (uniform on the m/2 elements of the parity class; 4-torsion {0,m/4,m/2,3m/4}: at q=+-3 mod 8 two of each
# parity, at q=+-1 mod 8 all four even).  Then a seeded simulation with explicit witnesses.
import sys, math, random
sys.path.insert(0, '.')
from lib14 import mq, split_primes
from sympy import primerange
import numpy as np

def run(M, h0, A, ncol, seed):
    X = 3 * M**A
    P = int(M**h0)
    pool = split_primes(P // 2, P)
    rng = random.Random(seed)
    ps = sorted(rng.sample(pool, ncol))
    qs = [q for q in primerange(3, X + 1)]
    qarr = np.array(qs, dtype=np.int64)
    marr = np.array([mq(q) for q in qs], dtype=np.float64)
    cls3 = (qarr % 8 == 3) | (qarr % 8 == 5)
    surv_exact = []
    surv_3only = []
    witnesses = []
    nrng = np.random.default_rng(seed)
    nviol = 0
    for p in ps:
        Y = math.sqrt(2 * p)
        sel = qarr > Y
        sel &= (qarr != p)
        # nu_j(q) = 1 iff p is a non-residue mod q (Euler)
        nuv = np.array([1 if pow(p % int(q), (int(q) - 1) // 2, int(q)) == int(q) - 1 else 0 for q in qarr[sel]], dtype=np.int64)
        qq = qarr[sel]; mm = marr[sel]; c3 = cls3[sel]
        prob = np.where(c3, 4.0 / mm, np.where(nuv == 0, 8.0 / mm, 0.0))
        surv_exact.append(float(np.exp(np.sum(np.log1p(-prob)))))
        surv_3only.append(float(np.exp(np.sum(np.log1p(-np.where(c3, 4.0 / mm, 0.0))))))
        # simulation: label uniform in parity class; hit iff it lands in the torsion part of the class
        u = nrng.random(len(qq))
        hit = u < prob
        if hit.any():
            nviol += 1
            idx = np.nonzero(hit)[0]
            k = idx[-1]                      # largest colliding prime
            q = int(qq[k]); m = int(mq(q))
            # pick the colliding label explicitly: torsion element of the right parity
            par = int(nuv[k]) if not c3[k] else None
            tors = [0, m // 4, m // 2, 3 * m // 4]
            if c3[k]:
                cand = [t for t in tors if t % 2 == (1 if pow(p % q, (q - 1) // 2, q) == q - 1 else 0)]
            else:
                cand = [t for t in tors if t % 2 == 0]
            lab = cand[int(u[k] * 1e9) % len(cand)]
            s = lab // (m // 4) if lab else 0
            # unit s of the collision depends on iota = log i; budget for odd s larger by sqrt2 -> use the max budget
            budget = 0.5 * math.log(2 * p)
            witnesses.append((p, q, lab, m, len(idx), math.log(q) - budget))
    return ps, surv_exact, surv_3only, nviol, witnesses, X

out = []
for (M, h0, A, ncol, seed) in [(30, 3, 3, 200, 1), (30, 3, 4, 200, 2), (50, 2.5, 3, 200, 3), (50, 2.5, 3.5, 200, 4)]:
    ps, se, s3, nv, wit, X = run(M, h0, A, ncol, seed)
    pred = (math.log(math.sqrt(2 * max(ps))) / math.log(X))
    line = (f"M={M} h0={h0} A={A}: X={X}, {ncol} columns, p_j in [{min(ps)},{max(ps)}]\n"
            f"   exact P(column survives norm-1 per-unit): mean {np.mean(se):.4f}, max {max(se):.4f}"
            f"  (q=+-3 mod 8 only: mean {np.mean(s3):.4f}; Mertens prediction (log sqrt(2p)/log X)^2 = {pred**2:.4f})\n"
            f"   expected violating columns {ncol*(1-np.mean(se)):.1f}; simulated (seeded) violating columns {nv} of {ncol}\n"
            f"   first 3 witnesses (p_j, q, label ell_j(q), m_q, #hits, log q - log sqrt(2p_j)):")
    print(line); out.append(line)
    for w in wit[:3]:
        s = f"      {w[0]}, {w[1]}, {w[2]}, {w[3]}, {w[4]}, excess {w[5]:.2f}"
        print(s); out.append(s)
    ex = [w[5] for w in wit]
    s = f"   excess over budget among violators: median {np.median(ex):.2f}, max {max(ex):.2f}; (A - h0/2) log M = {(A - h0/2)*math.log(M):.2f}"
    print(s); out.append(s)
