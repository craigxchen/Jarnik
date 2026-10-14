# Finite instances refuting literal L4(K) (Theorem 2 in the test range).
# System: genuine labels at every q <= q* (all levels), a random collision-free parity-respecting vector at
# every q in (q*, Y].  Goodness: genuine part has per-unit demand = q*-smooth part of D_s(Z_v) <= budget; free part
# has no collision (verified exhaustively at level 1, hence at all levels).  So the system is good up to norm K.
# Coherence: for each pair of columns (j,k) and each pair (zeta_j, zeta_k) of Gaussian integers of norm <= H
# (first-quadrant representatives, so units are quotiented, conjugates kept), the Q_7 log-mass of the free range
# on which (c_j, c_k) = u (L_j, L_k) mod n_q for a unit u.  Pair coherence is necessary for system coherence.
import sys, math, random, time
import numpy as np
from lib import *
from free import verify

def random_free(n, par, K, need_par, rng, tries=50):
    r = len(par)
    V = monos(r, K)
    for _ in range(tries):
        c = []
        ok = True
        for j in range(r):
            Vj = [v for v in V if v[j] != 0 and all(t == 0 for t in v[j + 1:])]
            cand = [x for x in range(n) if (not need_par or x % 2 == par[j])]
            rng.shuffle(cand)
            for x in cand:
                cc = c + [x]
                if all(sum(a * b for a, b in zip(v[:j + 1], cc)) % n for v in Vj):
                    c = cc; break
            else:
                ok = False; break
        if ok:
            return c
    return None

def zetas(H):
    out = []
    for a in range(1, math.isqrt(H) + 1):
        for b in range(0, math.isqrt(H) + 1):
            if a * a + b * b <= H:
                out.append((a, b))
    return out

def run(P, K, qstar, Y, H, seed, pairs):
    rng = random.Random(seed)
    r = len(P)
    pis = [gauss_prime(p) for p in P]
    Z = zetas(H); Nz = len(Z)
    Q7 = [q for q in primerange(qstar + 1, Y + 1) if q % 8 == 7]
    th_all = sum(math.log(q) for q in primerange(3, Y + 1) if q % 8 == 7)
    th_gen = sum(math.log(q) for q in primerange(3, qstar + 1) if q % 8 == 7)
    labels = {}
    for q in Q7:
        n = mq(q) // 4; need = (n % 2 == 0); par = nu(P, q)
        c = random_free(n, par, K, need, rng)
        if c is None:   # near the threshold: exhaustive DFS witness, scaled by a random unit (same collision set)
            from free import exists_free
            c0, _ = exists_free(n, par, K, need)
            us = [u for u in range(1, n) if math.gcd(u, n) == 1]
            u = rng.choice(us); c = [(u * x) % n for x in c0]
        assert c is not None and verify(c, n, K, par, need)
        labels[q] = c
    res = {}
    for (j, k) in pairs:
        mass = np.zeros((Nz, Nz), dtype=np.float32)
        gen_mass = 0.0
        for q in Q7:
            n = mq(q) // 4
            T = Torus(q)
            La = np.array([T.L(z) % n if (z[0] * z[0] + z[1] * z[1]) % q else -1 for z in Z])
            tab = np.zeros((n + 1, n + 1), dtype=bool)   # index n = 'undefined'
            cj, ck = labels[q][j] % n, labels[q][k] % n
            for u in range(1, n + 1):
                if math.gcd(u, n) == 1:
                    tab[(u * cj) % n, (u * ck) % n] = True
            Li = np.where(La < 0, n, La)
            mass += math.log(q) * tab[Li[:, None], Li[None, :]]
            # genuine pair for comparison
            Lg = (T.L(pis[j]) % n, T.L(pis[k]) % n)
            if tab[Lg]: gen_mass += math.log(q)
        i1, i2 = np.unravel_index(np.argmax(mass), mass.shape)
        res[(j, k)] = (float(mass[i1, i2]) / th_all, Z[i1], Z[i2], gen_mass / th_all)
    return th_gen / th_all, res, len(Q7), Nz

if __name__ == '__main__':
    allp = split_primes(5, 200)
    rng0 = random.Random(7)
    sets = {}
    for r in range(4, 9):
        sets[('first', r)] = allp[:r]; sets[('last', r)] = allp[-r:]; sets[('rand', r)] = sorted(rng0.sample(allp, r))
    qstars = {}
    import re
    for line in open('free_out.txt'):
        m = re.match(r"K=(\d) r=(\d) (\w+) P=.*?non-free q = (\[.*?\]);", line)
        if m: qstars[(int(m[1]), int(m[2]), m[3])] = max(eval(m[4]))
    H = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    todo = [(int(a.split(',')[0]), int(a.split(',')[1])) for a in sys.argv[2:]] or [(K, r) for K in (2, 3) for r in (4, 6, 8)]
    for K, r in todo:
        if True:
            for tag in ('first', 'last'):
                P = sets[(tag, r)]
                QK = 8 * K * (2 * r + 1) ** K + 1
                Y = min(QK, 4700)
                t0 = time.time()
                share_gen, res, nq, Nz = run(P, K, qstars[(K, r, tag)], Y, H, 1, [(0, 1), (r - 2, r - 1)])
                for (j, k), (mx, z1, z2, gm) in res.items():
                    print(f"K={K} r={r} {tag} P={P} Y={Y} (Q_K={QK}) H={H} Nz={Nz} #Q7 free={nq}: genuine-part share={share_gen:.4f}; "
                          f"cols ({j},{k}): max pair-coherent free share={mx:.4f} at zeta=({z1},{z2}); genuine pi pair on free part={gm:.4f}; "
                          f"=> best coherent share <= {share_gen + mx:.4f}  [{time.time()-t0:.0f}s]", flush=True)
