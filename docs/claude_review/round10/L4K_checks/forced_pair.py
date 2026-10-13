# For stored good systems (forced_coh pickles), max over pairs (zeta_j, zeta_k) of norm <= H of the log-mass share
# of forced moduli (Q_7 part, and all of F) on which columns (j,k) are pair-coherent (one unit u_q).
import sys, math, pickle
import numpy as np
from lib import *

def pair_share(system, P, F, j, k, H, only7):
    Zs = [(a, b) for a in range(1, math.isqrt(H) + 1) for b in range(0, math.isqrt(H) + 1) if a * a + b * b <= H]
    qs = [q for q in F if (q % 8 == 7 or not only7)]
    tot = sum(math.log(q) for q in qs)
    mass = np.zeros((len(Zs), len(Zs)), dtype=np.float32)
    for q in qs:
        T = Torus(q); n = T.m // 4
        La = np.array([T.L(z) % n if (z[0] ** 2 + z[1] ** 2) % q else n for z in Zs])
        tab = np.zeros((n + 1, n + 1), dtype=bool)
        cj, ck = system[q][j] % n, system[q][k] % n
        for u in range(1, n + 1):
            if math.gcd(u, n) == 1: tab[(u * cj) % n, (u * ck) % n] = True
        tab[n, :] = True; tab[:, n] = True     # zeta not a unit mod q: counted as coherent (conservative)
        mass += math.log(q) * tab[La[:, None], La[None, :]]
    i1, i2 = np.unravel_index(np.argmax(mass), mass.shape)
    return float(mass[i1, i2]) / tot, Zs[i1], Zs[i2], len(qs)

if __name__ == '__main__':
    f = sys.argv[1]; H = int(sys.argv[2])
    d = pickle.load(open(f, 'rb'))
    P, F = d['P'], d['F']; r = len(P)
    for ex in d['examples']:
        for (j, k) in [(0, 1), (r - 2, r - 1)]:
            for only7 in (True, False):
                sh, z1, z2, nq = pair_share(ex, P, F, j, k, H, only7)
                print(f"{f} H={H} cols ({j},{k}) {'Q_7 part' if only7 else 'all of F'} ({nq} moduli): max pair-coherent share {sh:.3f} at {z1},{z2}", flush=True)
