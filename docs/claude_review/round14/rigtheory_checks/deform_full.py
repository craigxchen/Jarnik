"""One-modulus deformations with the FULL genuine reference (round 14): every odd prime q != q0 not dividing N carries the
genuine labels at all levels (no free moduli; this is the all-moduli, norm-<=K slice of P3-all* with X >= max|D_s|).
The slack of (v,s) is then log(2N-part of D_s) + log(1/|sin|) only (Lemma 3.1).  Usage: python3 deform_full.py K "ps" q0...|0"""
import sys, math, time
import numpy as np
import deform
from rlib14 import *
from forced import forced_set
def setup_full(ps, K, F, q0, charged_levels=True):
    r = len(ps); pis = gaussian_primes(ps); N = math.prod(ps)
    V = VK_mod_sign(r, K); m, L = genuine_labels(ps, q0); n = m // 4
    alw = np.zeros((len(V), 4), np.uint8); gcoll = np.zeros(len(V), np.uint8); gsig = np.full(len(V), -1, np.int32)
    last = np.zeros(len(V), np.int32)
    for i, v in enumerate(V):
        z = Zv(pis, v); nz = z[0] ** 2 + z[1] ** 2
        for s in range(4):
            D = abs(Ds(z, s))
            # charged part = D with factors 2, p_j (j any) and q0 removed
            c = D
            while c % 2 == 0: c //= 2
            for p in ps:
                while c % p == 0: c //= p
            while c % q0 == 0: c //= q0
            pr = q0 * c
            alw[i, s] = 1 if pr * pr <= nz * (2 if s % 2 else 1) else 0
        t = sum(a * b for a, b in zip(v, L)) % m
        if t % n == 0: gcoll[i] = 1; gsig[i] = t // n
        last[i] = max(j for j in range(r) if v[j] != 0)
    return V, m, L, alw, gcoll, gsig, last
deform.setup = setup_full
if __name__ == "__main__":
    K = int(sys.argv[1]); ps = [int(x) for x in sys.argv[2].split(",")]; q0s = [int(x) for x in sys.argv[3:]]
    F = forced_set(ps, K, 60000)
    if q0s == [0]: q0s = F
    for q0 in q0s:
        t0 = time.time(); d = deform.run_q0(ps, K, F, q0)
        print(f"FULL K={K} q0={q0} n={d['n']} survivors={d['total']} inclass={d['inclass']} strict={d['strict']} "
              f"outside={d['outside']} [F2={d['F2']} F1={d['F1']}] allowed(v,s)={d['n_allowed']} ({time.time()-t0:.1f}s)", flush=True)
