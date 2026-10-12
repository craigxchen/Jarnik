"""Evidence for Section D (discrete realisability with equal weights).

Equal weights tau (primes in [P, P(1 + 1/r)]): the pair margin is (tau/2) * (-<s_x, s_y>) + O(1),
so a demand D_xy needs the integer correlation <s_x, s_y> <= -k_xy, k_xy = max(1, ceil(2 D_xy/tau)).
We generate D from the random independent residue model (unit-1 collisions at odd prime powers
<= X = 3M) and ask for the least number r of +-1 columns realising k:
  * lower bound (proved, Perron): r >= lam_max(k);
  * uniform-margin cost: (M-1) k_max  (Hadamard copies);
  * simulated annealing (anneal.c): an UPPER bound on the true minimum, both with all columns
    free and with K fixed Sylvester blocks (each giving exactly -1) plus free columns.
Calibration: uniform patterns k = k0 (true minimum (M-1) k0 by Plotkin + Hadamard).
Every reported feasible matrix is re-verified with integer arithmetic (in C and again here).
Evidence only.
"""
import math, subprocess, sys, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ANNEAL = os.path.join(HERE, 'anneal')
rng = np.random.default_rng(7)

def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0:2] = b'\x00\x00'
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]

def demand_random(M, X):
    D = np.zeros((M, M))
    for q in primes_upto(X):
        if q == 2:
            continue
        m = q - 1 if q % 4 == 1 else q + 1
        cls = rng.integers(0, m, size=M)
        Q = q
        while Q <= X:
            D += 2 * math.log(q) * (cls[:, None] == cls[None, :])
            cls = cls * q + rng.integers(0, q, size=M)
            Q *= q
    np.fill_diagonal(D, 0)
    return D

def sylvester(n):
    H = np.array([[1]], dtype=np.int64)
    while H.shape[0] < n:
        H = np.block([[H, H], [H, -H]])
    return H

def run(k, r, K, iters, seed, restarts=3):
    M = k.shape[0]
    inp = f"{M} {r} {K} {iters} {seed} {restarts}\n" + "\n".join(" ".join(str(int(v)) for v in row) for row in k)
    out = subprocess.run([ANNEAL], input=inp, capture_output=True, text=True).stdout.split("\n")
    viol = int(out[0].split()[0])
    if viol == 0:
        S = np.array([[1 if ch == '+' else -1 for ch in line] for line in out[2:2 + M]], dtype=np.int64)
        G = S @ S.T
        if K > 0:
            H = sylvester(M)[:, 1:]
            G = G + K * (H @ H.T)
            G[np.diag_indices(M)] = 0
            Gfull = G
        else:
            Gfull = G.copy(); Gfull[np.diag_indices(M)] = 0
        off = ~np.eye(M, dtype=bool)
        assert np.all(Gfull[off] <= -k[off]), "re-verification failed"
        return True
    return False

def min_r(k, K, start, iters, seed, step=1, limit=4000):
    r = max(1, start)
    while r <= limit:
        if run(k, r, K, iters, seed):
            return r
        r += step
    return None

def main():
    out = []
    M = 32
    # calibration on uniform patterns (truth: (M-1) k0)
    for k0 in (1, 2, 3):
        k = k0 * (np.ones((M, M), dtype=int) - np.eye(M, dtype=int))
        rf = min_r(k, 0, (M - 1) * k0 - 2, 3_000_000, 11, step=2)
        line = f"calibration M={M} uniform k0={k0}: true min r = {(M-1)*k0}, SA found r = {rf}"
        print(line); out.append(line); sys.stdout.flush()
    for M, kappa in ((32, 2.0), (32, 1.0), (64, 2.0)):
        for trial in range(2):
            D = demand_random(M, 3 * M)
            tau = kappa * math.log(M)
            k = np.maximum(1, np.ceil(2 * D / tau - 1e-12)).astype(int)
            np.fill_diagonal(k, 0)
            lam = max(np.linalg.eigvalsh(k.astype(float)))
            kmax = k.max(); rowmax = k.sum(axis=1).max(); avg = k.sum() / (M * (M - 1))
            unif = (M - 1) * kmax
            hdr = (f"M={M} tau={kappa}logM trial {trial}: D_max/logM={D.max()/math.log(M):.2f}, "
                   f"k: mean={avg:.2f} max={kmax}, lam_max(k)={lam:.1f}, max row sum={rowmax}, "
                   f"uniform cost (M-1)k_max={unif}")
            print(hdr); out.append(hdr); sys.stdout.flush()
            iters = 4_000_000 if M == 32 else 8_000_000
            rf = min_r(k, 0, int(lam) - 1, iters, 100 + trial, step=2)
            line = f"   all-free SA: r = {rf}   (ratio to lam_max {rf/lam:.3f}, to uniform {rf/unif:.3f})"
            print(line); out.append(line); sys.stdout.flush()
            best = None
            for K in range(1, kmax + 1):
                re = min_r(k, K, 0, iters, 200 + K, step=2)
                tot = K * (M - 1) + re
                line = f"   Hadamard base K={K}: free r_e = {re}, total r = {tot} (ratio to lam_max {tot/lam:.3f})"
                print(line); out.append(line); sys.stdout.flush()
                best = tot if best is None else min(best, tot)
            line = f"   best found {min(best, rf)} vs lam_max {lam:.1f} vs uniform {unif}"
            print(line); out.append(line); sys.stdout.flush()
    print("ALL FEASIBLE MATRICES RE-VERIFIED EXACTLY")

if __name__ == '__main__':
    main()
