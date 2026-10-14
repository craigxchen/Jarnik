"""Packing core (round 14): q is packing-rigid if every parity-respecting level-1 label at q whose every collision
(v,s) is individually affordable (q^2 <= Norm(Z_v) * (2 if s odd)) is genuine up to conjugation (u * (eps o L) mod n).
This uses only the single-collision budget, so it holds for every K-good system and every member of P3-all*
(Prop. 3.7(1)): packing-rigid moduli lie in every joint core, for all norms.  All v with Norm(Z_v) < 2 q^2 are
included (K_q = floor(log(2q^2)/log p_min)); vectors of larger norm are individually affordable and irrelevant.
Usage: python3 packing.py qmax p1,p2,...   (prints per q: #labels avoiding P(q), #outside every eps-class)"""
import sys, math, itertools, ctypes, os
import numpy as np
from rlib14 import *
from deform import U1_orbit_reps
HERE = os.path.dirname(os.path.abspath(__file__))
lib = ctypes.CDLL(os.path.join(HERE, "deformc.so"))
def packing_q(ps, q):
    r = len(ps); pis = gaussian_primes(ps)
    K = int(math.floor(math.log(2 * q * q) / math.log(min(ps))))
    V = [v for v in VK_mod_sign(r, K) if math.prod(p ** abs(x) for p, x in zip(ps, v)) < 2 * q * q]
    m, L = genuine_labels(ps, q); n = m // 4; nuv = nu(ps, q)
    alw = np.zeros((len(V), 4), np.uint8); gcoll = np.zeros(len(V), np.uint8); gsig = np.full(len(V), -1, np.int32)
    last = np.zeros(len(V), np.int32)
    for i, v in enumerate(V):
        nz = math.prod(p ** abs(x) for p, x in zip(ps, v))
        for s in range(4): alw[i, s] = 1 if q * q <= nz * (2 if s % 2 else 1) else 0
        t = sum(a * b for a, b in zip(v, L)) % m
        if t % n == 0: gcoll[i] = 1; gsig[i] = t // n
        last[i] = max(j for j in range(r) if v[j] != 0)
    units = [u for u in range(1, n) if math.gcd(u, n) == 1] if n > 1 else [0]
    epss = [e for e in itertools.product([1, -1], repeat=r) if e[0] == 1]
    cls = sorted({sum(((u * e[j] * L[j]) % n) * n ** (r - 1 - j) for j in range(r)) for u in units for e in epss}) if n > 1 else [0]
    reps, wts = U1_orbit_reps(m, nuv[0])
    Varr = np.ascontiguousarray(np.array(V, np.int32).reshape(-1)); clsarr = np.array(cls, np.int64)
    out = np.zeros(5, np.int64); ex = np.zeros(4 * r, np.int32); nex = ctypes.c_int(0)
    P = lambda a: ctypes.c_void_p(a.ctypes.data)
    alwf = np.ascontiguousarray(alw.reshape(-1)); nuarr = np.array(nuv, np.int32)
    reparr = np.array(reps, np.int32); wtarr = np.array(wts, np.int64)
    lib.run(r, m, len(V), P(Varr), P(alwf), P(last), P(gcoll), P(gsig), P(clsarr), len(cls),
            P(nuarr), len(reps), P(reparr), P(wtarr), P(out), 4, P(ex), ctypes.byref(nex))
    return n, len(V), int(out[0]), int(out[3] + out[4])
if __name__ == "__main__":
    qmax = int(sys.argv[1]); ps = [int(x) for x in sys.argv[2].split(",")]
    N = math.prod(ps); rig = []; first_flex = None
    for q in primes_upto(qmax):
        if q == 2 or N % q == 0: continue
        n, nv, tot, outside = packing_q(ps, q)
        print(f"ps={ps} q={q} n={n} |P(q)|={nv} labels avoiding P(q)={tot} outside every eps-class={outside}", flush=True)
        if outside == 0 and first_flex is None: rig.append(q)
        if outside > 0 and first_flex is None: first_flex = q
    print(f"SUMMARY ps={ps}: packing-rigid initial segment up to {rig[-1] if rig else None}; first non-rigid q={first_flex}")
