"""One-modulus deformations (round 14).  Usage: python3 deform.py K "p1,p2,..." q0 [q0 ...]
Charged set = exact forced set F_K^ex (genuine labels, all levels); moduli outside it carry collision-free labels
(zero demand, Theorem 1(b)); at q0 every parity-respecting level-1 label is tested (level-1 relaxation at q0).
Prints survivors / in genuine class / strict (genuine signature) / outside class: F2 (collision-subset) / F1 (extra)."""
import sys, ctypes, os, math
import numpy as np
from rlib14 import *
from forced import forced_set
HERE = os.path.dirname(os.path.abspath(__file__))
lib = ctypes.CDLL(os.path.join(HERE, "deformc.so"))
lib.run.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int] + [ctypes.c_void_p] * 5 + [ctypes.c_void_p, ctypes.c_int,
    ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p]
def U1_orbit_reps(m, par):
    U1 = [u for u in range(1, m) if math.gcd(u, m) == 1 and u % 4 == 1]
    seen = set(); reps = []; wts = []
    for x in range(par, m, 2):
        if x in seen: continue
        orb = {(u * x) % m for u in U1}
        seen |= orb; reps.append(x); wts.append(len(orb))
    return reps, wts
def setup(ps, K, F, q0, charged_levels=True):
    r = len(ps); pis = gaussian_primes(ps)
    V = VK_mod_sign(r, K)
    m, L = genuine_labels(ps, q0); n = m // 4
    others = [q for q in F if q != q0]
    alw = np.zeros((len(V), 4), np.uint8); gcoll = np.zeros(len(V), np.uint8); gsig = np.full(len(V), -1, np.int32)
    last = np.zeros(len(V), np.int32)
    for i, v in enumerate(V):
        z = Zv(pis, v); nz = z[0] ** 2 + z[1] ** 2
        for s in range(4):
            D = Ds(z, s); pr = q0
            for q in others: pr *= q ** vq(D, q)
            alw[i, s] = 1 if pr * pr <= nz * (2 if s % 2 else 1) else 0
        t = sum(a * b for a, b in zip(v, L)) % m
        if t % n == 0: gcoll[i] = 1; gsig[i] = t // n
        last[i] = max(j for j in range(r) if v[j] != 0)
    return V, m, L, alw, gcoll, gsig, last
def run_q0(ps, K, F, q0, maxex=20):
    r = len(ps)
    V, m, L, alw, gcoll, gsig, last = setup(ps, K, F, q0)
    n = m // 4; nuv = nu(ps, q0)
    units = [u for u in range(1, n) if math.gcd(u, n) == 1] if n > 1 else [0]
    cls = sorted({sum(((u * L[j]) % n) * n ** (r - 1 - j) for j in range(r)) for u in units}) if n > 1 else [0]
    reps, wts = U1_orbit_reps(m, nuv[0])
    Varr = np.ascontiguousarray(np.array(V, np.int32).reshape(-1))
    clsarr = np.array(cls, np.int64); nuarr = np.array(nuv, np.int32)
    repsarr = np.array(reps, np.int32); wtsarr = np.array(wts, np.int64)
    out = np.zeros(5, np.int64); ex = np.zeros(maxex * r, np.int32); nex = ctypes.c_int(0)
    P = lambda a: a.ctypes.data
    lib.run(r, m, len(V), P(Varr), P(np.ascontiguousarray(alw.reshape(-1))), P(last), P(gcoll), P(gsig), P(clsarr), len(cls),
            P(nuarr), len(reps), P(repsarr), P(wtsarr), P(out), maxex, P(ex), ctypes.byref(nex))
    exs = [list(ex[k * r:(k + 1) * r]) for k in range(nex.value)]
    return dict(q0=q0, m=m, n=n, L=L, total=int(out[0]), inclass=int(out[1]), strict=int(out[2]), F2=int(out[3]), F1=int(out[4]),
                outside=int(out[3] + out[4]), examples=exs, n_allowed=int(alw.sum()), n_gcoll=int(gcoll.sum()))
if __name__ == "__main__":
    K = int(sys.argv[1]); ps = [int(x) for x in sys.argv[2].split(",")]
    q0s = [int(x) for x in sys.argv[3:]]

    import time
    F = forced_set(ps, K, 60000 if K >= 6 else 12000)
    if q0s == [0]: q0s = F
    for q0 in q0s:
        t0 = time.time(); d = run_q0(ps, K, F, q0)
        print(f"K={K} q0={q0} n={d['n']} survivors={d['total']} inclass={d['inclass']} strict={d['strict']} "
              f"outside={d['outside']} [F2={d['F2']} F1={d['F1']}] allowed(v,s)={d['n_allowed']} gcoll={d['n_gcoll']} ({time.time()-t0:.1f}s)", flush=True)
