import sys, ctypes, math, itertools, numpy as np
from glib import *
from forced import forced_set
lib = ctypes.CDLL('./onemodc.so'); P = lambda a: a.ctypes.data_as(ctypes.c_void_p)
lib.run.argtypes = [ctypes.c_int]*3 + [ctypes.c_void_p]*6
def canon(v):
    for x in v:
        if x: return x > 0
def onemod(ps, K, F, q0, full=False):
    r = len(ps); pis = [two_sq(p) for p in ps]; V = [v for v in VK(r, K) if canon(v)]
    m, L = genuine(ps, q0); n = m//4
    others = [q for q in F if q != q0]
    allow = np.zeros(4*len(V), np.int32)
    for k, v in enumerate(V):
        z = Zv(pis, v); nz = z[0]**2+z[1]**2
        for s in range(4):
            D = Ds(z, s); g = 1
            if full:   # every odd prime != q0 not dividing N charged: g = |D| stripped of 2, p_j, q0
                g = abs(D)
                while g % 2 == 0: g //= 2
                for p in ps + [q0]:
                    while g % p == 0: g //= p
            else:
                for q in others: g *= q**vq(D, q)
            allow[4*k+s] = 1 if (q0*g)**2 <= nz*(2 if s % 2 else 1) else 0
    gc = np.array([1 if sum(a*b for a, b in zip(v, L)) % n == 0 else 0 for v in V], np.int32)
    cl = np.zeros(n**r, np.uint8)
    c = [x % n for x in L]
    for u in range(1, n+1):
        if math.gcd(u, n) == 1:
            idx = sum(((u*x) % n)*n**t for t, x in enumerate(c)); cl[idx] = 1
    par = np.array(nu(ps, q0), np.int32); Va = np.array(V, np.int32).ravel()
    res = np.zeros(4, np.int64); w = np.zeros(r, np.int32)
    lib.run(r, m, len(V), P(Va), P(par), P(allow), P(gc), P(cl), P(res), P(w))
    return [int(x) for x in res], [int(x) for x in w], L
if __name__ == "__main__":
    ps = [int(x) for x in sys.argv[1].split(",")]; K = int(sys.argv[2]); qmax = int(sys.argv[3])
    full = len(sys.argv) > 5 and sys.argv[5] == "full"
    F, _ = forced_set(ps, K, int(sys.argv[4]))
    for q0 in F:
        if q0 > qmax: break
        res, w, L = onemod(ps, K, F, q0, full)
        print(f"K={K} q0={q0} n={(q0-1 if q0%4==1 else q0+1)//4} surv={res[0]} nonclass={res[1]} F1={res[2]} F2={res[3]} wit={w if res[1] else '-'}", flush=True)
