# onemod_angle.py -- exhaustive one-modulus test at q0 (genuine labels at all other forced moduli, all levels;
# level-1 relaxation at q0), in two classes:
#   PU  : per-unit budgets only (the write-up's class);
#   ANG : per-unit budgets + per-v four-unit angle compatibility (a necessary condition of (1.2) in sharp.md Def 1.1).
# Reports survivors (all / outside the genuine class / outside every column-sign class) in each class.
import sys, json, math, itertools
import numpy as np
from reflib import *
from rverify import per_v_angle_ok

def run(name, K, q0, Fsrc):
    thr = json.load(open(Fsrc)); P = thr["P"]; r = len(P); F = sorted(set(thr["forced"]) | set(thr.get("aborted", [])))
    pis = [gprime(p) for p in P]
    m = m_of(q0); n = m // 4
    L, iota, _ = genuine_L(P, q0); k = (iota // n) % 4
    V = monos(r, K, half=True)
    G = np.zeros((len(V), 4)); Gi = [[1] * 4 for _ in V]; capi = [[capv(P, v, s) for s in range(4)] for v in V]
    for i, v in enumerate(V):
        Dv = Dvals(Zv(pis, v))
        for s in range(4):
            g = 1
            for q in F:
                if q != q0:
                    t = val(Dv[s], q)
                    if t: g *= q ** t
            Gi[i][s] = g
    # allowed collisions (v,s) at q0: (g*q0)^2 <= cap
    allowed = np.array([[(Gi[i][s] * q0) ** 2 <= capi[i][s] for s in range(4)] for i in range(len(V))])
    Vm = np.array(V, dtype=np.int64)
    nu = nu_par(P, q0)
    cols = [np.arange(nu[j], m, 2, dtype=np.int64) for j in range(r)]
    surv = []
    U1 = [u for u in range(1, m, 4) if math.gcd(u, m) == 1]
    seen = set(); reps = []
    for x in cols[0]:
        x = int(x)
        if x in seen: continue
        orb = set((u * x) % m for u in U1); seen |= orb; reps.append((x, len(orb)))
    rest = np.array(list(itertools.product(*cols[1:])), dtype=np.int64)
    for x0, wgt in reps:
        ell = np.concatenate([np.full((len(rest), 1), x0, dtype=np.int64), rest], axis=1)
        d = (ell @ Vm.T) % m  # (N, |V|)
        coll = (d % n) == 0
        s = ((d // n) * (1 if k == 1 else 3)) % 4
        ok = np.ones(len(ell), dtype=bool)
        for t in range(4):
            bad = coll & (s == t) & (~allowed[:, t][None, :])
            ok &= ~bad.any(axis=1)
        for idx in np.nonzero(ok)[0]:
            surv.append((tuple(int(x) for x in ell[idx]), wgt))
    gen = set(tuple((u * x) % n for x in L) for u in units(n))
    res = dict(PU=[0, 0, 0], ANG=[0, 0, 0])
    angsurv = []
    for e, wgt in surv:
        c = tuple(x % n for x in e)
        mask = class_mask(c, L, n, r)
        cls = (c not in gen, mask == 0)
        res["PU"][0] += wgt; res["PU"][1] += wgt * cls[0]; res["PU"][2] += wgt * cls[1]
        # angle compatibility with level-1 demand
        okang = True
        for i, v in enumerate(V):
            dd = sum(a * b for a, b in zip(v, e)) % m
            dem = list(Gi[i])
            if dd % n == 0:
                st = ((dd // n) * (1 if k == 1 else 3)) % 4
                dem[st] *= q0
            Nv = norm_v(P, v)
            if not angle_ok_exact(dem, Nv):
                okang = False; break
        if okang:
            res["ANG"][0] += wgt; res["ANG"][1] += wgt * cls[0]; res["ANG"][2] += wgt * cls[1]
            if cls[0]: angsurv.append(e)
    return res, angsurv[:3]

if __name__ == "__main__":
    name, K = sys.argv[1], int(sys.argv[2])
    Fsrc = sys.argv[4] if len(sys.argv) > 4 else f"../rigcompute_checks/thr_{name}_K{K}.json"
    qs = [int(x) for x in sys.argv[3].split(",")]
    if len(qs) == 2 and qs[0] < 0:
        F = sorted(json.load(open(Fsrc))["forced"]); qs = [q for q in F if -qs[0] <= q <= qs[1]]
    for q0 in qs:
      res, ex = run(name, K, q0, Fsrc)
      print(name, "K", K, "q0", q0, "survivors [all, outside genuine class, outside every eps-class]:",
          "per-unit", res["PU"], "| per-unit+angle", res["ANG"], "| example ANG non-genuine", ex, flush=True)
