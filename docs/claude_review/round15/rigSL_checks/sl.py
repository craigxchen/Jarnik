# sl.py -- round 15: Rig^SL (labels + ONE angle vector phi) on top of the round-14 PU harness (copied to base/).
# Pair ids e = 4*i + s  (i = index of the sign-normalised monomial V[i], s = unit).  (-v, s) == (v, -s).
# Budget of (v,s) at phi (log scale):  B_phi(e) = 0.5*log(N_v * 2^[s odd]) + log|sin(<v,phi> - s pi/4)|.
# (1.2) of sharp.md for (v,s)  <=>  log m_(v,s) <= B_phi(e).  At phi = arg pi_j: B_phi(e) = log|D_s(Z_v)|.
import sys, os, json, math, pickle, subprocess, time
HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(HERE, "base")
sys.path.insert(0, BASE)
from rl import *                       # noqa  (copied round-14 helpers)
from deform14 import Setup, modulus_data
import fj

CACHE = os.path.join(HERE, "cache"); os.makedirs(CACHE, exist_ok=True)

def load_set(name, K):
    fn = os.path.join(BASE, f"thr_{name}_K{K}.json")
    if not os.path.exists(fn):          # forced set not needed (e.g. pinning tests): primes only
        return json.load(open(os.path.join(BASE, f"thr_{name}_K4.json")))["P"], []
    thr = json.load(open(fn))
    return thr["P"], thr["forced"]

def setup(name, K):
    P, F = load_set(name, K)
    S = Setup(P, K, F)
    S.name = name
    S.N = [1] * len(S.V)
    for i, v in enumerate(S.V):
        N = 1
        for p, t in zip(P, v): N *= p ** abs(t)
        S.N[i] = N
    S.phig = [math.atan2(b, a) for (a, b) in S.pis]
    return S

def tables(S, qs):
    """PU level-1 tables (individually PU-admissible labels) + inclusion-minimal ids; cached."""
    tabs = {}
    for q in qs:
        cf = os.path.join(CACHE, f"l1_{S.name}_K{S.K}_{q}.pkl")
        if os.path.exists(cf):
            tabs[q] = pickle.load(open(cf, "rb")); continue
        T = fj.tab_q(S, q, None); T["min"] = fj.minimal_ids(T["sigs"])
        pickle.dump(T, open(cf, "wb")); tabs[q] = T
    return tabs

def deep_tab(S, q, a):
    cf = os.path.join(CACHE, f"deep_{S.name}_K{S.K}_{q}_{a}.pkl")
    if os.path.exists(cf): return pickle.load(open(cf, "rb"))
    T = fj.tab_q(S, q, None, deep=a); T["min"] = list(range(len(T["sigs"])))
    pickle.dump(T, open(cf, "wb")); return T

def xv(S, phi):
    return [sum(a * b for a, b in zip(v, phi)) for v in S.V]

def budgets_PU(S):
    return [0.5 * math.log(S.N[i] * (2 if s % 2 else 1)) for i in range(len(S.V)) for s in range(4)]

def budgets_phi(S, phi, margin=0.0):
    x = xv(S, phi); out = []
    for i in range(len(S.V)):
        for s in range(4):
            sn = abs(math.sin(x[i] - s * math.pi / 4))
            out.append(0.5 * math.log(S.N[i] * (2 if s % 2 else 1)) + (math.log(sn) if sn > 0 else -1e9) - margin)
    return out

def budgets_genuine(S):
    return [math.log(abs(S.D[i][s])) for i in range(len(S.V)) for s in range(4)]

def sup_abs_sin(lo, hi):
    """sup of |sin t| over [lo, hi] (float; caller adds a margin)."""
    if hi - lo >= math.pi: return 1.0
    k = math.ceil((lo - math.pi / 2) / math.pi)
    if math.pi / 2 + k * math.pi <= hi: return 1.0
    return max(abs(math.sin(lo)), abs(math.sin(hi)))

def budgets_box(S, lo, hi, margin=1e-9):
    """valid upper bounds of B_phi(e) over the box prod [lo_j, hi_j] (+margin in favour of feasibility)."""
    out = []
    for i, v in enumerate(S.V):
        a = sum(t * (lo[j] if t > 0 else hi[j]) for j, t in enumerate(v))
        b = sum(t * (hi[j] if t > 0 else lo[j]) for j, t in enumerate(v))
        for s in range(4):
            sn = sup_abs_sin(a - s * math.pi / 4, b - s * math.pi / 4)
            out.append(0.5 * math.log(S.N[i] * (2 if s % 2 else 1)) + (math.log(sn) if sn > 0 else -1e9) + margin)
    return out

def cpc_write(S, tabs, C, full, Bud, fn):
    F = sorted(tabs); opts = {}
    npairs = 4 * len(S.V)
    with open(fn, "w") as fh:
        fh.write(f"{npairs}\n" + " ".join(repr(b) for b in Bud) + "\n" + f"{full} {len(F)}\n")
        for q in F:
            T = tabs[q]
            if q in C: o = [(el, mu, mk, rep) for (el, mu, mks) in T["sigs"] for (mk, rep) in mks]
            else: o = [(T["sigs"][k][0], T["sigs"][k][1], 0, T["sigs"][k][2][0][1]) for k in T["min"]]
            opts[q] = o
            fh.write(f"{math.log(q)!r} {1 if q in C else 0} {len(o)}\n")
            for el, mu, mk, rep in o:
                mu = mu or (1,) * len(el)
                fh.write(f"{mk} {len(el)} " + " ".join(f"{e} {b}" for e, b in zip(el, mu)) + "\n")
    return opts

def run_cpc(S, tabs, C, full, Bud, nodecap=10 ** 6, maxsol=1, fn=None, timeout=None):
    fn = fn or os.path.join(CACHE, f"cpcin_{os.getpid()}.txt")
    opts = cpc_write(S, tabs, C, full, Bud, fn); F = sorted(tabs)
    out = subprocess.run(f"{os.path.join(BASE, 'cpc')} {nodecap} {maxsol} < {fn}", shell=True, capture_output=True,
                         text=True, timeout=timeout).stdout.splitlines()
    sols = []
    for ln in out:
        t = ln.split()
        if t and t[0] == "FOUND": sols.append({q: opts[q][int(k)] for q, k in zip(F, t[2:])})
    last = out[-1].split(); return last[0], int(last[1]), sols

def full_mask(S):
    return (1 << (2 ** (S.r - 1))) - 1

# ---------------- witnesses at a FIXED phi (caps replaced by the phi-caps) ----------------
import joint_ilp, jointcore
def set_caps(S, capmat):
    """replace S.cap by capmat[i][s] (squared budgets, PU-convention: m^2 <= cap) everywhere in the base code."""
    S.cap = capmat
    capd = {}
    for i, v in enumerate(S.V):
        for s in range(4):
            capd[(tuple(v), s)] = capmat[i][s]
            capd[(tuple(-t for t in v), (-s) % 4)] = capmat[i][s]
    joint_ilp.cap = lambda P, v, s: capd[(tuple(v), s % 4)]

def caps_genuine(S):
    return [[S.D[i][s] ** 2 for s in range(4)] for i in range(len(S.V))]

def caps_phi(S, phi, shrink=1 - 1e-9):
    x = xv(S, phi)
    return [[S.N[i] * (2 if s % 2 else 1) * math.sin(x[i] - s * math.pi / 4) ** 2 * shrink for s in range(4)]
            for i in range(len(S.V))]

def tables_tag(S, qs, tag):
    tabs = {}
    for q in qs:
        cf = os.path.join(CACHE, f"l1_{S.name}_K{S.K}_{q}_{tag}.pkl")
        if os.path.exists(cf): tabs[q] = pickle.load(open(cf, "rb")); continue
        T = fj.tab_q(S, q, None); T["min"] = fj.minimal_ids(T["sigs"])
        pickle.dump(T, open(cf, "wb")); tabs[q] = T
    return tabs

def deep_tag(S, q, a, tag):
    cf = os.path.join(CACHE, f"deep_{S.name}_K{S.K}_{q}_{a}_{tag}.pkl")
    if os.path.exists(cf): return pickle.load(open(cf, "rb"))
    T = fj.tab_q(S, q, None, deep=a); T["min"] = list(range(len(T["sigs"])))
    pickle.dump(T, open(cf, "wb")); return T

def budgets_from_caps(S, capmat):
    return [0.5 * math.log(capmat[i][s]) if capmat[i][s] > 0 else -1e9 for i in range(len(S.V)) for s in range(4)]

def build_witness(S, capmat, tag, Q, Qsub, dsets=(((3, 3), (5, 3), (7, 2)), ((3, 3), (5, 3), (7, 2), (11, 2), (13, 2))),
                  nodecap=10 ** 6, maxsol=30, log=print):
    """search a good (w.r.t. capmat) system on all of F, non-coherent on F ∩ [3,Q]; returns verify() dict or None."""
    set_caps(S, capmat)
    qs = [q for q in S.F if q <= Qsub]; tabs = tables_tag(S, qs, tag)
    Bud = budgets_from_caps(S, capmat); full = full_mask(S); C = set(q for q in qs if q <= Q)
    st, nodes, sols = run_cpc(S, tabs, C, full, Bud, nodecap)
    log(f"   level-1 on F<={Qsub}, Q={Q}: {st} nodes={nodes}")
    if st == "INFEASIBLE": return dict(status="INFEASIBLE")
    for ds in dsets:
        td = dict(tabs)
        for qd, a in ds:
            if qd in td: td[qd] = deep_tag(S, qd, a, tag)
        st, nodes, sols = run_cpc(S, td, C, full, Bud, nodecap, maxsol=maxsol)
        log(f"   deep {ds}: {st} nodes={nodes} sols={len(sols)}")
        if st == "INFEASIBLE": return dict(status="INFEASIBLE_DEEP", dset=ds)
        for sol in sols:
            ch, vt = fj.to_verify(S, td, sol)
            ext = jointcore.extend(S, dict(tabs=vt), ch)
            if ext is None: continue
            chf, tabsf = ext
            ver = joint_ilp.verify(S.P, S.V, tabsf, chf, None, "mix")
            if ver["ok"]:
                return dict(status="WITNESS", labels=ver["labels"], choice=chf,
                            deep={q: tabsf[q].get("deep") for q in tabsf if tabsf[q].get("deep")})
            else: log(f"      verify failed: {ver['why']}")
    return dict(status="UNDECIDED")
