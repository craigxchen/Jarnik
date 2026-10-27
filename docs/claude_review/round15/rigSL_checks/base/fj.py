# fj.py -- fast joint-core pipeline (same mathematics and certificates as jc2.py/jointcore.py; compact encodings).
# Signatures are sorted tuples of pair ids  e = 4*i + s  (monomial index i, unit s); deep signatures carry
# multiplicities.  Steps for (NAME, K, Qsub):
#   tables  : ./csp mode 3 (level 1, individually admissible labels, optional projection capmax) for F ∩ [3,Qsub];
#             ./deepsig at 3:3, 7:2, 11:2 for witness construction;
#   lower   : ./cpc on C = F ∩ [3,Q] (options: all (signature, eps-mask) pairs on C, minimal signatures elsewhere),
#             'INFEASIBLE' = exact certificate that every good system is coherent on C (Lemma R 2-5);
#   witness : ./cpc with deep tables, extension to all of F (jointcore.extend), exact verification (joint_ilp.verify)
#             and independent re-check (check_system.py) by the caller.
# Usage: python3 fj.py NAME K Qsub [nodecap] [cexp] [Qstart]
import sys, os, json, pickle, time, math, subprocess
from rl import *
from rl import _mk
from deform14 import Setup, modulus_data
from sigtab import class_masks
from jointcore import extend
from joint_ilp import verify

def tab_q(S, q, capmax=None, deep=None):
    m, n, L, iota = modulus_data(S, q); eps = 1 if iota == n else 3
    masks = class_masks(L, n, S.r); par = nu(S.P, q)
    if deep is None:
        M = m; forb = [[(s * iota) % m for s in range(4) if q * q > S.cap[i][s]] for i in range(len(S.V))]
        rec = [sum(1 << t for t in range(4) if capmax is None or S.cap[i][(t * eps) % 4] <= capmax) for i in range(len(S.V))]
        lines = [f"{m} {n} {S.r} 3 0 0", " ".join(map(str, par))]
        col0 = orbit_reps_U1(m, par[0]); lines.append(str(len(col0))); lines += [f"{x} {w}" for x, w in col0]
        lines.append(str(len(S.V)))
        for i, v in enumerate(S.V):
            lines.append(" ".join(map(str, v)) + f" {len(forb[i])} " + " ".join(map(str, forb[i])) + f" {rec[i]}")
        prog = "csp"
    else:
        a = deep; M = m * q ** (a - 1); col0 = orbit_reps_U1(M, par[0])
        lines = [f"{M} {n} {S.r} 3 0 0 {q} {n}", " ".join(map(str, par)), str(len(col0))] + [f"{x} {w}" for x, w in col0]
        lines.append(str(len(S.V)))
        for i, v in enumerate(S.V):
            f = [(s * iota) % m for s in range(4) if q * q > S.cap[i][s]]
            bm = []
            for t in range(4):
                s = (t * eps) % 4; b = 0
                while b < a and q ** (2 * (b + 1)) <= S.cap[i][s]: b += 1
                bm.append(b)
            lines.append(" ".join(map(str, v)) + f" {len(f)} " + " ".join(f"{x} {m}" for x in f) + " " + " ".join(map(str, bm)))
        prog = "deepsig"
    lines.append(str(len(masks))); lines += [" ".join(map(str, c)) + f" {mk}" for c, mk in sorted(masks.items())]
    out = subprocess.run([os.path.join(HERE, prog)], input="\n".join(lines) + "\n", capture_output=True, text=True).stdout
    sigs = []
    for ln in out.splitlines():
        if not ln.startswith("S "): continue
        t = ln.split(); Lr = int(t[5]); r = S.r
        mk = _mk(t, 6 + 2 * r, r, Lr)
        if deep is None:
            el = tuple(sorted(4 * (e // 4) + ((e % 4) * eps) % 4 for e in mk["el"])); mult = None
        else:
            d = {}
            for e in mk["el"]:
                k, rest = divmod(e, 64); tt, b = divmod(rest, 16); d[4 * k + (tt * eps) % 4] = b
            el = tuple(sorted(d)); mult = tuple(d[x] for x in el)
        sigs.append((el, mult, mk["mks"]))
    return dict(q=q, m=m, n=n, L=L, iota=iota, M=M, deep=deep, capmax=capmax, sigs=sigs)

def minimal_ids_deep(sigs):
    """multiplicity-minimal deep signatures (numpy dominance test, preallocated buffer)."""
    import numpy as np
    keys = sorted(set(e for el, mu, _ in sigs for e in el)); ix = {e: i for i, e in enumerate(keys)}
    A = np.zeros((len(sigs), len(keys)), dtype=np.int8)
    for k, (el, mu, _) in enumerate(sigs):
        for e, b in zip(el, mu): A[k, ix[e]] = b
    order = np.argsort(A.sum(axis=1), kind="stable"); keep = []
    Kp = np.zeros((len(sigs), len(keys)), dtype=np.int8); nk = 0
    for k in order:
        if nk and (Kp[:nk] <= A[k][None, :]).all(axis=1).any(): continue
        keep.append(int(k)); Kp[nk] = A[k]; nk += 1
    return keep

def minimal_ids(sigs):
    """indices of inclusion-minimal level-1 signatures (multiplicity-minimal for deep)."""
    order = sorted(range(len(sigs)), key=lambda k: (sum(sigs[k][1]) if sigs[k][1] else len(sigs[k][0])))
    keep = []; index = {}
    for k in order:
        el, mu, _ = sigs[k]; dk = dict(zip(el, mu)) if mu else None; dom = False; sel = set(el)
        for e in el:
            for j in index.get(e, ()):
                ej, mj, _ = sigs[j]
                if (mj is None and sel.issuperset(ej)) or (mj is not None and all(x in sel and m2 <= dk[x] for x, m2 in zip(ej, mj))):
                    dom = True; break
            if dom: break
        if dom: continue
        if not el: return [k]
        keep.append(k); index.setdefault(el[0], []).append(k)
    return keep

def cpc_input(S, tabs, C, full, margin=0.0):
    F = sorted(tabs); opts = {}
    for q in F:
        T = tabs[q]
        if q in C: opts[q] = [(el, mu, mk, rep) for (el, mu, mks) in T["sigs"] for (mk, rep) in mks]
        else: opts[q] = [(T["sigs"][k][0], T["sigs"][k][1], 0, T["sigs"][k][2][0][1]) for k in T["min"]]
    pairs = sorted(set(e for q in F for o in opts[q] for e in o[0])); pid = {e: i for i, e in enumerate(pairs)}
    lines = [str(len(pairs)), " ".join(repr(math.log(S.cap[e // 4][e % 4]) / 2 - margin) for e in pairs), f"{full} {len(F)}"]
    for q in F:
        lines.append(f"{math.log(q)!r} {1 if q in C else 0} {len(opts[q])}")
        for el, mu, mk, rep in opts[q]:
            mu = mu or (1,) * len(el)
            lines.append(f"{mk} {len(el)} " + " ".join(f"{pid[e]} {b}" for e, b in zip(el, mu)))
    return opts, "\n".join(lines) + "\n"

def run_cpc(S, tabs, C, full, nodecap, maxsol=1, margin=0.0):
    opts, inp = cpc_input(S, tabs, C, full, margin); F = sorted(tabs)
    out = subprocess.run([os.path.join(HERE, "cpc"), str(nodecap), str(maxsol)], input=inp, capture_output=True, text=True).stdout.splitlines()
    sols = []
    for ln in out:
        t = ln.split()
        if t[0] == "FOUND": sols.append({q: opts[q][int(k)] for q, k in zip(F, t[2:])})
    last = out[-1].split(); return last[0], int(last[1]), sols

def to_verify(S, tabs, sol):
    """convert a compact solution to the dicts used by jointcore.extend / joint_ilp.verify."""
    ch = {}; vt = {}
    for q, (el, mu, mk, rep) in sol.items():
        T = tabs[q]
        d = dict(sig=frozenset((e // 4, e % 4) for e in el), rep=tuple(rep), repN=tuple(rep))
        if mu: d["mult"] = {(e // 4, e % 4): b for e, b in zip(el, mu)}
        ch[q] = d
        vt[q] = dict(q=q, m=T["m"], n=T["n"], L=T["L"], iota=T["iota"], capmax=T["capmax"])
        if T["deep"]: vt[q].update(M=T["M"], deep=T["deep"])
    return ch, vt

if __name__ == '__main__':
    name, K, Qsub = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    nodecap = int(sys.argv[4]) if len(sys.argv) > 4 else 10 ** 6
    cexp = float(sys.argv[5]) if len(sys.argv) > 5 and sys.argv[5] != "0" else None
    Qstart = int(sys.argv[6]) if len(sys.argv) > 6 else 0
    capmax = int(Qsub ** cexp) if cexp else None
    thr = json.load(open(f"thr_{name}_K{K}.json")); P, F = thr["P"], thr["forced"]
    S = Setup(P, K, F); full = (1 << (2 ** (len(P) - 1))) - 1
    Fs = [q for q in F if q <= Qsub]; t0 = time.time()
    tabs = {}; tabsd = {}
    for q in Fs:
        cf = f"l1cache_{name}_K{K}_{q}_{capmax}.pkl"
        if os.path.exists(cf): tabs[q] = pickle.load(open(cf, "rb"))
        else:
            tabs[q] = tab_q(S, q, capmax); tabs[q]["min"] = minimal_ids(tabs[q]["sigs"])
            pickle.dump(tabs[q], open(cf, "wb"))
    print(f"tables built for {len(Fs)} moduli <= {Qsub} [{time.time()-t0:.0f}s]; sizes " +
          " ".join(f"{q}:{len(tabs[q]['sigs'])}/{len(tabs[q]['min'])}" for q in Fs), flush=True)
    lower = None; upper = None; log = []
    DSETS = (((3, 2), (5, 2), (7, 2)), ((3, 3), (5, 3), (7, 2)), ((3, 3), (5, 3), (7, 2), (11, 2), (13, 2)))
    def deep_tabs(dset):
        td = dict(tabs)
        for qd, a in dset:
            if qd in tabs:
                cf = f"deepcache_{name}_K{K}_{qd}_{a}.pkl"
                if os.path.exists(cf): td[qd] = pickle.load(open(cf, "rb"))
                else:
                    td[qd] = tab_q(S, qd, None, deep=a)
                    td[qd]["min"] = list(range(len(td[qd]["sigs"])))   # no dominance reduction (exactness unaffected)
                    pickle.dump(td[qd], open(cf, "wb"))
                print(f"   deep table q={qd} a={a}: {len(td[qd]['sigs'])} signatures, {len(td[qd]['min'])} minimal", flush=True)
        return td
    cur = tabs; level = 0; dcache = {}
    def try_witness(tb, C, Q):
        st2, nodes2, sols2 = run_cpc(S, tb, C, full, nodecap, maxsol=30)
        for sol in sols2:
            ch, vt = to_verify(S, tb, sol)
            ext = extend(S, dict(tabs=vt), ch)
            if ext is None: continue
            chf, tabsf = ext
            ver = verify(P, S.V, tabsf, chf, None, "mix")
            if ver["ok"]:
                pickle.dump(dict(q0=Q, Q=Q, target="mix", P=P, K=K, choice=chf, labels=ver["labels"],
                                 deep={q: tabsf[q]["deep"] for q in tabsf if tabsf[q].get("deep")}),
                            open(f"jwmix_{name}_K{K}_{Q}.pkl", "wb"))
                return True, len(sols2)
        return False, len(sols2)
    for Q in Fs:
        if Q < Qstart: continue
        t1 = time.time(); C = set(q for q in Fs if q <= Q); done = False
        while True:
            st, nodes, sols = run_cpc(S, cur, C, full, nodecap)
            if st == "INFEASIBLE":
                lower = Q; log.append((Q, "INFEASIBLE", nodes, level))
                print(f"Q={Q}: coherent (exact, depth level {level}) nodes={nodes} [{time.time()-t1:.0f}s]", flush=True)
                break
            ok, nsol = try_witness(cur, C, Q)
            if ok:
                upper = Q; log.append((Q, "WITNESS", nodes, level)); done = True
                print(f"Q={Q}: explicit verified non-coherent good system on all of F (depth level {level}) [{time.time()-t1:.0f}s]", flush=True)
                break
            if level >= len(DSETS):
                log.append((Q, "UNDECIDED-" + st, nodes, level)); done = True
                print(f"Q={Q}: UNDECIDED ({nsol} solutions at full depth, none extended+verified) [{time.time()-t1:.0f}s]", flush=True)
                break
            print(f"Q={Q}: {st} at depth level {level} (nodes={nodes}), no verified witness among {nsol}; deepening", flush=True)
            if level not in dcache: dcache[level] = deep_tabs(DSETS[level])
            cur = dcache[level]; level += 1
        if done: break
    json.dump(dict(name=name, K=K, Qsub=Qsub, capmax=capmax, P=P, core_lower=lower, first_noncoherent=upper, log=log),
              open(f"jointcore_{name}_K{K}.json", "w"))
    print(f"RESULT {name} K={K}: joint core >= {lower} (exact); non-coherent witness at Q={upper}")

# ---------------- streaming variant (memory-light): tables are loaded one at a time from the caches ----------------
def table_file(name, K, q, deep=None, capmax=None):
    return f"deepcache_{name}_K{K}_{q}_{deep}.pkl" if deep else f"l1cache_{name}_K{K}_{q}_{capmax}.pkl"

def run_cpc_stream(S, name, K, qs, C, full, nodecap, maxsol=1, deepmap=None, margin=0.0, tmp="cpc_in.txt"):
    """same search as run_cpc, input streamed to a file; pair ids are e = 4i+s (all pairs listed)."""
    deepmap = deepmap or {}
    npairs = 4 * len(S.V)
    with open(tmp, "w") as fh:
        fh.write(f"{npairs}\n")
        fh.write(" ".join(repr(math.log(S.cap[e // 4][e % 4]) / 2 - margin) for e in range(npairs)) + "\n")
        fh.write(f"{full} {len(qs)}\n")
        for q in qs:
            T = pickle.load(open(table_file(name, K, q, deepmap.get(q)), "rb"))
            if q in C: opts = [(el, mu, mk) for (el, mu, mks) in T["sigs"] for (mk, rep) in mks]
            else: opts = [(T["sigs"][k][0], T["sigs"][k][1], 0) for k in T["min"]]
            fh.write(f"{math.log(q)!r} {1 if q in C else 0} {len(opts)}\n")
            for el, mu, mk in opts:
                mu = mu or (1,) * len(el)
                fh.write(f"{mk} {len(el)} " + " ".join(f"{e} {b}" for e, b in zip(el, mu)) + "\n")
            del T, opts
    out = subprocess.run(f"{os.path.join(HERE, 'cpc')} {nodecap} {maxsol} < {tmp}", shell=True, capture_output=True, text=True).stdout.splitlines()
    sols = []
    for ln in out:
        t = ln.split()
        if t[0] == "FOUND": sols.append(dict(zip(qs, map(int, t[2:]))))
    last = out[-1].split(); return last[0], int(last[1]), sols

def resolve(S, name, K, sol, C, deepmap=None):
    """option indices -> compact (el, mu, mk, rep) per modulus (reloading only the needed tables)."""
    deepmap = deepmap or {}; out = {}; tabs = {}
    for q, k in sol.items():
        T = pickle.load(open(table_file(name, K, q, deepmap.get(q)), "rb"))
        if q in C:
            opts = [(el, mu, mk, rep) for (el, mu, mks) in T["sigs"] for (mk, rep) in mks]
        else:
            opts = [(T["sigs"][j][0], T["sigs"][j][1], 0, T["sigs"][j][2][0][1]) for j in T["min"]]
        out[q] = opts[k]
        tabs[q] = {x: T[x] for x in ("q", "m", "n", "L", "iota", "M", "deep", "capmax")}
    return out, tabs
