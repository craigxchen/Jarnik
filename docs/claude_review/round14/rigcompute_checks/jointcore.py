# jointcore.py -- exact JOINT rigid core of a (prime set, K) (round 14).
#   Q_core = largest forced Q such that EVERY good system (all forced moduli free, every level) is system-level
#            coherent with one column-sign conjugate of the genuine labels on F ∩ [3, Q].
# Lower bound: cp_mix INFEASIBLE on C = F ∩ [3,Q] using only the moduli F ∩ [3, Qsub] (dropping moduli is a
#   relaxation; level-1 relaxation) -> exact certificate.
# Upper bound: an explicit good system on ALL of F (verified exactly at every level by joint_ilp.verify and then
#   independently by check_system.py) that is not coherent on F ∩ [3, Q'] -> Q_core < Q'.  It is built from a
#   cp_mix solution on F ∩ [3, Qsub] (deep tables at the smallest moduli) extended modulus by modulus to the rest
#   of F by an exact DFS for a label whose collisions fit the remaining budget (./csp mode 2).
# Usage: python3 jointcore.py NAME K Qsub [nodecap]
import sys, os, json, pickle, time, math, random
from rl import *
from deform14 import Setup, modulus_data
from sigtab import table
from deep import deep_table
from ls_joint import prepare
from cp_mix import CPmix
from joint_ilp import verify

def build(name, K, Qsub, deep=((3, 3), (7, 2), (11, 2)), capmax=None):
    fn = f"sigj_{name}_K{K}_{Qsub}_{capmax}.pkl"
    if os.path.exists(fn): return pickle.load(open(fn, "rb"))
    thr = json.load(open(f"thr_{name}_K{K}.json")); P, F = thr["P"], thr["forced"]
    S = Setup(P, K, F); tabs = {}
    for q in F:
        if q > Qsub: break
        t0 = time.time(); tabs[q] = table(S, q, capmax)
        print(f"   table q={q}: {len(tabs[q]['sigs'])} signatures [{time.time()-t0:.1f}s]", flush=True)
    D1 = dict(P=P, K=K, F=F, tabs=tabs)
    tabsd = dict(tabs)
    for q, a in deep:
        if q in tabsd: tabsd[q] = deep_table(S, q, a)
    Dd = dict(P=P, K=K, F=F, tabs=tabsd)
    out = dict(D1=D1, Dd=Dd); pickle.dump(out, open(fn, "wb")); return out

def extend(S, Dd, choice, nodecap=10 ** 6):
    """extend a choice on F ∩ [3,Qsub] to all of F; returns (full choice, full tabs) or None."""
    P, K = S.P, S.K; V = S.V
    tabs = dict(Dd["tabs"]); ch = dict(choice)
    # level-1 (+deep multiplicity) loads of the partial system
    prod = {}
    for q, d in ch.items():
        if "mult" in d:
            for e, b in d["mult"].items(): prod[e] = prod.get(e, 1) * q ** b
        else:
            for e in d["sig"]: prod[e] = prod.get(e, 1) * q
    for q in S.F:
        if q in ch: continue
        m, n, L, iota = modulus_data(S, q)
        # leave room for one extra level of each collision (q^2 per collision) to ease lifting
        forb = [[(s * iota) % m for s in range(4) if (prod.get((i, s), 1) * q) ** 2 > S.cap[i][s]] for i in range(len(V))]
        par = nu(P, q); col0 = orbit_reps_U1(m, par[0])
        res, w = run_csp(m, n, S.r, 2, 0, nodecap, par, col0, V, forb, [])
        if res["raw"] == 0: return None
        # mode 2 records no witness line; recover the label by re-running in mode 1 with an empty class set
        res, w = run_csp(m, n, S.r, 1, 1, nodecap, par, col0, V, forb, [])
        ell = list(w[0]); coll = []
        for i, v in enumerate(V):
            t = sum(a * b for a, b in zip(v, ell)) % m
            if t % n == 0:
                s = [u for u in range(4) if (u * iota) % m == t][0]; coll.append((i, s))
                prod[(i, s)] = prod.get((i, s), 1) * q
        tabs[q] = dict(q=q, m=m, n=n, L=L, iota=iota)
        ch[q] = dict(sig=frozenset(coll), rep=tuple(ell), repN=tuple(ell))
    return ch, tabs

if __name__ == '__main__':
    name, K, Qsub = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    nodecap = int(sys.argv[4]) if len(sys.argv) > 4 else 20000
    cexp = float(sys.argv[5]) if len(sys.argv) > 5 else None      # capmax = Qsub^cexp (projection), None = full
    capmax = int(Qsub ** cexp) if cexp else None
    T = build(name, K, Qsub, capmax=capmax)
    D1 = prepare(T["D1"]); Dd = prepare(T["Dd"])
    P, F = D1["P"], D1["F"]; S = Setup(P, K, F)
    Csel = [q for q in F if q <= Qsub]
    lower = None; upper = None; log = []
    for Q in Csel:
        t0 = time.time(); cp = CPmix(D1, Q); res = cp.solve(nodecap)
        if res is None and not cp.aborted and cp.l1sol == 0:
            lower = Q; log.append((Q, "INFEASIBLE", cp.nodes)); print(f"Q={Q}: coherent (exact) nodes={cp.nodes} [{time.time()-t0:.0f}s]", flush=True)
            continue
        print(f"Q={Q}: level-1 relaxation on F<= {Qsub} has a non-coherent system (or aborted={cp.aborted}); building witness", flush=True)
        # witness: deep tables, then extension
        cpd = CPmix(Dd, Q); cpd.maxl1 = 50; resd = cpd.solve(nodecap)
        ok = False
        if resd is not None:
            ext = extend(S, Dd, resd)
            if ext is not None:
                ch, tabs = ext
                ver = verify(P, S.V, tabs, ch, None, "mix")
                if ver["ok"]:
                    ok = True
                    W = dict(q0=Q, Q=Q, target="mix", P=P, K=K, choice=ch, labels=ver["labels"],
                             deep={q: tabs[q]["deep"] for q in tabs if "deep" in tabs[q]})
                    pickle.dump(W, open(f"jwmix_{name}_K{K}_{Q}.pkl", "wb"))
        log.append((Q, "WITNESS" if ok else "UNDECIDED", cp.nodes))
        print(f"Q={Q}: {'explicit verified non-coherent good system on all of F' if ok else 'UNDECIDED'} [{time.time()-t0:.0f}s]", flush=True)
        if ok: upper = Q
        break
    json.dump(dict(name=name, K=K, Qsub=Qsub, capmax=capmax, P=P, core_lower=lower, first_noncoherent=upper, log=log),
              open(f"jointcore_{name}_K{K}.json", "w"))
    print(f"RESULT {name} K={K}: joint core >= {lower} (exact); non-coherent witness at Q={upper}")
