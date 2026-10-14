# fjs.py -- memory-light version of fj.py (same certificates; tables streamed from the caches one at a time).
# Usage: python3 fjs.py NAME K Qsub Qstart Qend [nodecap]
#   tests Q in F ∩ [Qstart, Qend]; certificate = cpc INFEASIBLE (Lemma R); witness verified exactly.
import sys, os, json, pickle, time, math
from fj import *

def ensure_tables(S, name, K, qs):
    for q in qs:
        cf = table_file(name, K, q)
        if not os.path.exists(cf):
            t0 = time.time(); T = tab_q(S, q); T["min"] = minimal_ids(T["sigs"]); pickle.dump(T, open(cf, "wb"))
            print(f"   table q={q}: {len(T['sigs'])} sigs, {len(T['min'])} minimal [{time.time()-t0:.0f}s]", flush=True)
            del T

def ensure_deep(S, name, K, dset):
    for q, a in dset:
        cf = table_file(name, K, q, a)
        if not os.path.exists(cf):
            T = tab_q(S, q, None, deep=a); T["min"] = list(range(len(T["sigs"]))); pickle.dump(T, open(cf, "wb"))
            print(f"   deep table q={q} a={a}: {len(T['sigs'])} sigs", flush=True); del T

if __name__ == '__main__':
    name, K, Qsub, Qstart, Qend = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    nodecap = int(sys.argv[6]) if len(sys.argv) > 6 else 10 ** 6
    thr = json.load(open(f"thr_{name}_K{K}.json")); P, F = thr["P"], thr["forced"]
    S = Setup(P, K, F); full = (1 << (2 ** (len(P) - 1))) - 1
    qs = [q for q in F if q <= Qsub]; ensure_tables(S, name, K, qs)
    DSETS = (((3, 2), (5, 2), (7, 2)), ((3, 3), (5, 3), (7, 2)), ((3, 3), (5, 3), (7, 2), (11, 2), (13, 2)))
    level = 0; deepmap = {}; lower = None; upper = None; log = []
    for Q in qs:
        if Q < Qstart or Q > Qend: continue
        t1 = time.time(); C = set(q for q in qs if q <= Q); done = False
        while True:
            st, nodes, sols = run_cpc_stream(S, name, K, qs, C, full, nodecap, deepmap=deepmap, tmp=f"cpcin_{name}_{K}.txt")
            if st == "INFEASIBLE":
                lower = Q; log.append((Q, "INFEASIBLE", nodes, level))
                print(f"Q={Q}: coherent (exact, depth level {level}) nodes={nodes} [{time.time()-t1:.0f}s]", flush=True); break
            st2, n2, sols2 = run_cpc_stream(S, name, K, qs, C, full, nodecap, maxsol=30, deepmap=deepmap, tmp=f"cpcin_{name}_{K}.txt")
            ok = False
            for sol in sols2:
                comp, vt = resolve(S, name, K, sol, C, deepmap)
                ch, vt2 = to_verify(S, vt, comp)
                ext = extend(S, dict(tabs=vt2), ch)
                if ext is None: continue
                chf, tabsf = ext
                ver = verify(P, S.V, tabsf, chf, None, "mix")
                if ver["ok"]:
                    pickle.dump(dict(q0=Q, Q=Q, target="mix", P=P, K=K, choice=chf, labels=ver["labels"],
                                     deep={q: tabsf[q]["deep"] for q in tabsf if tabsf[q].get("deep")}),
                                open(f"jwmix_{name}_K{K}_{Q}.pkl", "wb"))
                    ok = True; break
            if ok:
                upper = Q; log.append((Q, "WITNESS", nodes, level)); done = True
                print(f"Q={Q}: explicit verified non-coherent good system on all of F (depth level {level}) [{time.time()-t1:.0f}s]", flush=True); break
            if level >= len(DSETS):
                log.append((Q, "UNDECIDED", nodes, level)); done = True
                print(f"Q={Q}: UNDECIDED [{time.time()-t1:.0f}s]", flush=True); break
            print(f"Q={Q}: {st} at depth level {level}; deepening", flush=True)
            dset = [(q, a) for q, a in DSETS[level] if q in qs]; ensure_deep(S, name, K, dset)
            deepmap = dict(dset); level += 1
        if done: break
    json.dump(dict(name=name, K=K, Qsub=Qsub, P=P, core_lower=lower, first_noncoherent=upper, log=log),
              open(f"jointcore_s_{name}_K{K}_{Qsub}.json", "w"))
    print(f"RESULT {name} K={K} Qsub={Qsub}: coherent up to {lower} (exact, tested from {Qstart}); witness at {upper}")
