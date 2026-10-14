# jc2.py -- joint rigid core via the C search ./cpc (same certificates as jointcore.py; see its header).
# Options: moduli outside C -> inclusion-minimal signatures (deep tables: multiplicity-minimal); moduli in C -> every
# realised (signature, eps-mask) pair (no dominance reduction needed for exactness).
# Usage: python3 jc2.py NAME K Qsub [nodecap] [cexp]
import sys, os, json, pickle, time, math, subprocess
from rl import *
from deform14 import Setup
from ls_joint import prepare
from jointcore import build, extend
from joint_ilp import verify

def options(D, C):
    """per modulus list of option dicts (with 'mask'); C = set of segment moduli."""
    out = {}
    for q in sorted(D["tabs"]):
        T = D["tabs"][q]
        if q in C:
            ent = []
            for d in T["sigs"]:
                for mk, rep in d["mks"]:
                    e = dict(d); e["mask"] = mk; e["rep"] = rep; e["repN"] = rep; ent.append(e)
            out[q] = ent
        else:
            out[q] = [dict(d, mask=-1) for d in D["fams"][q]]
    return out

def run_cpc(D, opts, full, nodecap, maxsol=1):
    P, K = D["P"], D["K"]; V = monos_norm(len(P), K); F = sorted(opts)
    pairs = sorted(set(e for q in F for d in opts[q] for e in d["sig"])); pid = {e: i for i, e in enumerate(pairs)}
    lines = [str(len(pairs)), " ".join(repr(math.log(cap(P, V[i], s)) / 2) for (i, s) in pairs), f"{full} {len(F)}"]
    for q in F:
        lines.append(f"{math.log(q)!r} {1 if opts[q][0]['mask'] != -1 else 0} {len(opts[q])}")
        for d in opts[q]:
            mult = d.get("mult") or {e: 1 for e in d["sig"]}
            lines.append(f"{d['mask'] if d['mask'] != -1 else 0} {len(mult)} " + " ".join(f"{pid[e]} {b}" for e, b in mult.items()))
    out = subprocess.run([os.path.join(HERE, "cpc"), str(nodecap), str(maxsol)], input="\n".join(lines) + "\n",
                         capture_output=True, text=True).stdout.splitlines()
    sols = []
    for ln in out:
        t = ln.split()
        if t[0] == "FOUND": sols.append({q: opts[q][int(k)] for q, k in zip(F, t[2:])})
    last = out[-1].split()
    return last[0], int(last[1]), sols

if __name__ == '__main__':
    name, K, Qsub = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    nodecap = int(sys.argv[4]) if len(sys.argv) > 4 else 10 ** 6
    cexp = float(sys.argv[5]) if len(sys.argv) > 5 and sys.argv[5] != "0" else None
    capmax = int(Qsub ** cexp) if cexp else None
    T = build(name, K, Qsub, capmax=capmax)
    D1 = prepare(T["D1"]); Dd = prepare(T["Dd"]); P, F = D1["P"], D1["F"]; S = Setup(P, K, F)
    full = (1 << (2 ** (len(P) - 1))) - 1
    Csel = [q for q in F if q <= Qsub]; lower = None; upper = None; log = []
    for Q in Csel:
        t0 = time.time(); C = set(q for q in Csel if q <= Q)
        st, nodes, sols = run_cpc(D1, options(D1, C), full, nodecap)
        if st == "INFEASIBLE":
            lower = Q; log.append((Q, "INFEASIBLE", nodes))
            print(f"Q={Q}: coherent (exact) nodes={nodes} [{time.time()-t0:.0f}s]", flush=True); continue
        print(f"Q={Q}: {st} on F<={Qsub} (nodes={nodes}); building witness", flush=True)
        ok = False
        st2, nodes2, sols2 = run_cpc(Dd, options(Dd, C), full, nodecap, maxsol=30)
        for sol in sols2:
            ext = extend(S, Dd, sol)
            if ext is None: continue
            ch, tabs = ext
            ver = verify(P, S.V, tabs, ch, None, "mix")
            if ver["ok"]:
                ok = True
                pickle.dump(dict(q0=Q, Q=Q, target="mix", P=P, K=K, choice=ch, labels=ver["labels"],
                                 deep={q: tabs[q]["deep"] for q in tabs if "deep" in tabs[q]}),
                            open(f"jwmix_{name}_K{K}_{Q}.pkl", "wb"))
                break
        log.append((Q, "WITNESS" if ok else ("UNDECIDED-" + st), nodes))
        print(f"Q={Q}: {'explicit verified non-coherent good system on all of F' if ok else 'UNDECIDED'} [{time.time()-t0:.0f}s]", flush=True)
        if ok: upper = Q
        break
    json.dump(dict(name=name, K=K, Qsub=Qsub, capmax=capmax, P=P, core_lower=lower, first_noncoherent=upper, log=log),
              open(f"jointcore_{name}_K{K}.json", "w"))
    print(f"RESULT {name} K={K}: joint core >= {lower} (exact); non-coherent witness at Q={upper}")
