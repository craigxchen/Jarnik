# cp_mix.py -- EXACT test of system-level coherence on an initial segment C = F ∩ [3, Q] (round 14).
# Question: is there a good system (every forced modulus free) that is NOT coherent on C with any column-sign
# conjugate of the genuine labels, i.e. with AND_{q in C} mask(ell(q)) = 0, where mask(ell) = set of eps (eps_0 = +1;
# the global sign is a unit) with ell mod n in class(u * eps * L)?
# Options at q in C: pairs (signature, mask) realised by some label, reduced by dominance (sigma' <= sigma and
# mask' <= mask); elsewhere minimal signatures.  Complete DFS with capacity propagation (as cp_joint.py) plus:
# prune when some eps-bit is contained in the mask of every viable option at every q in C (it can no longer be
# removed).  INFEASIBLE => every good system is coherent on C with one conjugate genuine system (exact; level-1
# relaxation, sound tolerances).  FOUND => explicit system, verified exactly at all levels (joint_ilp.verify).
# Usage: python3 cp_mix.py SIGPKL Qlist [nodecap]
import sys, os, pickle, math, time
import numpy as np
from ls_joint import prepare, minimal_deep
from sigtab import minimal
from cp_joint import CP, TOL
from joint_ilp import verify

def mask_options(T, deep):
    """expand signatures into (sig, mask) entries, keep dominance-minimal ones."""
    ent = []
    for d in T["sigs"]:
        for mk, rep in d["mks"]:
            e = dict(d); e["mask"] = mk; e["rep"] = rep; e["repN"] = rep; ent.append(e)
    # dominance: same-or-smaller signature and submask
    ent.sort(key=lambda e: (len(e["sig"]), bin(e["mask"]).count("1")))
    keep = []
    for e in ent:
        dom = False
        for x in keep:
            if (x["mask"] & ~e["mask"]) == 0 and (x["sig"] <= e["sig"]) and \
               (not deep or all(x["mult"].get(k, 0) <= e["mult"].get(k, 0) for k in x["mult"])):
                dom = True; break
        if not dom: keep.append(e)
    return keep

class CPmix(CP):
    def __init__(self, D, Q, Qsub=None):
        tabs = D["tabs"]; F = sorted(q for q in tabs if Qsub is None or q <= Qsub)
        self.C = [q for q in F if q <= Q]
        D2 = dict(D); D2["fams"] = dict(D["fams"])
        for q in self.C:
            D2["fams"][q] = mask_options(tabs[q], "deep" in tabs[q])
        D2["fams_mix"] = D2["fams"]
        super().__init__(D2, self.C[0], "mix", Fsub=F)
        self.masks = {q: [d.get("mask", -1) for d in self.fams[q]] for q in self.F}
        self.r = len(D["P"]); self.full = (1 << (2 ** (self.r - 1))) - 1

    def dfs(self, viable):
        self.nodes += 1
        if self.nodes > self.nodecap: self.aborted = True; return None
        viable = self.propagate(viable)
        if viable is None: return None
        # bits that every viable option at some q in C ... must be killable: bit e killable iff some q in C has a
        # viable option without e
        killable = 0
        for q in self.C:
            for k in viable[q]: killable |= (~self.masks[q][k]) & self.full
        if killable != self.full: return None
        open_ = [q for q in self.F if len(viable[q]) > 1]
        if not open_:
            andm = self.full
            for q in self.C: andm &= self.masks[q][viable[q][0]]
            if andm: return None
            ch = {q: self.fams[q][viable[q][0]] for q in self.F}
            self.l1sol += 1
            ver = verify(self.P, self.V, self.tabs, ch, None, "mix")
            if ver["ok"]: self.ver = ver; return ch
            if self.l1sol >= self.maxl1: self.aborted = True
            return None
        # branch first on segment moduli (mask choice), smallest domain first
        openC = [q for q in open_ if q in self.C]
        pool = openC if openC else open_
        q = min(pool, key=lambda x: (len(viable[x]), -x))
        for k in viable[q]:
            nv = dict(viable); nv[q] = [k]
            res = self.dfs(nv)
            if res is not None or self.aborted: return res
        return None

if __name__ == '__main__':
    D = prepare(pickle.load(open(sys.argv[1], "rb")))
    Qs = [int(x) for x in sys.argv[2].split(",")]; nodecap = int(sys.argv[3]) if len(sys.argv) > 3 else 10 ** 5
    Qsub = int(sys.argv[4]) if len(sys.argv) > 4 else None
    for Q in Qs:
        t0 = time.time(); cp = CPmix(D, Q, Qsub); res = cp.solve(nodecap)
        if res is None:
            st = ("ABORTED" if cp.l1sol == 0 else "UNLIFTED-L1-ONLY") if cp.aborted else \
                 ("INFEASIBLE (exact): every good system is coherent on C" if cp.l1sol == 0 else "ONLY-UNLIFTED-L1")
            print(f"Q={Q} |C|={len(cp.C)}: {st} nodes={cp.nodes} l1sol={cp.l1sol} [{time.time()-t0:.0f}s]", flush=True)
        else:
            print(f"Q={Q} |C|={len(cp.C)}: NON-COHERENT good system FOUND, verified={cp.ver['ok']} nodes={cp.nodes} [{time.time()-t0:.0f}s]", flush=True)
            pickle.dump(dict(q0=None, Q=Q, target="mix", P=cp.P, K=cp.K, choice=res, labels=cp.ver["labels"],
                             deep={q: cp.tabs[q]["deep"] for q in cp.F if "deep" in cp.tabs[q]}),
                        open(f"jwmix_{os.path.basename(sys.argv[1])[:-4]}_{Q}.pkl", "wb"))
