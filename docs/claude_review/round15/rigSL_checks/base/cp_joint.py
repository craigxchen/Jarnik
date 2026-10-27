# cp_joint.py -- EXACT joint test (round 14): is there a good system, every forced modulus free, whose label at q0 is
# a target label?  Complete DFS over per-modulus signatures with capacity propagation:
#   options(q) = inclusion-minimal achievable signatures at q (q0: minimal among target-achievable ones);
#   mandatory load: every element common to all still-viable options of an unassigned modulus is charged;
#   an option is dropped when one of its elements cannot absorb log q on top of the load of the OTHER moduli.
# Soundness of 'INFEASIBLE': (1) level-1 relaxation (deeper levels only add demand); (2) minimal signatures suffice
# (demand is monotone in the signature); (3) capacity tests use log sums with tolerance +1e-9 in favour of
# feasibility, so every truly feasible choice survives every test; (4) the DFS is complete.
# A 'FOUND' choice is converted to labels and verified exactly (joint_ilp.verify).
# Usage: python3 cp_joint.py SIGPKL target q0list [nodecap]
import sys, os, pickle, math, time
import numpy as np
from ls_joint import prepare
from rl import cap, monos_norm
from joint_ilp import verify

TOL = 1e-9

class CP:
    def __init__(self, D, q0, target, Fsub=None):
        self.P, self.K, self.tabs = D["P"], D["K"], D["tabs"]
        self.V = monos_norm(len(self.P), self.K)
        self.F = sorted(Fsub or self.tabs)
        fams = {q: D["fams"][q] for q in self.F}; fams[q0] = D["fams_" + target][q0]
        self.fams = fams; self.q0 = q0
        pairs = sorted(set(e for q in self.F for d in fams[q] for e in d["sig"]))
        self.pid = {e: i for i, e in enumerate(pairs)}; self.np = len(pairs)
        self.B = np.array([math.log(cap(self.P, self.V[i], s)) / 2 for (i, s) in pairs])
        self.lq = {q: math.log(q) for q in self.F}
        self.opts = {q: [sum(1 << self.pid[e] for e in d["sig"]) for d in fams[q]] for q in self.F}
        # deep moduli: full multiplicity matrices (options x pairs)
        self.deepM = {}
        for q in self.F:
            if "deep" in self.tabs[q]:
                A = np.zeros((len(fams[q]), self.np))
                for k, d in enumerate(fams[q]):
                    for e, b in d["mult"].items(): A[k, self.pid[e]] = b
                self.deepM[q] = A
        self.nodes = 0; self.l1sol = 0; self.maxl1 = 200; self.target = target; self.ver = None

    def vec(self, mask):
        out = np.zeros(self.np); i = 0
        while mask:
            if mask & 1: out[i] = 1
            mask >>= 1; i += 1
        return out

    def tomask(self, boolarr):
        idx = np.flatnonzero(boolarr); m = 0
        for i in idx: m |= 1 << int(i)
        return m

    def propagate(self, viable):
        """viable: dict q -> list of option indices.  Returns pruned dict or None."""
        viable = {q: list(v) for q, v in viable.items()}
        while True:
            common = {}
            load = np.zeros(self.np)
            own = {}
            for q in self.F:
                if q in self.deepM:
                    own[q] = self.lq[q] * self.deepM[q][viable[q]].min(axis=0); load += own[q]; common[q] = 1; continue
                c = -1
                for k in viable[q]: c &= self.opts[q][k]
                common[q] = c if viable[q] else 0
                own[q] = self.lq[q] * self.vec(common[q]) if common[q] else 0
                if common[q]: load += own[q]
            if (load > self.B + TOL).any(): return None
            changed = False
            for q in self.F:
                if len(viable[q]) <= 1:
                    continue
                if q in self.deepM:
                    other = load - own[q]
                    ok = ((other[None, :] + self.lq[q] * self.deepM[q][viable[q]]) <= self.B[None, :] + TOL).all(axis=1)
                    nv = [k for k, o in zip(viable[q], ok) if o]
                else:
                    over = self.tomask(load - own[q] + self.lq[q] > self.B + TOL)
                    nv = [k for k in viable[q] if not (self.opts[q][k] & over)]
                if not nv: return None
                if len(nv) != len(viable[q]): viable[q] = nv; changed = True
            if not changed: return viable

    def solve(self, nodecap=10 ** 7):
        self.nodecap = nodecap; self.aborted = False
        v0 = {q: list(range(len(self.opts[q]))) for q in self.F}
        return self.dfs(v0)

    def dfs(self, viable):
        self.nodes += 1
        if self.nodes > self.nodecap: self.aborted = True; return None
        viable = self.propagate(viable)
        if viable is None: return None
        open_ = [q for q in self.F if len(viable[q]) > 1]
        if not open_:
            ch = {q: self.fams[q][viable[q][0]] for q in self.F}
            self.l1sol += 1
            ver = verify(self.P, self.V, self.tabs, ch, self.q0, self.target)
            if ver["ok"]: self.ver = ver; return ch
            if self.l1sol >= self.maxl1: self.aborted = True
            return None
        q = min(open_, key=lambda x: (len(viable[x]), -x))
        for k in viable[q]:
            nv = dict(viable); nv[q] = [k]
            res = self.dfs(nv)
            if res is not None or self.aborted: return res
        return None

if __name__ == '__main__':
    D = prepare(pickle.load(open(sys.argv[1], "rb"))); target = sys.argv[2]
    q0s = [int(x) for x in sys.argv[3].split(",")]; nodecap = int(sys.argv[4]) if len(sys.argv) > 4 else 10 ** 6
    for q0 in q0s:
        t0 = time.time()
        if not D["fams_" + target][q0]:
            print(f"q0={q0} target={target}: NO target label at all (local)", flush=True); continue
        cp = CP(D, q0, target); res = cp.solve(nodecap)
        if res is None:
            st = ("ABORTED" if cp.l1sol == 0 else "LEVEL1-SOLUTIONS-UNLIFTED") if cp.aborted else ("INFEASIBLE (exact)" if cp.l1sol == 0 else "ONLY-UNLIFTED-L1-SOLUTIONS")
            print(f"q0={q0} target={target}: {st} nodes={cp.nodes} l1sol={cp.l1sol} [{time.time()-t0:.0f}s]", flush=True)
        else:
            ver = cp.ver
            print(f"q0={q0} target={target}: FOUND nodes={cp.nodes} verified={ver['ok']} {ver.get('why','')} [{time.time()-t0:.0f}s]", flush=True)
            if ver["ok"]:
                pickle.dump(dict(q0=q0, target=target, P=cp.P, K=cp.K, choice=res, labels=ver["labels"],
                                 deep={q: cp.tabs[q]["deep"] for q in cp.F if "deep" in cp.tabs[q]}),
                            open(f"jw_{os.path.basename(sys.argv[1])[:-4]}_{q0}_{target}.pkl", "wb"))
