# ls_joint.py -- min-conflicts local search for JOINT witnesses (round 14): a good system (every forced modulus free,
# level-1 relaxation during search) whose label at q0 is a target label ('none': outside every column-sign class).
# Moves: switch one modulus to another achievable (minimal) signature, weighted min-conflicts with noise.
# Every solution is converted to explicit labels and verified EXACTLY (joint_ilp.verify: integer budgets over all
# forced moduli, collision sets recomputed from the labels, level-2 lifts with no level-2 collision).
# A failed search proves nothing (evidence only).
# Usage: python3 ls_joint.py SIGPKL target q0list [restarts] [steps] [seed]
import sys, pickle, math, time, random
import numpy as np
from scipy.sparse import csr_matrix
from sigtab import minimal
from rl import cap, monos_norm
from joint_ilp import verify

def minimal_deep(sigs, key=lambda d: True):
    """minimal elements under componentwise multiplicity dominance (deep tables)."""
    cand = [d for d in sigs if key(d)]
    if not cand: return []
    keys = sorted(set(e for d in cand for e in d["mult"])); ix = {e: i for i, e in enumerate(keys)}
    A = np.zeros((len(cand), len(keys)), dtype=np.int16)
    for k, d in enumerate(cand):
        for e, b in d["mult"].items(): A[k, ix[e]] = b
    order = np.argsort(A.sum(axis=1), kind="stable"); mins = []
    for k in order:
        if mins and (A[mins] <= A[k][None, :]).all(axis=1).any(): continue
        mins.append(k)
    return [cand[k] for k in mins]

def prepare(D):
    if "fams" in D: return D
    tabs = D["tabs"]; F = sorted(tabs)
    mn = lambda q, key=(lambda d: True): (minimal_deep if "deep" in tabs[q] else minimal)(tabs[q]["sigs"], key)
    D["fams"] = {q: mn(q) for q in F}
    D["fams_none"] = {q: mn(q, lambda d: d["none"]) for q in F}
    D["fams_notL"] = {q: mn(q, lambda d: d["notL"]) for q in F}
    from deform14 import Setup
    S = Setup(D["P"], D["K"], F)
    D["gsig"] = {q: frozenset(S.val[q].keys()) for q in F}
    for q in F:
        assert any(d["sig"] <= D["gsig"][q] for d in D["fams"][q]), q
    return D

def search(D, q0, target="none", restarts=20, steps=20000, seed=0, verbose=False):
    P, K, tabs = D["P"], D["K"], D["tabs"]; V = monos_norm(len(P), K); F = sorted(tabs)
    fams = dict(D["fams"]); fams[q0] = D["fams_" + target][q0]
    if not fams[q0]: return dict(status="no-target-label")
    pairs = sorted(set(e for q in F for d in fams[q] for e in d["sig"]))
    pid = {e: i for i, e in enumerate(pairs)}; npair = len(pairs)
    Bv = np.array([math.log(cap(P, V[i], s)) / 2 for (i, s) in pairs]) - 1e-9
    mats = {}
    for q in F:
        rows, cols = [], []
        for k, d in enumerate(fams[q]):
            for e in d["sig"]: rows.append(k); cols.append(pid[e])
        mats[q] = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(fams[q]), npair))
    lq = {q: math.log(q) for q in F}
    # genuine-like start: cheapest signature by plain cost
    rng = random.Random(seed); nrs = np.random.default_rng(seed)
    best = None
    for rs in range(restarts):
        choice = {}
        load = np.zeros(npair)
        for q in F:
            if q != q0:
                gs = D["gsig"][q]
                ks = [k for k, d in enumerate(fams[q]) if d["sig"] <= gs]
                k = rng.choice(ks)
            else:
                k = rs % len(fams[q]) if rs < len(fams[q]) else rng.randrange(len(fams[q]))
            choice[q] = k
            load += lq[q] * mats[q][k].toarray().ravel()
        w = np.ones(npair)
        for it in range(steps):
            over = load - Bv
            viol = np.flatnonzero(over > 0)
            if len(viol) == 0:
                ch = {q: fams[q][choice[q]] for q in F}
                ver = verify(P, V, tabs, ch, q0, target)
                return dict(status="found", restart=rs, step=it, choice=ch, verified=ver)
            e = int(rng.choice(viol))
            # moduli contributing to e
            contrib = [q for q in F if mats[q][choice[q], e] != 0]
            q = rng.choice(contrib)
            A = mats[q]; base = load - lq[q] * A[choice[q]].toarray().ravel()
            f = w * (np.maximum(0, base + lq[q] - Bv) - np.maximum(0, base - Bv))
            score = A @ f
            if rng.random() < 0.05: k = rng.randrange(A.shape[0])
            else:
                mn = score.min(); cand = np.flatnonzero(score <= mn + 1e-12); k = int(rng.choice(list(cand)))
            choice[q] = k; load = base + lq[q] * A[k].toarray().ravel()
            w[viol] += 0.05
        nv = int((load - Bv > 0).sum())
        if verbose: print(f"  restart {rs}: {nv} violated after {steps} steps", flush=True)
        if best is None or nv < best: best = nv
    return dict(status="not-found", best_violations=best)

if __name__ == '__main__':
    D = prepare(pickle.load(open(sys.argv[1], "rb"))); target = sys.argv[2]
    q0s = [int(x) for x in sys.argv[3].split(",")]
    restarts = int(sys.argv[4]) if len(sys.argv) > 4 else 10; steps = int(sys.argv[5]) if len(sys.argv) > 5 else 5000
    seed = int(sys.argv[6]) if len(sys.argv) > 6 else 0
    for q0 in q0s:
        t0 = time.time(); o = search(D, q0, target, restarts, steps, seed)
        v = o.get("verified")
        print(f"q0={q0} target={target} {o['status']} {('restart %d step %d' % (o['restart'], o['step'])) if o['status']=='found' else o.get('best_violations')} "
              f"verified={v['ok'] if v else None}{(' maxload=%.4f' % v['max_load']) if v and v['ok'] else ''} [{time.time()-t0:.0f}s]", flush=True)
        if v and v["ok"]:
            pickle.dump(dict(q0=q0, target=target, P=D["P"], K=D["K"], choice=o["choice"], labels=v["labels"]),
                        open(f"jw_{sys.argv[1][4:-4]}_{q0}_{target}.pkl", "wb"))
