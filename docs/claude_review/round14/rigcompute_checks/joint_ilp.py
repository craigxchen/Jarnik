# joint_ilp.py -- search (HiGHS MILP) for a good system, everything free, with a target label at q0; any solution is
# turned into explicit labels (the stored representatives realise the chosen signatures exactly) and verified
# EXACTLY: integer budgets (prod q)^2 <= cap(v,s) over all forced moduli, plus level-2 lifts without level-2 collisions.
# A MILP 'infeasible' answer is NOT a certificate (reported as evidence only); certificates come from joint_lp.py
# or exhaustive search.
import sys, pickle, math, time, random
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix, csr_matrix
from sigtab import minimal
from rl import cap, monos_norm, nu

def solve(D, q0, target, time_limit=600, fams=None, scale=1.0):
    P, K, tabs = D["P"], D["K"], D["tabs"]; V = monos_norm(len(P), K); F = sorted(tabs)
    if fams is None:
        fams = {}
        for q in F:
            if q == q0:
                key = {"none": lambda d: d["none"], "notL": lambda d: d["notL"]}[target]
                fams[q] = minimal(tabs[q]["sigs"], key)
            else:
                fams[q] = minimal(tabs[q]["sigs"])
    if not fams[q0]: return dict(status="no-target-label")
    pairs = sorted(set(e for q in F for d in fams[q] for e in d["sig"]))
    pid = {e: i for i, e in enumerate(pairs)}
    cols = [(q, d) for q in F for d in fams[q]]
    nx = len(cols)
    A = lil_matrix((len(F) + len(pairs), nx)); lo = []; hi = []
    qidx = {q: i for i, q in enumerate(F)}
    for j, (q, d) in enumerate(cols):
        A[qidx[q], j] = 1
        for e in d["sig"]: A[len(F) + pid[e], j] = math.log(q)
    lo = [1] * len(F) + [-np.inf] * len(pairs)
    hi = [1] * len(F) + [math.log(cap(P, V[i], s)) / 2 * scale + 1e-9 for (i, s) in pairs]
    t0 = time.time()
    res = milp(c=np.zeros(nx), constraints=LinearConstraint(csr_matrix(A), lo, hi), integrality=np.ones(nx),
               bounds=Bounds(0, 1), options=dict(time_limit=time_limit, disp=False))
    out = dict(q0=q0, target=target, milp_status=res.status, msg=res.message, t=time.time() - t0, ncols=nx)
    if res.x is None: return out
    choice = {}
    for j in np.where(res.x > 0.5)[0]:
        q, d = cols[j]; choice[q] = d
    out["choice"] = choice
    out["verified"] = verify(P, V, tabs, choice, q0, target)
    return out

def verify(P, V, tabs, choice, q0, target, rng=None, maxlev=12):
    """exact: labels = stored reps; recompute level-1 collisions from scratch; then lift level by level, choosing at
    each level a lift t in (Z/q)^r (exhaustive if q^r <= 20000, else 4000 random tries) that minimises the number of
    deeper collisions that the remaining budget cannot absorb; every deeper collision is charged exactly.
    Final test: for every (v,s), (prod of all colliding prime powers)^2 <= cap(v,s)."""
    import itertools
    rng = rng or random.Random(1)
    r = len(P); prod = {}; labels = {}; coll1 = {}; deep = {}
    for q, d in choice.items():
        T = tabs[q]; m, n, iota = T["m"], T["n"], T["iota"]; Mq = T.get("M", m)
        ell = list(d["repN"] if (q == q0) else d["rep"])
        assert all(ell[j] % 2 == x for j, x in enumerate(nu(P, q)))
        coll = []; deepact = []
        for i, v in enumerate(V):
            dd = sum(a * b for a, b in zip(v, ell)) % Mq; t = dd % m
            if t % n == 0:
                s = [u for u in range(4) if (u * iota) % m == t][0]; coll.append((i, s))
                b = 1; md = n * q
                while md <= Mq // 4 and dd % md == 0: b += 1; md *= q
                prod[(i, s)] = prod.get((i, s), 1) * q ** b
                if "M" in T:
                    assert d["mult"][(i, s)] == b, "deep rep does not realise multiplicities"
                    if md > Mq // 4: deepact.append((i, s))
        cm = T.get("capmax")
        assert frozenset(e for e in coll if cm is None or cap(P, V[e[0]], e[1]) <= cm) == d["sig"], "rep does not realise signature"
        coll1[q] = (ell, coll if "M" not in T else deepact, Mq)
    for q, (ell, coll, Mq) in sorted(coll1.items()):
        T = tabs[q]; m, n = T["m"], T["n"]
        cur = list(ell); mod = Mq; active = list(coll); levels = [cur]
        for lev in range(2 + T.get('deep', 1) - 1, maxlev + 1):
            if not active: break
            nmod = mod * q
            cands = itertools.product(range(q), repeat=r) if q ** r <= 20000 else (tuple(rng.randrange(q) for _ in range(r)) for _ in range(4000))
            best = None
            for tt in cands:
                l2 = [cur[j] + mod * tt[j] for j in range(r)]
                hit = [(i, s) for (i, s) in active if sum(a * b for a, b in zip(V[i], l2)) % (nmod // 4) == 0]
                bad = sum(1 for (i, s) in hit if (prod[(i, s)] * q) ** 2 > cap(P, V[i], s))
                key = (bad, len(hit))
                if best is None or key < best[0]: best = (key, l2, hit)
                if key == (0, 0): break
            _, cur, hit = best; levels.append(cur)
            for k in hit: prod[k] *= q; deep.setdefault(k, []).append((q, lev))
            active = hit; mod = nmod
        if active: return dict(ok=False, why=f"collisions persist to level {maxlev} at {q}")
        labels[q] = dict(ell=ell, lifts=levels)
    bad = [(k, p) for k, p in prod.items() if p * p > cap(P, V[k[0]], k[1])]
    if bad: return dict(ok=False, why=f"{len(bad)} budget violations", bad=[(V[k[0]], k[1], deep.get(k)) for k, p in bad])
    return dict(ok=True, labels=labels, max_load=max(math.log(p) / (math.log(cap(P, V[k[0]], k[1])) / 2)
                                                   for k, p in prod.items()))

if __name__ == '__main__':
    D = pickle.load(open(sys.argv[1], "rb")); target = sys.argv[2]
    q0s = [int(x) for x in sys.argv[3].split(",")]; tl = float(sys.argv[4]) if len(sys.argv) > 4 else 600
    for q0 in q0s:
        o = solve(D, q0, target, tl)
        v = o.get("verified")
        print(f"q0={q0} target={target} milp_status={o.get('milp_status')} ({o.get('msg')}) cols={o.get('ncols')} "
              f"t={o.get('t', 0):.1f}s verified={v['ok'] if v else None} {'' if not v or v['ok'] else v['why']}"
              f"{' maxload=%.4f' % v['max_load'] if v and v['ok'] else ''}", flush=True)
        if v and v["ok"]:
            pickle.dump(o, open(f"jointwit_{sys.argv[1][4:-4]}_{q0}_{target}.pkl", "wb"))
