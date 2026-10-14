# repair.py -- extend a verified joint witness from norm K-1 to norm K by exact local repair (round 14).
# Start: the labels (with lifts) of jwmix_<set>_K<K-1>_<Q>.pkl on F_(K-1); add the moduli of F_K \ F_(K-1).
# Loop: compute every collision of every v (||v||_1 <= K, sign-normalised) at every level from the stored top labels;
# if some (v,s) is over budget, pick a modulus q contributing to it (never one whose label must stay as it is) and
# re-choose its label by exact DFS (./csp) among labels whose level-1 collisions fit the budget left by the OTHER
# moduli (plus one spare level), lift it greedily, repeat.  The target modulus Q keeps a label outside every
# column-sign class (re-chosen, if needed, among such labels).  Success = a system passing joint_ilp.verify-style exact
# checks; it is then saved as jwmix_<set>_K<K>_<Q>.pkl and must be re-checked with check_system.py.
# A failure proves nothing.
# Usage: python3 repair.py NAME K Q [maxiter] [seed]
import sys, json, pickle, random, math, itertools
from rl import *
from deform14 import Setup, modulus_data
from sigtab import class_masks

def collisions_top(S, q, top, A):
    """per (i,s): number of levels b at which v_i collides, from the top label mod m q^(A-1); also whether it persists."""
    m = mq(q); n = m // 4; T = S.T[q]; iota = T[3]; M = m * q ** (A - 1); out = {}; persist = False
    for i, v in enumerate(S.V):
        d = sum(a * b for a, b in zip(v, top)) % M
        if d % n: continue
        s = [u for u in range(4) if (u * iota) % m == d % m][0]; b = 1
        while b < A and d % (n * q ** b) == 0: b += 1
        if b == A and d % (n * q ** (A - 1)) == 0 and A >= 1:
            # collision at the deepest stored level: deeper levels undetermined
            persist = True
        out[(i, s)] = b
    return out, persist

def lift_label(S, q, ell, slack_fn, rng, maxlev=12, tries=3000):
    """greedy level-by-level lift of a level-1 label; returns (top, A) or None."""
    m = mq(q); n = m // 4; r = S.r; cur = list(ell); mod = m; A = 1
    for lev in range(2, maxlev + 1):
        active = [i for i, v in enumerate(S.V) if sum(a * b for a, b in zip(v, cur)) % (mod // 4) == 0]
        if not active: return cur, A
        best = None
        cands = itertools.product(range(q), repeat=r) if q ** r <= 20000 else (tuple(rng.randrange(q) for _ in range(r)) for _ in range(tries))
        for tt in cands:
            l2 = [cur[j] + mod * tt[j] for j in range(r)]
            hit = [i for i in active if sum(a * b for a, b in zip(S.V[i], l2)) % (mod * q // 4) == 0]
            if best is None or len(hit) < len(best[1]): best = (l2, hit)
            if not hit: break
        cur, hit = best; mod *= q; A = lev
        if not hit: return cur, A
    return None

def main():
    name, K, Q = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    maxiter = int(sys.argv[4]) if len(sys.argv) > 4 else 300; rng = random.Random(int(sys.argv[5]) if len(sys.argv) > 5 else 1)
    W = pickle.load(open(f"jwmix_{name}_K{K-1}_{Q}.pkl", "rb"))
    thr = json.load(open(f"thr_{name}_K{K}.json")); P, F = thr["P"], thr["forced"]
    S = Setup(P, K, F); S.T = {}
    for q in F:
        m, n, L, iota = modulus_data(S, q); S.T[q] = (m, n, L, iota)
    lab = {}
    for q in F:
        if q in W["labels"]:
            l = W["labels"][q]; A = len(l["lifts"]) + (W.get("deep", {}).get(q, 1) - 1); lab[q] = (list(l["lifts"][-1]), A)
    def loads():
        prod = {}; per = {}
        for q, (top, A) in lab.items():
            c, _ = collisions_top(S, q, top, A); per[q] = c
            for k, b in c.items(): prod[k] = prod.get(k, 1) * q ** b
        return prod, per
    def choose(q, prod, per, need_none):
        m, n, L, iota = S.T[q]
        other = dict(prod)
        for k, b in per.get(q, {}).items(): other[k] //= q ** b
        par = nu(P, q); col0 = orbit_reps_U1(m, par[0]); rng.shuffle(col0)
        cls = sorted(class_masks(L, n, S.r).keys()) if need_none else []
        w = None
        for spare in (q, 1, None):
            if spare is None:
                forb = [[] for i in range(len(S.V))]
            else:
                forb = [[(s * iota) % m for s in range(4) if (other.get((i, s), 1) * q * spare) ** 2 > S.cap[i][s]] for i in range(len(S.V))]
            res, w = run_csp(m, n, S.r, 1, 1, 2 * 10 ** 6, par, col0, S.V, forb, cls)
            if w: break
        if not w: return None
        lifted = lift_label(S, q, list(w[0]), None, rng)
        return lifted
    for q in F:                                      # new forced moduli
        if q not in lab:
            prod, per = loads(); x = choose(q, prod, per, False)
            if x is None: print(f"cannot place new modulus {q}"); return
            lab[q] = x
    for it in range(maxiter):
        prod, per = loads()
        bad = [k for k, p in prod.items() if p * p > S.cap[k[0]][k[1]]]
        if not bad:
            break
        k = rng.choice(bad); qs = [q for q in lab if k in per[q]]
        q = rng.choice(qs)
        x = choose(q, prod, per, q == Q)
        if x is not None: lab[q] = x
        if it % 20 == 0: print(f"iter {it}: {len(bad)} violated pairs", flush=True)
    prod, per = loads()
    bad = [k for k, p in prod.items() if p * p > S.cap[k[0]][k[1]]]
    if bad: print(f"FAILED: {len(bad)} violated pairs after {maxiter} iterations"); return
    labels = {q: dict(ell=lab[q][0][:], lifts=[lab[q][0]]) for q in lab}
    deep = {q: lab[q][1] for q in lab}
    pickle.dump(dict(q0=Q, Q=Q, target="mix", P=P, K=K, labels=labels, deep=deep, how="repair from K-1"),
                open(f"jwmix_{name}_K{K}_{Q}.pkl", "wb"))
    print(f"SUCCESS: saved jwmix_{name}_K{K}_{Q}.pkl (re-check with check_system.py)")

if __name__ == '__main__':
    main()
