# sigtab.py -- per-modulus tables of achievable collision signatures (level 1) for the joint tests (round 14).
# At each forced q: enumerate EVERY parity-respecting ell in (Z/m)^r (column 0 up to U1-orbits, exact weights) whose
# level-1 collisions (v,s) are all individually admissible (q^2 <= cap(v,s)); labels with an inadmissible collision
# can never occur in a good system.  Record each distinct signature sigma = {(v,s)} with:
#   orm  = OR of class masks of its labels (bit e: c in class(u * eps_e * L), eps_e = column sign vector, eps_0 = +1)
#   notL = some label with this signature is outside class(L);  none = some label is outside every eps-class
#   rep / repN = a representative label (repN outside class(L) when notL).
# Usage: python3 sigtab.py THRJSON [qmax]  -> sig_<name>_K<K>.pkl
import sys, os, json, time, pickle, math
from rl import *
from deform14 import Setup, modulus_data

def eps_vectors(r):
    return [[1] + [(-1 if (b >> (j - 1)) & 1 else 1) for j in range(1, r)] for b in range(2 ** (r - 1))]

def class_masks(L, n, r):
    Lc = [x % n for x in L]; masks = {}
    for e, eps in enumerate(eps_vectors(r)):
        for c in genuine_class([(a * x) % n for a, x in zip(eps, Lc)], n):
            masks[c] = masks.get(c, 0) | (1 << e)
    return masks

def table(S, q, capmax=None):
    m, n, L, iota = modulus_data(S, q)
    eps = 1 if iota == n else 3
    forb = [[(s * iota) % m for s in range(4) if q * q > S.cap[i][s]] for i in range(len(S.V))]
    par = nu(S.P, q); col0 = orbit_reps_U1(m, par[0])
    masks = class_masks(L, n, S.r)
    rec = None
    if capmax is not None:
        rec = [sum(1 << t for t in range(4) if S.cap[i][(t * eps) % 4] <= capmax) for i in range(len(S.V))]
    res, sigs = run_csp(m, n, S.r, 3, 0, 0, par, col0, S.V, forb, masks, rec=rec)
    assert res["status"] == "complete"
    out = []
    for sg in sigs:
        el = frozenset((e // 4, ((e % 4) * eps) % 4) for e in sg["el"])   # residue index t -> unit s = t*eps mod 4
        out.append(dict(sig=el, orm=sg["orm"], none=sg["none"], notL=sg["notL"], cnt=sg["cnt"], rep=sg["rep"],
                        repN=sg["repN"], mks=sg["mks"]))
    return dict(q=q, m=m, n=n, L=L, iota=iota, sigs=out, leaves=res["leaves"], capmax=capmax)

def minimal(sigs, key=lambda d: True):
    """inclusion-minimal elements (by signature set) among the entries satisfying key."""
    cand = [d for d in sigs if key(d)]
    cand.sort(key=lambda d: len(d["sig"]))
    mins = []; index = {}          # element -> list of minimal sigs whose 'anchor' (min element) is it
    for d in cand:
        s = d["sig"]; dom = False
        for e in s:
            for x in index.get(e, ()):
                if x <= s: dom = True; break
            if dom: break
        if dom: continue
        mins.append(d)
        if s: index.setdefault(min(s), []).append(s)
        else: index.setdefault(None, []).append(s)
    if None in index:                # empty signature present: it dominates everything
        return [d for d in mins if not d["sig"]]
    return mins

if __name__ == '__main__':
    thr = json.load(open(sys.argv[1])); P, K, F = thr["P"], thr["K"], thr["forced"]
    qmax = int(sys.argv[2]) if len(sys.argv) > 2 else 10 ** 9
    S = Setup(P, K, F); tabs = {}
    for q in F:
        if q > qmax: break
        t0 = time.time(); T = table(S, q); tabs[q] = T
        mn = minimal(T["sigs"]); mnN = minimal(T["sigs"], lambda d: d["notL"])
        print(f"q={q} n={T['n']} admissible labels={T['leaves']} signatures={len(T['sigs'])} minimal={len(mn)} "
              f"minimal-notL={len(mnN)} min|sig|={min(len(d['sig']) for d in mn)} [{time.time()-t0:.1f}s]", flush=True)
    name = os.path.basename(sys.argv[1])[4:-5]
    pickle.dump(dict(P=P, K=K, F=F, tabs=tabs), open(f"sig_{name}.pkl", "wb"))
