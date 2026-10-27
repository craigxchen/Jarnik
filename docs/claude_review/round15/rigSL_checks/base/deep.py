# deep.py -- multi-level signature tables at small moduli (round 14).  A 'deep' label at q is ell mod M = m q^(a-1);
# its signature is the multiset {(v,s): b} where b = number of levels <= a at which v collides (cost b log q).
# Leaves are kept only if q^(2b) <= cap(v,s) for every collision (individual admissibility, exact).
# Usage: python3 deep.py SIGPKL q1:a1,q2:a2,...   -> rewrites SIGPKL with tabs[q]['deep'] = a (backup kept)
import sys, os, pickle, subprocess, time, shutil
from rl import *
from rl import _mk
from deform14 import Setup, modulus_data
from sigtab import class_masks

def deep_table(S, q, a):
    m, n, L, iota = modulus_data(S, q); eps = 1 if iota == n else 3
    M = m * q ** (a - 1)
    par = nu(S.P, q); col0 = orbit_reps_U1(M, par[0])
    masks = class_masks(L, n, S.r)
    lines = [f"{M} {n} {S.r} 3 0 0 {q} {n}", " ".join(map(str, par)), str(len(col0))]
    lines += [f"{x} {w}" for x, w in col0]
    lines.append(str(len(S.V)))
    for i, v in enumerate(S.V):
        f = [(s * iota) % m for s in range(4) if q * q > S.cap[i][s]]
        bm = []
        for t in range(4):
            s = (t * eps) % 4; b = 0
            while b < a and q ** (2 * (b + 1)) <= S.cap[i][s]: b += 1
            bm.append(b)
        lines.append(" ".join(map(str, v)) + f" {len(f)} " + " ".join(f"{x} {m}" for x in f) + " " + " ".join(map(str, bm)))
    lines.append(str(len(masks)))
    lines += [" ".join(map(str, c)) + f" {mk}" for c, mk in sorted(masks.items())]
    out = subprocess.run([os.path.join(HERE, "deepsig")], input="\n".join(lines) + "\n", capture_output=True, text=True).stdout
    sigs = []
    for ln in out.splitlines():
        if not ln.startswith("S "): continue
        t = ln.split(); Lr = int(t[5]); r = S.r
        rep = tuple(int(x) for x in t[6:6 + r]); repN = tuple(int(x) for x in t[6 + r:6 + 2 * r])
        mk = _mk(t, 6 + 2 * r, r, Lr)
        el = {}
        for e in mk["el"]:
            k, rest = divmod(e, 64); tt, b = divmod(rest, 16)
            el[(k, (tt * eps) % 4)] = b
        sigs.append(dict(sig=frozenset(el), mult=el, orm=int(t[1]), none=int(t[2]), notL=int(t[3]), cnt=int(t[4]),
                         rep=rep, repN=repN, mks=mk["mks"]))
    return dict(q=q, m=m, n=n, L=L, iota=iota, M=M, deep=a, sigs=sigs)

if __name__ == '__main__':
    fn = sys.argv[1]; D = pickle.load(open(fn, "rb"))
    if not os.path.exists(fn + ".level1"): shutil.copy(fn, fn + ".level1")
    S = Setup(D["P"], D["K"], D["F"])
    for tok in sys.argv[2].split(","):
        q, a = map(int, tok.split(":")); t0 = time.time()
        T = deep_table(S, q, a); D["tabs"][q] = T
        print(f"q={q} depth={a} M={T['M']} signatures={len(T['sigs'])} [{time.time()-t0:.1f}s]", flush=True)
    for k in ("fams", "fams_none", "fams_notL", "gsig"): D.pop(k, None)
    pickle.dump(D, open(fn, "wb"))
