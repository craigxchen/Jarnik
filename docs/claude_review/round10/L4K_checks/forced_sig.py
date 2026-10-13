# Forced region, exhaustive at level 1: per forced prime q, all parity-respecting label vectors ell in (Z/m)^r,
# grouped by collision signature sig(ell) = {(v, s): <v,ell> = s*iota mod m}, v in V_K (one of each +-pair).
# Cost of a signature = log q on each (v,s) it contains.  A system on the forced primes is good iff for every (v,s)
# the summed cost <= B(v,s) = w(v)/2 + [s odd] log2/2.
import sys, math, itertools, time, pickle, re
import numpy as np
from lib import *

def sig_table(P, K, q):
    r = len(P); T = Torus(q); m = T.m; iota = T.iota; n = m // 4
    V = monos(r, K); Va = np.array(V, dtype=np.int64)          # (nv, r)
    par = nu(P, q)
    cols = [np.arange(par[j], m, 2, dtype=np.int64) for j in range(r)]
    grids = np.meshgrid(*cols, indexing='ij')
    ells = np.stack([g.ravel() for g in grids], axis=1)        # (N, r)
    dots = (ells @ Va.T) % m                                    # (N, nv)
    # unit index s with s*iota = dot  (iota = n or 3n)
    s_of = {(s * iota) % m: s for s in range(4)}
    coll = (dots % n == 0)
    sidx = np.full(dots.shape, -1, dtype=np.int64)
    for val, s in s_of.items():
        sidx[dots == val] = s
    code = np.where(coll, sidx, -1)                              # (N, nv): -1 or s
    # group
    keys, inv, counts = np.unique(code, axis=0, return_inverse=True, return_counts=True)
    L = np.array([T.L(gauss_prime(p)) for p in P], dtype=np.int64)
    gl = ((L @ Va.T) % m)
    gcode = np.where(gl % n == 0, np.array([s_of[x] if x in s_of else -1 for x in gl]), -1)
    gk = int(np.where((keys == gcode).all(axis=1))[0][0])
    return dict(q=q, m=m, n=n, iota=iota, V=V, keys=keys, counts=counts, ells=ells, inv=inv.ravel(), genuine_key=gk, L=L)

def budgets(P, V):
    B = np.zeros((len(V), 4))
    for i, v in enumerate(V):
        for s in range(4):
            B[i, s] = budget(P, v, s)
    return B

def count_good(P, K, F, tables, B, limit_systems=None):
    nv = len(tables[0]['V'])
    # cost of each key: matrix (nv,4) of log q
    costs = []
    for t in tables:
        lq = math.log(t['q']); C = []
        for key in t['keys']:
            c = np.zeros((nv, 4))
            for i, s in enumerate(key):
                if s >= 0: c[i, s] = lq
            C.append(c)
        costs.append(C)
    order = sorted(range(len(tables)), key=lambda i: len(tables[i]['keys']))
    good_tuples = []  # list of tuples of key indices (in table order)
    total = [0]
    eps = 1e-9
    def dfs(k, load, choice, mult):
        if k == len(order):
            total[0] += mult
            good_tuples.append((dict(choice), mult))
            return
        i = order[k]; t = tables[i]
        for kk, c in enumerate(costs[i]):
            nl = load + c
            if (nl <= B + eps).all():
                choice[i] = kk
                dfs(k + 1, nl, choice, mult * int(t['counts'][kk]))
                del choice[i]
    dfs(0, np.zeros((nv, 4)), {}, 1)
    return total[0], good_tuples

if __name__ == '__main__':
    K = int(sys.argv[1]); r = int(sys.argv[2]); tag = sys.argv[3]
    allp = split_primes(5, 200)
    import random
    rng0 = random.Random(7); sets = {}
    for rr in range(4, 9):
        sets[('first', rr)] = allp[:rr]; sets[('last', rr)] = allp[-rr:]; sets[('rand', rr)] = sorted(rng0.sample(allp, rr))
    P = sets[(tag, r)]
    for line in open('free_out.txt'):
        mm = re.match(r"K=(\d) r=(\d) (\w+) P=.*?non-free q = (\[.*?\]);", line)
        if mm and int(mm[1]) == K and int(mm[2]) == r and mm[3] == tag: F = eval(mm[4])
    t0 = time.time()
    tables = [sig_table(P, K, q) for q in F]
    V = tables[0]['V']; B = budgets(P, V)
    print(f"K={K} r={r} {tag} P={P} forced primes F={F}")
    for t in tables:
        print(f"  q={t['q']} m={t['m']} vectors={len(t['ells'])} signatures={len(t['keys'])} genuine-signature class size={t['counts'][t['genuine_key']]}")
    tot, good = count_good(P, K, F, tables, B)
    print(f"  good level-1 systems on F: {tot} (signature tuples: {len(good)}) [{time.time()-t0:.1f}s]")
    gen_tuple = {i: t['genuine_key'] for i, t in enumerate(tables)}
    n_gen = sum(m for ch, m in good if ch == gen_tuple)
    print(f"  of which with the genuine signature at every q: {n_gen}")
    # per-q: number of good systems whose signature at q is non-genuine
    for i, t in enumerate(tables):
        nn = sum(m for ch, m in good if ch[i] != t['genuine_key'])
        print(f"    q={t['q']}: good systems with non-genuine signature at q: {nn}")
    pickle.dump(dict(P=P, K=K, F=F, good=good, tables=[{k: t[k] for k in ('q','m','n','iota','keys','counts','genuine_key','L')} for t in tables]),
                open(f'forced_sig_{K}_{r}_{tag}.pkl', 'wb'))
