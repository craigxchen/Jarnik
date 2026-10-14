# Forced region: generate good systems (level 1 + exact lifts through all levels), classify coherence.
# Usage: python3 forced_coh.py K r tag nsys H [seed]
import sys, math, random, re, time, pickle, itertools
import numpy as np
from lib import *

def load_set(K, r, tag):
    allp = split_primes(5, 200); rng0 = random.Random(7); sets = {}
    for rr in range(4, 9):
        sets[('first', rr)] = allp[:rr]; sets[('last', rr)] = allp[-rr:]; sets[('rand', rr)] = sorted(rng0.sample(allp, rr))
    P = sets[(tag, r)]
    for line in open('free_out.txt'):
        mm = re.match(r"K=(\d) r=(\d) (\w+) P=.*?non-free q = (\[.*?\]);", line)
        if mm and int(mm[1]) == K and int(mm[2]) == r and mm[3] == tag: F = eval(mm[4])
    for line in open('free_K4_out.txt'):
        mm = re.match(r"K=(\d) r=(\d) (\w+) P=.*?forced primes q<=1200: (\[.*?\]);", line)
        if mm and int(mm[1]) == K and int(mm[2]) == r and mm[3] == tag: F = eval(mm[4])
    return P, F

class Mod:
    def __init__(self, q, P, K, V):
        self.q = q; self.T = Torus(q); self.m = self.T.m; self.n = self.m // 4; self.iota = self.T.iota
        self.par = nu(P, q); self.V = V
        self.Va = np.array(V, dtype=np.int64)
        self.L = [self.T.L(gauss_prime(p)) for p in P]
        self.units_s = {(s * self.iota) % self.m: s for s in range(4)}
    def collisions(self, ell):
        d = (self.Va @ np.array(ell, dtype=np.int64)) % self.m
        out = []
        for i, x in enumerate(d):
            if x % self.n == 0: out.append((i, self.units_s[int(x)]))
        return out

def lift_levels(mod, ell, coll, load, B, lq, maxlev=12, branch=6):
    """DFS over lifts: at each level choose t in (Z/q)^r (exhaustive if q^r <= 3e5, else 1e5 random), trying the
    `branch` best t (fewest new collisions) that keep every (v,s) within budget.  Success when a level has no
    collision (deeper levels then have none).  Returns (added load dict, last colliding level) or None."""
    q = mod.q; r = len(ell); nodes = [0]
    def rec(cur, m_a, C, add, lev):
        if not C: return add, lev - 1
        nodes[0] += 1
        if lev > maxlev or nodes[0] > 300: return None
        n_a = m_a // 4
        rows = []
        for (i, s) in C:
            v = mod.V[i]; x = sum(a_ * b_ for a_, b_ in zip(v, cur))
            assert x % n_a == 0
            rows.append((i, s, v, x // n_a))
        rng = np.random.default_rng(lev + 7 * q)
        Vr = np.array([v for (_, _, v, _) in rows], dtype=np.int64) % q
        dr = np.array([d for (_, _, _, d) in rows], dtype=np.int64) % q
        if q ** r <= 300000:
            Tt = np.array(list(itertools.product(range(q), repeat=r)), dtype=np.int64)
        else:
            Tt = rng.integers(0, q, size=(100000, r), dtype=np.int64)
        hit = ((dr[None, :] + 4 * (Tt @ Vr.T)) % q == 0)
        # budget feasibility of each t
        idx_rows = [(rows[k][0], rows[k][1]) for k in range(len(rows))]
        slack = np.array([B[i, s] - load[i, s] - add.get((i, s), 0.0) - lq for (i, s) in idx_rows])
        feas = ~(hit & (slack[None, :] < -1e-9)).any(axis=1)
        # also two rows may share (i,s)? each v collides with one unit, rows distinct (i) -> fine
        cand = np.where(feas)[0]
        if len(cand) == 0: return None
        cnt = hit[cand].sum(axis=1)
        tried = set()
        for k in cand[np.argsort(cnt, kind='stable')]:
            key = tuple(np.where(hit[k])[0])
            if key in tried: continue
            tried.add(key)
            if len(tried) > branch: break
            newC = [idx_rows[j] for j in key]
            add2 = dict(add)
            for (i, s) in newC: add2[(i, s)] = add2.get((i, s), 0.0) + lq
            cur2 = [a_ + m_a * int(b_) for a_, b_ in zip(cur, Tt[k])]
            res = rec(cur2, m_a * q, newC, add2, lev + 1)
            if res is not None: return res
        return None
    return rec(list(ell), mod.m, coll, {}, 2)

def gen_good(mods, B, rng, maxtries=400):
    nv = B.shape[0]
    load = np.zeros_like(B); system = {}
    order = list(mods)
    rng.shuffle(order)
    for mod in order:
        lq = math.log(mod.q)
        for _ in range(maxtries):
            ell = [rng.randrange(mod.par[j], mod.m, 2) for j in range(len(mod.par))]
            coll = mod.collisions(ell)
            ok = all(load[i, s] + lq <= B[i, s] + 1e-9 for i, s in coll)
            if not ok: continue
            nl = load.copy()
            for i, s in coll: nl[i, s] += lq
            lifted = lift_levels(mod, ell, coll, nl, B, lq)
            if lifted is None: continue
            add, levs = lifted
            for (i, s), x in add.items(): nl[i, s] += x
            load = nl; system[mod.q] = (ell, levs); break
        else:
            return None
    return system, load

def gen_good2(mods, B, rng, per_q=8, sample=6000, node_cap=400):
    """DFS over moduli in increasing q.  Candidate pool: all parity vectors if <= 2e5, else a random sample;
    vectorised level-1 feasibility; candidates ordered by (max load ratio after adding, random); each tried
    candidate must pass the exact all-level lift search.  Global node cap; returns None on failure."""
    mods = sorted(mods, key=lambda md: md.q)
    nv = B.shape[0]; ar = np.arange(nv)[None, :]
    nodes = [0]
    def pool(md):
        r = len(md.par); half = md.m // 2
        if half ** r <= 200000:
            E = np.array(list(itertools.product(*[range(md.par[j], md.m, 2) for j in range(r)])), dtype=np.int64)
        else:
            E = np.stack([2 * np.array([rng.randrange(half) for _ in range(sample)]) + md.par[j] for j in range(r)], axis=1).astype(np.int64)
        D = (E @ md.Va.T) % md.m
        S = np.full(D.shape, -1, dtype=np.int64)
        for val, s_ in md.units_s.items(): S[D == val] = s_
        return E, S
    depth = [0]
    def rec(k, load, system):
        depth[0] = max(depth[0], k)
        gen_good2.depth = depth[0]
        if k == len(mods): return system, load
        nodes[0] += 1
        if nodes[0] > node_cap: return None
        md = mods[k]; lq = math.log(md.q)
        E, S = pool(md)
        coll = S >= 0
        Sg = np.where(coll, S, 0)
        newload = load[ar, Sg] + lq
        ratio = np.where(coll, newload / B[ar, Sg], 0.0)
        worst = ratio.max(axis=1)
        feas = np.where(worst <= 1 + 1e-9)[0]
        if len(feas) == 0: return None
        noise = np.array([rng.random() for _ in range(len(feas))])
        order = feas[np.lexsort((noise, np.round(worst[feas], 6)))]
        # keep diversity: take the best per_q among a random half of the 4*per_q best
        top = list(order[:4 * per_q]); rng.shuffle(top); top = top[:per_q]
        for c in top:
            ell = [int(x) for x in E[c]]
            cl = [(i, int(S[c, i])) for i in np.where(coll[c])[0]]
            nl = load.copy()
            for i, s_ in cl: nl[i, s_] += lq
            lifted = lift_levels(md, ell, cl, nl, B, lq)
            if lifted is None: continue
            add, levs = lifted
            for (i, s_), x in add.items(): nl[i, s_] += x
            system[md.q] = (ell, levs)
            res = rec(k + 1, nl, system)
            if res is not None: return res
            del system[md.q]
            if nodes[0] > node_cap: return None
        return None
    return rec(0, np.zeros_like(B), {})

def sample_from_pickle(pk, mods, B, rng):
    # exact sampling: signature tuple with prob proportional to multiplicity, then a uniform vector in each class
    good = pk['good']; weights = [m for _, m in good]
    ch, _ = rng.choices(good, weights=weights)[0]
    load = np.zeros_like(B); system = {}
    for i, md in enumerate(mods):
        t = pk['tables_full'][i]
        idx = np.where(t['inv'] == ch[i])[0]
        ell = [int(x) for x in t['ells'][rng.choice(idx)]]
        coll = md.collisions(ell); lq = math.log(md.q)
        for a, s_ in coll: load[a, s_] += lq
        system[md.q] = (ell, 1)
    assert (load <= B + 1e-9).all()
    return system, load

def coherent(system, mods, Zs, Lz, qset):
    """system-level coherence with some zeta-tuple (norm <= H) at every q in qset (level 1, mod n_q)."""
    r = len(next(iter(system.values()))[0])
    qs = sorted(qset, key=lambda q: mods_by_q[q].n)
    buckets = {}
    for q in qs:
        n = mods_by_q[q].n; b = {}
        for zi, x in enumerate(Lz[q]):
            if x >= 0: b.setdefault(x % n, set()).add(zi)
        buckets[q] = b
    allz = set(range(len(Zs)))
    def dfs(k, cands):
        if k == len(qs): return True
        q = qs[k]; n = mods_by_q[q].n; c = [x % n for x in system[q][0]]
        for u in range(1, n + 1):
            if n > 1 and math.gcd(u, n) != 1: continue
            if n == 1 and u > 1: break
            new = []
            for j in range(r):
                s = cands[j] & buckets[q].get((u * c[j]) % n, set())
                if not s: break
                new.append(s)
            else:
                if dfs(k + 1, new): return True
        return False
    return dfs(0, [allz] * r)

if __name__ == '__main__':
    K, r, tag, nsys, H = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
    seed = int(sys.argv[6]) if len(sys.argv) > 6 else 1
    P, F = load_set(K, r, tag)
    V = monos(r, K)
    B = np.array([[budget(P, v, s) for s in range(4)] for v in V])
    mods = [Mod(q, P, K, V) for q in F]
    mods_by_q = {md.q: md for md in mods}
    Zs = [(a, b) for a in range(1, math.isqrt(H) + 1) for b in range(0, math.isqrt(H) + 1) if a * a + b * b <= H]
    Lz = {md.q: [md.T.L(z) if (z[0] ** 2 + z[1] ** 2) % md.q else -1 for z in Zs] for md in mods}
    rng = random.Random(seed)
    # genuine system check (sanity): genuine signature must be good
    gen = {md.q: (md.L, 0) for md in mods}
    print(f"K={K} r={r} {tag} P={P} F={F} H={H} (#zeta={len(Zs)})")
    print(f"  genuine coherent (sanity): all F: {coherent(gen, mods, Zs, Lz, F)}")
    F7 = [q for q in F if q % 8 == 7]
    stats = dict(good=0, fail=0, coh_all=0, coh_7=0, gen_sig_all=0)
    t0 = time.time(); examples = []
    import os
    pkf = f'forced_sig_{K}_{r}_{tag}.pkl'
    pk = None
    if os.path.exists(pkf):
        from forced_sig import sig_table
        pk = pickle.load(open(pkf, 'rb'))
        pk['tables_full'] = [sig_table(P, K, q) for q in F]
        print("  mode: exact weighted sampling from all good level-1 systems, then lifts through all levels")
    else:
        print("  mode: randomized greedy generation (budget-respecting), lifts through all levels")
    stats['lift_fail'] = 0
    while stats['good'] < nsys:
        if pk is not None:
            system, load = sample_from_pickle(pk, mods, B, rng)
            okl = True
            for md in mods:
                ell = system[md.q][0]; coll = md.collisions(ell)
                lifted = lift_levels(md, ell, coll, load, B, math.log(md.q))
                if lifted is None: okl = False; break
                add, levs = lifted
                for (i, s_), x in add.items(): load[i, s_] += x
                system[md.q] = (ell, levs)
            if not okl: stats['lift_fail'] += 1; continue
            g = (system, load)
        else:
            g = gen_good2(mods, B, rng)
        if g is None: stats['fail'] += 1; continue
        system, load = g
        stats['good'] += 1
        ca = coherent(system, mods, Zs, Lz, F); c7 = coherent(system, mods, Zs, Lz, F7)
        stats['coh_all'] += ca; stats['coh_7'] += c7
        gs = all(sorted(md.collisions(system[md.q][0])) == sorted(md.collisions(md.L)) for md in mods)
        stats['gen_sig_all'] += gs
        if not ca and len(examples) < 3: examples.append({q: system[q][0] for q in system})
    print(flush=True); print(f"  generated good systems (all levels): {stats['good']} (failed restarts {stats['fail']}, level-1-good samples with no good lift {stats['lift_fail']}); "
          f"coherent on all forced q: {stats['coh_all']}; coherent on forced Q_7 q {F7}: {stats['coh_7']}; "
          f"genuine signature everywhere: {stats['gen_sig_all']}  [{time.time()-t0:.0f}s]")
    for ex in examples[:1]:
        print("  example non-coherent good system (level-1 labels):", ex)
    pickle.dump(dict(P=P, F=F, K=K, examples=examples, stats=stats), open(f'forced_coh_{K}_{r}_{tag}_H{H}.pkl', 'wb'))
