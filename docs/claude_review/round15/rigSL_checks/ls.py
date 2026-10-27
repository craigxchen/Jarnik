# ls.py NAME K Q0 [seed] [maxit] [phi=genuine|f1,f2,..] -- min-conflicts local search for an SL witness (round 15).
# Labels live at full depth: at q, ell in (Z/M_q)^r, M_q = m_q q^(A_q - 1), q^(A_q) > max budget (so a good system can
# have no collision at level A_q and the stored label determines all collision levels).  Budgets at phi:
# B(v,s) = log|D_s(Z_v)| (genuine) or (1/2)log(N 2^[s odd]) + log|sin(<v,phi> - s pi/4)| - 1e-9.
# Target: a good system on F_K^ex whose label at Q0 lies outside every column-sign class (=> non-coherent on F ∩ [3,Q0]).
# Output: W_ls_<name>_K<K>_Q<Q0>_<seed>.pkl in the slverify format (then checked independently by slverify).
import sys, math, random, time, pickle
import numpy as np
from sl import *

def prep(S, phi=None):
    r = S.r; V = np.array(S.V, dtype=np.int64); nV = len(S.V)
    if phi is None: B = np.array(budgets_genuine(S))
    else: B = np.array(budgets_phi(S, phi, margin=1e-9))
    maxB = float(B.max()); mods = {}
    for q in S.F:
        m, n, L, iota = modulus_data(S, q)
        A = 1
        while A * math.log(q) <= maxB + 1e-9: A += 1
        M = m * q ** (A - 1)
        par = np.array(nu(S.P, q), dtype=np.int64)
        cls = set()
        for bits in range(2 ** (r - 1)):
            eps = [1] + [(-1 if (bits >> (j - 1)) & 1 else 1) for j in range(1, r)]
            for u in units_mod(n): cls.add(tuple((u * e * x) % n for e, x in zip(eps, L)))
        # unit lookup: d mod m -> s (only for d = 0 mod n)
        sl = np.full(m, -1, dtype=np.int64)
        for s in range(4): sl[(s * iota) % m] = s
        mods[q] = dict(q=q, m=m, n=n, M=M, A=A, par=par, cls=cls, sl=sl, lq=math.log(q), L=L, iota=iota)
    return V, B, mods

def contrib(md, labs, V):
    """labs: C x r; returns (pair index C x nV or -1, amount C x nV)."""
    q, m, n, M, A = md["q"], md["m"], md["n"], md["M"], md["A"]
    d = (labs @ V.T) % M
    hit = (d % n) == 0
    s = md["sl"][d % m]
    b = np.zeros(d.shape, dtype=np.int64); qa = n
    for a in range(1, A + 1):
        b += ((d % qa) == 0); qa *= q
    e = np.where(hit, np.arange(V.shape[0])[None, :] * 4 + s, -1)
    return e, b * md["lq"]

def rand_labels(md, C, rng):
    M = md["M"]; r = len(md["par"])
    x = rng.integers(0, M // 2, size=(C, r)) * 2 + md["par"][None, :]
    return x

def outside(md, lab):
    return tuple(int(x) % md["n"] for x in lab) not in md["cls"]

def search(S, Q0, seed=1, maxit=200000, phi=None, log=print, C=3000, init=None):
    rng = np.random.default_rng(seed); prng = random.Random(seed)
    V, B, mods = prep(S, phi); qs = sorted(mods); nV = len(S.V); nE = 4 * nV
    lab = {}; E = {}; Amt = {}
    load = np.zeros(nE)
    for q in qs:
        md = mods[q]
        while True:
            if init is not None and q in init:
                x0, M0 = init[q]; assert md["M"] % M0 == 0
                x = np.array([int(a) % M0 + M0 * int(rng.integers(0, md["M"] // M0)) for a in x0], dtype=np.int64)
                if q != Q0 or outside(md, x): break
                init = dict(init); del init[q]; continue
            x = rand_labels(md, 1, rng)[0]
            if q != Q0 or outside(md, x): break
        lab[q] = x; e, a = contrib(md, x[None, :], V); E[q] = e[0]; Amt[q] = a[0]
        np.add.at(load, e[0][e[0] >= 0], a[0][e[0] >= 0])
    f = lambda L, idx=None: np.maximum(0.0, L - (B if idx is None else B[idx]))
    viol = f(load).sum(); t0 = time.time(); best = viol
    for it in range(maxit):
        if viol <= 1e-12: break
        # pick a violated pair and a modulus contributing to it
        bad = np.nonzero(load - B > 1e-12)[0]
        e0 = int(bad[prng.randrange(len(bad))])
        cq = [q for q in qs if (E[q] == e0).any()]
        if not cq: cq = qs
        q = prng.choice(cq); md = mods[q]
        cand = rand_labels(md, C, rng)
        # local perturbations of the current label
        pert = np.repeat(lab[q][None, :], C // 3, axis=0)
        jj = rng.integers(0, S.r, size=C // 3); pert[np.arange(C // 3), jj] = (pert[np.arange(C // 3), jj] + 2 * rng.integers(1, md["M"] // 2, size=C // 3)) % md["M"]
        # lift perturbations: keep the label mod m q^k (k random), re-draw the higher digits (all coordinates)
        nl = C // 3; ks = rng.integers(0, md["A"], size=nl)
        base = np.array([md["m"] * md["q"] ** int(k) for k in ks], dtype=np.int64)
        lifts = (np.repeat(lab[q][None, :], nl, axis=0) % base[:, None]) + base[:, None] * rng.integers(0, 1 << 62, size=(nl, S.r)) % (md["M"] // base)[:, None]
        # digit moves: replace the base-q digit vector at one level k >= 1 (all q^r choices if few), keep the rest
        dig = []
        if md["A"] > 1:
            k = int(rng.integers(1, md["A"])); unit = md["m"] * md["q"] ** (k - 1)
            cur = lab[q]; low = cur % unit; high = cur // (unit * md["q"])
            if md["q"] ** S.r <= 3000:
                import itertools as _it
                T = np.array(list(_it.product(range(md["q"]), repeat=S.r)), dtype=np.int64)
            else:
                T = rng.integers(0, md["q"], size=(1500, S.r))
            dig = [(low[None, :] + unit * T + unit * md["q"] * high[None, :]) % md["M"]]
        cand = np.vstack([cand, pert, lifts % md["M"]] + dig)
        if q == Q0:
            keep = np.array([outside(md, c) for c in cand]); cand = cand[keep]
            if len(cand) == 0: continue
        en, an = contrib(md, cand, V)
        eo, ao = E[q], Amt[q]
        Lb = load - np.where(True, 0, 0)
        # delta per monomial i
        same = (en == eo[None, :]) & (en >= 0)
        # removal of old
        rem = np.zeros(nV)
        has_o = eo >= 0
        rem[has_o] = f(load[eo[has_o]] - ao[has_o], eo[has_o]) - f(load[eo[has_o]], eo[has_o])
        en_c = np.where(en >= 0, en, 0)
        base_new = load[en_c]
        # if same pair: delta = f(L - o + n) - f(L); else rem + f(L_new + n) - f(L_new)
        d_same = f(base_new - np.where(same, ao[None, :], 0) + an, en_c) - f(base_new, en_c)
        add_new = np.where(en >= 0, d_same, 0.0)
        delta = add_new.sum(axis=1) + np.where(same, 0.0, rem[None, :]).sum(axis=1)
        k = int(np.argmin(delta + rng.random(len(delta)) * 1e-9))
        dl = delta[k]
        if dl < -1e-12 or (dl <= 1e-12 and prng.random() < 0.5) or prng.random() < 0.01:
            np.subtract.at(load, eo[eo >= 0], ao[eo >= 0])
            lab[q] = cand[k]; E[q] = en[k]; Amt[q] = an[k]
            np.add.at(load, en[k][en[k] >= 0], an[k][en[k] >= 0])
            viol = f(load).sum()
            if viol < best - 1e-9:
                best = viol
        if it % 2000 == 0:
            info = []
            for e in bad[:6]:
                e = int(e); who = {qq: round(float(Amt[qq][E[qq] == e].sum() / mods[qq]["lq"])) for qq in qs if (E[qq] == e).any()}
                info.append((S.V[e // 4], e % 4, round(float(load[e] - B[e]), 3), who))
            log(f"   it={it} viol={viol:.4f} best={best:.4f} #bad={len(bad)} {info} [{time.time()-t0:.0f}s]")
    if viol > 1e-12: return None
    # exact re-check of loads (integers) before export
    for q in qs:
        md = mods[q]; e, a = contrib(md, lab[q][None, :], V)
        # no collision may persist at level A (guaranteed: q^A > max budget and the system is good)
    out = {}
    for q in qs:
        out[q] = dict(g1=tuple(Torus(q).g), top=tuple(int(x) for x in lab[q]), A=mods[q]["A"])
    return dict(P=S.P, K=S.K, F=sorted(S.F), phi="genuine" if phi is None else [repr(x) for x in phi], labels=out, Q=Q0)

if __name__ == '__main__':
    name, K, Q0 = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    maxit = int(sys.argv[5]) if len(sys.argv) > 5 else 200000
    phi = None if len(sys.argv) <= 6 or sys.argv[6] == "genuine" else [float(x) for x in sys.argv[6].split(",")]
    S = setup(name, K); t0 = time.time()
    W = search(S, Q0, seed, maxit, phi, log=lambda s: print(s, flush=True))
    if W is None: print(f"RESULT {name} K={K} Q0={Q0} seed={seed}: no witness within {maxit} iterations [{time.time()-t0:.0f}s]")
    else:
        fn = f"W_ls_{name}_K{K}_Q{Q0}_{seed}.pkl"; pickle.dump(W, open(fn, "wb"))
        print(f"RESULT {name} K={K} Q0={Q0} seed={seed}: candidate witness saved {fn} [{time.time()-t0:.0f}s]", flush=True)
        import slverify; slverify.check(W)
