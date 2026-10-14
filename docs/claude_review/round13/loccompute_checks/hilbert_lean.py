"""Memory-lean isotypic Hilbert function (mod PR) for large models (Syl(16)).
For component s: matrix M^T (selected point-orbit reps x monomial-orbit sums of degree <= maxdeg, degree-ordered);
rref; pivot columns give the cumulative rank profile.  Using a SUBSET of point reps only lowers ranks:
rank_p(subset) <= rank_p(full) <= rank_Q(full), so N_{<=i} - (cum rank) stays an upper bound on dim_Q J(A)_{<=i}."""
import numpy as np, time, random
from flint import nmod_mat
from loclib import PR

def hilbert_lean(cr, A, maxdeg, factor=1.3, seed=1, verbose=True):
    ring = cr.ring; q = cr.q; d = ring.d
    A = np.asarray(A, dtype=np.int8)
    seen = set(); reps = []; rsize = []
    for k in A.tolist():
        k = tuple(k)
        if k in seen: continue
        orb = {k}; cur = k
        for t in range(1, q):
            cur = cur[-1:] + cur[:-1]; orb.add(cur)
        seen |= orb; reps.append(k); rsize.append(len(orb))
    R = np.array(reps, dtype=np.int64); rsize = np.array(rsize)
    nrow_total = 1 + sum(cr.orb[i].shape[0] for i in range(1, maxdeg + 1))
    bounds = [1] + [1 + sum(cr.orb[j].shape[0] for j in range(1, i + 1)) for i in range(1, maxdeg + 1)]
    rnd = random.Random(seed)
    H = [0] * (maxdeg + 1); percomp = []
    for s in range(q):
        t0 = time.time()
        okcols = np.nonzero(np.array([pow(cr.w, (s * int(r)) % q, PR) == 1 for r in rsize]))[0]
        ncol = min(len(okcols), int(factor * nrow_total) + 50)
        sel = np.array(sorted(rnd.sample(list(okcols), ncol)))
        Pts = R[sel]
        M = nmod_mat(ncol, nrow_total, PR)
        # constant monomial
        if s == 0:
            for a in range(ncol): M[a, 0] = 1
        col = 1
        for i in range(1, maxdeg + 1):
            O = cr.orb[i]; Mons = ring.mons[i]
            for c0 in range(0, O.shape[0], 256):
                blk = O[c0:c0 + 256]                      # (b, q)
                acc = np.zeros((blk.shape[0], ncol), dtype=np.int64)
                for t in range(q):
                    E = Mons[blk[:, t]]                    # (b, d)
                    vals = np.ones((blk.shape[0], ncol), dtype=np.int64)
                    for x in range(d):
                        e = E[:, x][:, None]; kx = Pts[:, x][None, :]
                        vals *= np.where(e == 0, 1, np.where(e == 1, kx, kx * kx))
                    acc = (acc + (vals % PR) * int(cr.Wbwd[s][t])) % PR
                for b in range(acc.shape[0]):
                    row = acc[b]
                    for a in np.nonzero(row)[0]:
                        M[int(a), col + b] = int(row[a])
                col += acc.shape[0]
        Rr, rk = M.rref()
        piv = []; j = 0
        for r in range(rk):
            while int(Rr[r, j]) == 0: j += 1
            piv.append(j); j += 1
        piv = np.array(piv)
        cum = [int((piv < bounds[i]).sum()) for i in range(maxdeg + 1)]
        hs = [cum[0]] + [cum[i] - cum[i - 1] for i in range(1, maxdeg + 1)]
        percomp.append(hs)
        for i in range(maxdeg + 1): H[i] += hs[i]
        if verbose: print(f"    comp {s}: cols {ncol}/{len(okcols)} rank {rk} H_s={hs} [{time.time()-t0:.0f}s]", flush=True)
    return H, percomp
