# kill_search: own realisation solver for a type (rows = sign vectors). Base form at row 0, free orientation,
# square system: moments of independent sets + 2 gauges + random slices; LM from random starts.
import sys, pickle, itertools, numpy as np
from scipy.optimize import least_squares
def realise(rows, c, starts, rng, require_orient=True):
    rows = np.array(rows); M, D = rows.shape
    K = (D + 2*c)//2; L = D-1-K
    Fs = np.array([(rows[i] != rows[0]).astype(float) for i in range(1, M)])
    # independent subset
    basis = []
    for F in Fs:
        if np.linalg.matrix_rank(np.array(basis + [F])) > len(basis): basis.append(F)
    B = np.array(basis); r = len(B)
    nfree = 2*D - r*L - 2
    s0 = rows[0].astype(float)
    hits = []
    for st in range(starts):
        x0 = rng.normal(0, 0.7, D); y0 = s0*np.abs(rng.normal(0.8, 0.4, D))
        z0 = np.concatenate([x0, y0])
        idx = rng.choice(2*D, size=max(nfree, 0), replace=False); vals = z0[idx]
        def res(z):
            q = z[:D] + 1j*z[D:]
            out = [B @ np.imag(q**m) for m in range(1, L+1)]
            out.append([np.sum(z[:D]), np.sum(np.abs(q)**2) - D])
            out.append(z[idx] - vals)
            return np.concatenate([np.ravel(o) for o in out])
        meth = 'lm' if r*L + 2 + len(idx) >= 2*D else 'trf'
        try:
            s = least_squares(res, z0, method=meth, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
        except Exception:
            continue
        if np.linalg.norm(s.fun) > 1e-11: continue
        q = s.x[:D] + 1j*s.x[D:]
        allr = np.concatenate([q, q.conj()])
        sep = min(abs(allr[i]-allr[j]) for i in range(2*D) for j in range(i))
        if sep < 1e-3: continue
        orient_ok = bool(np.all(np.sign(q.imag) == s0))
        if require_orient and not orient_ok: continue
        hits.append((sep, q, orient_ok))
        if len(hits) >= 1: break
    return hits, r, nfree
if __name__ == "__main__":
    fn, c, starts = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    types = pickle.load(open(fn, 'rb'))
    rng = np.random.default_rng(11)
    nreal = 0; nreal_any = 0
    for t, rows in enumerate(types):
        h1, r, nf = realise(rows, c, starts, rng, True)
        h2 = h1 if h1 else realise(rows, c, starts//2, rng, False)[0]
        nreal += bool(h1); nreal_any += bool(h2)
        print("type", t, "rank", r, "free dims", nf, "realised(orient)", bool(h1), "realised(any orient)", bool(h2),
              ("sep %.3g" % h2[0][0]) if h2 else "", flush=True)
    print("realised with prescribed orientation:", nreal, "of", len(types), "; with any orientation:", nreal_any)
