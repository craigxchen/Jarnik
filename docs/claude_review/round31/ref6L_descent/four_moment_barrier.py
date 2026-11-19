"""K-moment Im system on a real line with a collision barrier: minimize |F|^2 + eps * sum 1/|pairwise distance|^2
(normalized by scale), then report the smallest K-moment residual reached among strict points.
Evidence only."""
import numpy as np, sys
from itertools import combinations
from scipy.optimize import least_squares
edges = list(combinations(range(5), 2))
Inc = np.array([[int(i in e) for e in edges] for i in range(5)], dtype=float)
E = Inc[:4] - Inc[4]
K = int(sys.argv[1]); trials = int(sys.argv[2]); eps = float(sys.argv[3])
rng = np.random.default_rng(11)
pairs = list(combinations(range(20), 2))
def unpack(p):
    x = np.concatenate([[0.0], p[:9]]); y = np.concatenate([[1.0], p[9:]])
    return x + 1j*y
def F(p):
    r = unpack(p)
    sc = max(1.0, np.max(np.abs(r)))
    return np.concatenate([E @ np.imag((r/sc)**k) for k in range(1, K+1)])
def Fb(p):
    r = unpack(p); pts = np.concatenate([r, np.conj(r)])
    sc = max(1.0, np.max(np.abs(r)))
    bar = np.array([eps*sc/abs(pts[a]-pts[b]) for a, b in pairs])
    return np.concatenate([F(p), bar])
out = []
for t in range(trials):
    p0 = rng.normal(size=18)*2
    s = least_squares(Fb, p0, max_nfev=3000)
    r = unpack(s.x); pts = np.concatenate([r, np.conj(r)])
    sc = max(1.0, np.max(np.abs(r)))
    st = min(abs(pts[a]-pts[b]) for a, b in pairs)/sc
    out.append((np.max(np.abs(F(s.x))), st))
out.sort()
print('K=%d eps=%g: best (scaled residual, strictness):' % (K, eps), [(f'{a:.1e}', f'{b:.1e}') for a, b in out[:8]])
