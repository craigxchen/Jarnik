"""Five-point subconfigurations of the Part V cluster: reduced exponent delta_S, McKinnon-Roth exponent
alpha_S = 4/(1+2 delta_S), the 15 cut weights (balance), H_Sigma = sum of reduced-chord degrees.
Usage: python3 five.py m L [nsample]
"""
import sys, pickle, itertools, random, json
from fractions import Fraction as Fr

def parity(x):
    return bin(x).count('1') & 1

m = int(sys.argv[1]); L = int(sys.argv[2])
ns = int(sys.argv[3]) if len(sys.argv) > 3 else None
S = pickle.load(open(f'cluster_m{m}_L{L}.pkl', 'rb'))
d, Q, D, forms, H, sp, orient, bigdeg, vz, vN, ddiff, gd = (S[k] for k in
    ['d', 'Q', 'D', 'forms', 'H', 'sp', 'orient', 'bigdeg', 'vz', 'vN', 'ddiff', 'gd'])
from f24 import deg4
Hs = {a: set(H[a]) for a in forms}
layers = []   # (weight, {a: marked?}) as list of (weight, frozenset of forms marked)
for v in range(2, Q, 2):
    layers.append((2 * bigdeg[v], frozenset(a for a in forms if v in Hs[a])))
for f in orient:
    for l in range(1, vN[f] + 1):
        layers.append((2 * deg4(f), frozenset(a for a in forms if vz[(f, a)] >= l)))
subsets = list(itertools.combinations(forms, 5))
if ns and ns < len(subsets):
    random.seed(2)
    subsets = random.sample(subsets, ns)
def dd(a, b):
    return ddiff[(a, b)] if (a, b) in ddiff else ddiff[(b, a)]
def gg(a, b):
    return gd[(a, b)] if (a, b) in gd else gd[(b, a)]
stats = []
liouv_ok = True
for Sx in subsets:
    cut = {}
    tot = 0
    for w, mk in layers:
        A = frozenset(a for a in Sx if a in mk)
        if len(A) in (0, 5):
            continue
        key = min(A, frozenset(Sx) - A, key=lambda s: sorted(s))
        cut[key] = cut.get(key, 0) + w
        tot += w
    Dp = tot // 2            # reduced point degree D'' (sum of layers = deg N'' = 2 D'')
    gS = D - Dp
    pr = list(itertools.combinations(Sx, 2))
    maxd = max(dd(a, b) for a, b in pr) - gS
    HS = sum(dd(a, b) - gg(a, b) for a, b in pr)
    # Liouville-type identity: deg(zeta_a - zeta_b) = D'' - d_ab/2 + deg c_ab
    for a, b in pr:
        dab = sum(w for A, w in cut.items() if (a in A) != (b in A))
        if 2 * (dd(a, b) - gS) != 2 * Dp - dab + 2 * (dd(a, b) - gg(a, b)):
            liouv_ok = False
    ws = [cut.get(k, 0) for k in cut]
    ncut = len(cut)
    delta = Fr(1, 2) - Fr(maxd, Dp)
    stats.append((delta, ncut, min(ws) if ncut == 15 else 0, max(ws), HS, Dp))
deltas = [s[0] for s in stats]
out = dict(m=m, L=L, D=D, subsets=len(subsets), identity_ok=liouv_ok)
out['max_delta_S'] = float(max(deltas))
out['min_delta_S'] = float(min(deltas))
out['frac_delta_ge_1_40'] = sum(1 for x in deltas if x >= Fr(1, 40)) / len(deltas)
out['frac_delta_ge_1_31'] = sum(1 for x in deltas if x >= Fr(1, 31)) / len(deltas)
out['any_delta_gt_1_30'] = any(x > Fr(1, 30) for x in deltas)
best = max(stats, key=lambda s: s[0])
out['best'] = dict(delta=str(best[0]), delta_float=float(best[0]), alpha=float(4 / (1 + 2 * best[0])),
                   ncuts=best[1], minw=best[2], maxw=best[3], H_Sigma=best[4], Dred=best[5],
                   mean_w=2 * best[5] / 15)
full = [s for s in stats if s[1] == 15]
out['frac_all_15_cuts_present'] = len(full) / len(stats)
if full:
    out['exactly_balanced_frac'] = sum(1 for s in full if s[2] == s[3]) / len(stats)
    out['max_balance_ratio'] = max(s[3] / s[2] for s in full)
    out['max_HSigma_over_meanw'] = max(s[4] / (2 * s[5] / 15) for s in full)
print(json.dumps(out, indent=1))
json.dump(out, open(f'five_m{m}_L{L}.json', 'w'), indent=1)
from collections import Counter
for nc in sorted(set(s[1] for s in stats)):
    ds = [s[0] for s in stats if s[1] == nc]
    print("ncuts", nc, "count", len(ds), "delta min %.5f max %.5f" % (float(min(ds)), float(max(ds))),
          "frac>=1/40 %.3f" % (sum(1 for x in ds if x >= Fr(1, 40)) / len(ds)),
          "frac>=1/31 %.3f" % (sum(1 for x in ds if x >= Fr(1, 31)) / len(ds)),
          "frac>1/30 %.3f" % (sum(1 for x in ds if x > Fr(1, 30)) / len(ds)))
    hs = [s[4] / (2 * s[5] / 15) for s in stats if s[1] == nc]
    print("    H_Sigma/mean_w max %.4f ; balance ratio max %s" % (max(hs), max((s[3] / s[2]) if s[2] else float('inf') for s in stats if s[1] == nc)))
def rank2(vs):
    basis = []
    for x in vs:
        for b in basis:
            x = min(x, x ^ b)
        if x:
            basis.append(x)
    return len(basis)
print("--- by rank of the 5 forms ---")
rk = [rank2(Sx) for Sx in subsets]
for r_ in sorted(set(rk)):
    idx = [i for i in range(len(subsets)) if rk[i] == r_]
    ds = [stats[i][0] for i in idx]
    hs = [stats[i][4] / (2 * stats[i][5] / 15) for i in idx]
    br = [stats[i][3] / stats[i][2] if stats[i][2] else float('inf') for i in idx]
    print("rank", r_, "count", len(idx), "delta in [%.5f, %.5f]" % (float(min(ds)), float(max(ds))),
          "frac>=1/40 %.3f frac>=1/31 %.3f frac>1/30 %.3f" % (
              sum(1 for x in ds if x >= Fr(1, 40)) / len(ds), sum(1 for x in ds if x >= Fr(1, 31)) / len(ds),
              sum(1 for x in ds if x > Fr(1, 30)) / len(ds)),
          "H_Sigma/mean_w max %.4f" % max(hs), "balance ratio max %.4f" % max(br),
          "alpha range [%.4f, %.4f]" % (float(4 / (1 + 2 * max(ds))), float(4 / (1 + 2 * min(ds)))))
