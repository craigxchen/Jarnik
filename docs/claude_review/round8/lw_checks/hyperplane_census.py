"""Census of heavy hyperplanes (equivalently low-support characters) for several profiles.

For a list of +-1 type vectors s_j in Z^M (equal weights), a primitive zero-sum c defines the
hyperplane c^perp; load(c) = #{j : <c,s_j> = 0}, s(c) = #types - load(c).  Every hyperplane with
load >= M-2 independent types is spanned by M-2 of its types plus the all-ones vector, so
enumerating all (M-2)-subsets finds every c with s(c) <= #types - (M-2).  Float SVD is used to
find candidates; each reported c is recomputed and its support re-verified in exact arithmetic.

Profiles: 'dcore' (Paley skew double core), 'rand2' / 'rand3' (2 or 3 random row permutations of
the Paley core, nonconstant columns), 'single' (single Paley core, for reference).
Usage: python3 hyperplane_census.py q profile seed slack
  prints the number of distinct primitive c (up to sign) with s(c) = smin, smin+1, ..., smin+slack.
"""
import itertools, sys, time
import numpy as np
from dcore import paley_skew, double_core_types, random_permuted_bases, image, nullvec_int

q = int(sys.argv[1]); prof = sys.argv[2]; seed = int(sys.argv[3]); slack = int(sys.argv[4])
K = paley_skew(q)
M = q + 1
if prof == 'dcore':
    f1, f2 = double_core_types(K)
    types = f1 + f2
elif prof.startswith('rand'):
    t = int(prof[4:])
    types = sum(random_permuted_bases(K, t, seed), [])
elif prof == 'single':
    types = random_permuted_bases(K, 1, seed)[0]
n = len(types)
T = np.array(types, dtype=float)
ones = np.ones(M)
t0 = time.time()
cands = {}
B = 20000
it = itertools.combinations(range(n), M - 2)
total = 0
best = -1
while True:
    chunk = list(itertools.islice(it, B))
    if not chunk:
        break
    total += len(chunk)
    chunk = np.array(chunk)
    A = np.concatenate([T[chunk], np.broadcast_to(ones, (len(chunk), 1, M))], axis=1)
    U, S, Vt = np.linalg.svd(A)
    ok = S[:, -1] > 1e-8
    c = Vt[:, -1, :]
    loads = (np.abs(c @ T.T) < 1e-7).sum(axis=1)
    loads[~ok] = -1
    best = max(best, loads.max())
    thr = best - slack
    for j in np.nonzero(loads >= thr)[0]:
        key = tuple(np.round(c[j] / np.max(np.abs(c[j])), 6))
        key2 = tuple(-x for x in key)
        if key in cands or key2 in cands:
            continue
        cands[key] = tuple(chunk[j])
# exact verification
exact = {}
for key, sub in cands.items():
    rows = [types[j] for j in sub] + [[1] * M]
    cv = nullvec_int(rows, M)
    if cv is None:
        continue
    s = sum(1 for x in image(cv, types) if x != 0)
    exact[tuple(cv)] = s
smin = min(exact.values())
hist = {}
for cv, s in exact.items():
    hist[s] = hist.get(s, 0) + 1
lines = [f"q={q} M={M} profile={prof} seed={seed} types={n}: enumerated {total} hyperplanes in {time.time()-t0:.0f}s",
         f"  min s(c) = {smin} ; s(c)/#types = {smin/n:.4f} ; (M/2)/#types = {(M/2)/n:.4f}"]
for s in sorted(hist):
    if s <= smin + slack:
        ex = [cv for cv, ss in exact.items() if ss == s][:2]
        lines.append(f"  s = {s}: {hist[s]} primitive characters (up to sign); e.g. {[list(e) for e in ex]}")
print('\n'.join(lines))
with open(f'hyperplane_census_{prof}_q{q}_s{seed}_out.txt', 'w') as fh:
    fh.write('\n'.join(lines) + '\n')
