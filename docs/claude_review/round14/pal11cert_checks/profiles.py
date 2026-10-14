"""Enumerate the support/core-touch profiles {(sigma_i, t_i)} of 8 generators allowed by
sum sigma <= 26 (incidence), sum t <= 44 (core counts), t >= f(sigma) (hadflats2.py, exact), 2 <= sigma <= 12, t <= 11.
Only the multiset of sigma's and the 'exceptional' part (t <= 5) is reported."""
import itertools
from collections import Counter
f = {2: 6, 3: 9, 4: 6, 5: 7, 6: 2, 7: 5, 8: 4, 9: 3, 10: 4, 11: 4, 12: 1}
types = [(s, t) for s in range(2, 13) for t in range(f[s], 12)]
types.sort()
res = Counter()
def rec(start, k, ss, tt, cur):
    if k == 8:
        exc = tuple(sorted((s, t) for s, t in cur if t <= 5))
        sig = tuple(sorted(s for s, t in cur))
        res[(exc, sig)] += 1
        return
    for i in range(start, len(types)):
        s, t = types[i]
        if ss + s + 2 * (7 - k) > 26: continue
        if tt + t + 1 * (7 - k) > 44: continue
        # remaining 7-k generators each have t >= 1; sigma >= 2
        cur.append((s, t)); rec(i, k + 1, ss + s, tt + t, cur); cur.pop()
rec(0, 0, 0, 0, [])
ex = Counter()
for (exc, sig), c in res.items():
    ex[exc] += 1
for e, c in sorted(ex.items()):
    print('exceptional part', e, ':', c, 'support profiles')
allowed = Counter()
for (exc, sig), c in res.items():
    pass
# recompute keeping full multisets
full = []
def rec2(start, k, ss, tt, cur):
    if k == 8:
        full.append(tuple(cur)); return
    for i in range(start, len(types)):
        s, t = types[i]
        if ss + s + 2 * (7 - k) > 26: continue
        if tt + t + 1 * (7 - k) > 44: continue
        cur.append((s, t)); rec2(i, k + 1, ss + s, tt + t, cur); cur.pop()
rec2(0, 0, 0, 0, [])
nonexc = Counter()
for prof in full:
    for (s, t) in set(prof):
        if t >= 6: nonexc[(s, t)] += 1
print('number of feasible (sigma,t) profiles', len(full))
print('non-exceptional (sigma,t) types occurring:', sorted(nonexc))
for prof in full:
    if all(t >= 6 or True for s, t in prof):
        pass
print([p for p in full if sum(1 for s, t in p if t <= 5) == 1][:20])
