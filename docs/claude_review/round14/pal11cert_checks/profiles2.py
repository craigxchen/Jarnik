"""Feasible (sigma,t) profiles of 8 generators with exact attainable (sigma,t) pairs for sigma <= 8 (typeflats.py),
t >= f(sigma) for sigma = 9..12 (hadflats2.py), sum sigma <= 26, sum t <= 44.  Output: all profiles and the set of
(sigma,t) types occurring."""
from collections import Counter
att = {2: [6], 3: [9], 4: [6, 8, 9, 10, 11], 5: [7, 8, 9, 10, 11], 6: [2, 4, 6, 7, 8, 9, 10, 11],
       7: [5, 6, 7, 8, 9, 10, 11], 8: [4, 5, 6, 7, 8, 9, 10, 11],
       9: list(range(3, 12)), 10: list(range(4, 12)), 11: list(range(4, 12)), 12: list(range(1, 12))}
types = sorted((s, t) for s in att for t in att[s])
full = []
def rec(start, k, ss, tt, cur):
    if k == 8:
        full.append(tuple(cur)); return
    for i in range(start, len(types)):
        s, t = types[i]
        if ss + s + 2 * (7 - k) > 26: continue
        if tt + t + 1 * (7 - k) > 44: continue
        cur.append((s, t)); rec(i, k + 1, ss + s, tt + t, cur); cur.pop()
rec(0, 0, 0, 0, [])
print('feasible profiles:', len(full))
occ = Counter()
for p in full:
    for ty in set(p): occ[ty] += 1
print('types occurring:', sorted(occ))
for p in full: print(p)
