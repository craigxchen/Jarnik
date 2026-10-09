"""Re-tally non-admissible trees of anchor_check.py with an exact (brute-force) two-chain test
instead of the greedy chains_ok (tally only; not used in any proof)."""
import sys, itertools, random
sys.path.insert(0, '../L7')
import importlib.util
src = open('../L7/anchor_check.py').read().split('random.seed(1)')[0]
ns = {}; exec(src, ns)
def exact_ok(pos):
    n = len(pos)
    for mask in range(1 << n):
        c1 = sorted([pos[k] for k in range(n) if mask >> k & 1], key=len)
        c2 = sorted([pos[k] for k in range(n) if not mask >> k & 1], key=len)
        if all(x <= y for x, y in zip(c1, c1[1:])) and all(x <= y for x, y in zip(c2, c2[1:])):
            if not (c1 and c2 and (c1[-1] & c2[-1])): return True
    return False
random.seed(1)
greedy = exact = disagree = 0
for trial in range(400):
    pts = [(1, 0), (0, 1), (1, 1)]
    while len(pts) < 7:
        q = ns['prim']((random.randint(-60, 60), random.randint(1, 60)))
        if all(ns['det'](q, r) != 0 for r in pts): pts.append(q)
    w = pts
    allp = set()
    for i, j in itertools.combinations(range(7), 2): allp |= ns['primes_of'](ns['det'](w[i], w[j]))
    for p in allp:
        g = {(i, j): ns['vp'](ns['det'](w[i], w[j]), p) for i, j in itertools.combinations(range(7), 2)}
        l = ns['edges'](g, list(range(7)))
        pos = [T for T in ns['Ts'] if l[T] > 0]
        a = not ns['chains_ok'](pos); b = not exact_ok(pos)
        greedy += a; exact += b; disagree += (a != b)
print('non-admissible: greedy', greedy, 'exact', exact, 'disagreements', disagree)
