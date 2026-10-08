import itertools, sys
M = [1, 3, 1]
L = int(sys.argv[1]); B = int(sys.argv[2])
def mul(Q):
    out = [0]*(len(Q)+2)
    for i, a in enumerate(Q):
        for j, b in enumerate(M):
            out[i+j] += a*b
    return tuple(out)
cand = []
for Q in itertools.product(range(-B, B+1), repeat=L):
    if sum(Q) % 2: continue
    v = mul(Q)
    if 0 < sum(abs(x) for x in v) <= 12:
        cand.append(v)
cand = sorted(set(cand))
print("candidates", len(cand))
n = len(cand[0])
zero = tuple([0]*n)
# pairwise compatibility: l1 distance <= 12
import collections
ok = [[sum(abs(a-b) for a, b in zip(cand[i], cand[j])) <= 12 for j in range(len(cand))] for i in range(len(cand))]
found = collections.Counter(); ex = {}
def dfs(cur, lo, hi, start):
    k = len(cur) + 1
    if k >= 4:
        E = sum(h - l for l, h in zip(lo, hi))
        found[(k, E)] += 1
        ex.setdefault((k, E), list(cur))
    if k == 7: return
    for idx in range(start, len(cand)):
        if not all(ok[idx][c] for c in cur): continue
        v = cand[idx]
        nlo = tuple(min(a, b) for a, b in zip(lo, v)); nhi = tuple(max(a, b) for a, b in zip(hi, v))
        if sum(h - l for l, h in zip(nlo, nhi)) <= 12:
            dfs(cur + [idx], nlo, nhi, idx + 1)
dfs([], zero, zero, 0)
for key in sorted(found):
    print("rows=%d E=%d count=%d example=%s" % (key[0], key[1], found[key], [zero] + [cand[i] for i in ex[key]]))
