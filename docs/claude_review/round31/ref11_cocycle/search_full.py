# COMPLETE search (all bases, all divisors) of integral-cofactor configurations z = w + A/v on q_P, one-unit marks.
# A set S (|S|>=2) is a configuration at dilation 1/v iff v | T_ij and |A_j-A_i|^4 <= v^4 C^4 q for all pairs in S.
# (lambda = u/v with |u|>1 is never better, so v ranges over positive integers.)  Max clique per (base, v).
import sys, json, time
from collections import defaultdict
from sympy import divisors
from common import *
P = int(sys.argv[1]); r = int(sys.argv[2]); m = int(sys.argv[3])
Cs = [1, 4, 8]
ps, gs, q, s, marks = setup(P, (r, m)); n = 1 << s
print("P", P, "class", r, "mod", m, "s", s, "primes", ps, "q^{1/4}=%.4g" % q**0.25, flush=True)
divcache = {}
def okC(D2, v, C): return D2 * D2 <= v**4 * C**4 * q
best = {C: (1, None) for C in Cs}; count3 = {C: 0 for C in Cs}; sizes = {C: defaultdict(int) for C in Cs}
t0 = time.time()
def maxclique(nodes, adj):
    bestc = []
    def rec(R, Pset):
        nonlocal bestc
        if not Pset:
            if len(R) > len(bestc): bestc = R[:]
            return
        if len(R) + len(Pset) <= len(bestc): return
        for x in sorted(Pset):
            rec(R + [x], Pset & adj[x]); Pset = Pset - {x}
            if len(R) + len(Pset) <= len(bestc): return
    rec([], set(nodes)); return bestc
for i in range(n):
    Ai = marks[i]
    for C in Cs:
        lists = defaultdict(list)
        for j in range(i + 1, n):          # configurations whose least index is i
            T = Tof(Ai, marks[j]); D2 = d2(Ai, marks[j])
            if T not in divcache: divcache[T] = divisors(T)
            for v in divcache[T]:
                if okC(D2, v, C): lists[v].append(j)
        for v, L in lists.items():
            if len(L) < 2: continue
            adj = {x: set() for x in L}
            for a in range(len(L)):
                for b in range(a + 1, len(L)):
                    j, k = L[a], L[b]
                    if Tof(marks[j], marks[k]) % v == 0 and okC(d2(marks[j], marks[k]), v, C):
                        adj[j].add(k); adj[k].add(j)
            e = sum(len(adj[x]) for x in L) // 2
            if e == 0: continue
            count3[C] += sum(len(adj[x] & adj[y]) for x in L for y in adj[x] if y > x) // 3 * 0 + e
            cl = maxclique(L, adj)
            M = 1 + len(cl); sizes[C][M] += 1
            if M > best[C][0]: best[C] = (M, {"marks": [i] + cl, "v": v})
    if i % max(1, n // 16) == 0: print(i, round(time.time() - t0), {C: best[C][0] for C in Cs}, flush=True)
out = {"P": P, "cls": [r, m], "s": s, "best": {C: best[C] for C in Cs},
       "triples_through_least_index(edges in lists)": count3,
       "max clique size histogram per (least index, v)": {C: dict(sizes[C]) for C in Cs}}
json.dump(out, open(f"full_P{P}_{r}mod{m}.json", "w"), default=str)
print(json.dumps(out, default=str)); print("time", round(time.time() - t0))
