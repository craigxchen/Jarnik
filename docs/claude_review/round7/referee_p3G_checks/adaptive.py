"""Adaptive (greedy) prime selection against the pinning closure (evidence only).

Accept a candidate prime iff, together with the columns already accepted, NO certificate of Lemma 2.4
relative to Pi_0 = {0} exists among weight-<=2 relations of the following kinds:
  * same-v pairs (k=2, n=(1,1)) for singles, doubles and pairs (the author's 'pair' certificates);
  * double/single combinations 2e_a vs e_a with n=(1,2);
  * ANY closed walk (any length, multiplicities allowed) of pair and single relations in the
    signed double cover, with total arc weight sum n_i (pi - 2 d_i) < pi.
If no such certificate exists, Pi_1 = {0}, hence the corresponding sub-closure is {0}: all columns free.
usage: python3 adaptive.py M r K seed maxcand
"""
import sys, math, random, heapq, json, time
from rcore import is_prime, gauss, Lprime, gmul, conj, Ds, sigma_exact, sin2d

def tight_units(A, L):
    out = []
    for s, D in enumerate(Ds(A)):
        if D != 0 and L % abs(D) == 0:
            s2 = sin2d(A, s)
            out.append((s, sigma_exact(A, s), s2, math.pi - 2 * math.asin(math.sqrt(float(s2)))))
    return out

def seed_in(units):
    for i in range(len(units)):
        for j in range(i + 1, len(units)):
            if units[i][1] != units[j][1] and units[i][2] + units[j][2] > 1: return True
    return False

class Greedy:
    def __init__(self, L):
        self.L = L; self.gp = []; self.adj = {}  # node -> list of (node, w)
        self.G = ('G',)
    def edges_for(self, k, z, cand_pairs, singles):
        E = []
        for (j, t, units) in cand_pairs:          # v = e_j + t e_k
            for (s, sg, s2, w) in units:
                E.append(((j, -sg), (k, sg * t), w)); E.append(((k, -sg * t), (j, sg), w))
        for (s, sg, s2, w) in singles:            # v = e_k
            E.append((self.G, (k, sg), w)); E.append(((k, -sg), self.G, w))
        return E
    def shortest_cycle_through(self, x, adj):
        dist = {x: 0.0}; pq = [(0.0, x)]; best = math.inf
        while pq:
            d, u = heapq.heappop(pq)
            if d > dist.get(u, math.inf) or d >= math.pi: continue
            for (y, w) in adj.get(u, ()):
                nd = d + w
                if y == x: best = min(best, nd); continue
                if nd < dist.get(y, math.inf) and nd < math.pi:
                    dist[y] = nd; heapq.heappush(pq, (nd, y))
        return best
    def try_add(self, p):
        z = gauss(p); k = len(self.gp); L = self.L
        singles = tight_units(z, L); doubles = tight_units(gmul(z, z), L)
        if seed_in(singles) or seed_in(doubles): return 'seed1'
        for (s1, g1, a1, w1) in singles:          # n=(1,2): sigma_d 2e + 2 sigma_s e = 0
            for (s2, g2, a2, w2) in doubles:
                if g1 != g2 and 2 * w1 + w2 < math.pi - 1e-9: return 'dbl'
        cand_pairs = []
        for j, zj in enumerate(self.gp):
            for t, A in ((1, gmul(zj, z)), (-1, gmul(zj, conj(z)))):
                u = tight_units(A, L)
                if u:
                    if seed_in(u): return 'seed2'
                    cand_pairs.append((j, t, u))
        E = self.edges_for(k, z, cand_pairs, singles)
        adj = {a: list(b) for a, b in self.adj.items()}
        for (a, b, w) in E: adj.setdefault(a, []).append((b, w))
        for x in ((k, 1), (k, -1)):
            c = self.shortest_cycle_through(x, adj)
            if c < math.pi - 1e-9: return 'walk'
            if c < math.pi + 1e-9: return 'tie'
        self.adj = adj; self.gp.append(z); return 'ok'

if __name__ == '__main__':
    M, r, K, seed, maxc = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    P = int(M ** K); X = 2 * M; L = Lprime(X, set())
    rng = random.Random(seed); g = Greedy(L); seen = set(); stats = {}; ncand = 0; t = time.time(); ps = []
    while len(g.gp) < r and ncand < maxc:
        n = rng.randrange(P, 2 * P) | 1
        if n % 4 != 1 or n in seen or not is_prime(n): continue
        seen.add(n); ncand += 1
        res = g.try_add(n); stats[res] = stats.get(res, 0) + 1
        if res == 'ok': ps.append(n)
    print(json.dumps(dict(M=M, K=K, seed=seed, target_r=r, got=len(ps), candidates=ncand, outcomes=stats,
                          secs=round(time.time() - t, 1))), flush=True)
    with open(f'adaptive_primes_M{M}_K{K}_s{seed}.json', 'w') as f: json.dump(ps, f)
