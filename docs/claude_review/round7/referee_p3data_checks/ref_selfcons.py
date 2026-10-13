"""Referee: §5.1 self-consistency of the M=8 adapted actual-prime profiles, with and without winding freedom.
For each profile (q_results_M8: prime_cols and prime_sep_cols), each random conjugation choice:
compute z_x = prod pi_j^{a_xj} conj(pi_j)^{1-a_xj}; for each pair and each unit eta in {1,i,-1,-i} the collision
modulus m^eta_xy = prod_q q^{level} (z_x = eta z_y mod q^level, q^level <= X). With windings k:
unit-1 demand of pair = m^{i^{k_x-k_y}}. Count conj choices feasible (all pairs Pi_sep >= m^4 Pi_same)
(a) with k = 0, (b) for some k in (Z/4)^M (k_0 = 0)."""
import json, math, random, itertools, sys
import numpy as np
from sympy import primerange
P8 = '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round7/p3data_checks/q_results_M8.jsonl'
M = 8
def dec(mask): return [-1 if (mask >> x) & 1 else 1 for x in range(M)]
def two_sq(p):
    for a in range(1, int(math.isqrt(p)) + 1, 1):
        b2 = p - a*a; b = math.isqrt(b2)
        if b*b == b2: return a, b
def gmul(a, b, Q): return ((a[0]*b[0] - a[1]*b[1]) % Q, (a[0]*b[1] + a[1]*b[0]) % Q)
UNITS = [(1, 0), (0, 1), (-1, 0), (0, -1)]  # i^0..i^3
K = np.array(list(itertools.product(range(4), repeat=M-1)))
K = np.hstack([np.zeros((K.shape[0], 1), dtype=int), K])
ntr = int(sys.argv[1]) if len(sys.argv) > 1 else 200
rng = random.Random(12345)
tot = {}
for line in open(P8):
    R = json.loads(line); X = R['X']
    for key in ('prime_cols', 'prime_sep_cols'):
        cols = R[key]
        S = np.array([dec(m) for m, _ in cols]).T; a = (1 + S)//2
        prims = [pl for _, pl in cols]
        # physical primes per column (a column may carry a product: then each prime is its own physical column with same cut)
        phys = [(j, p) for j, pl in enumerate(prims) for p in pl]
        Pcol = [math.prod(pl) for pl in prims]; N = math.prod(Pcol)
        sepP = {}; 
        for x in range(M):
            for y in range(x+1, M):
                sep = math.prod(Pcol[j] for j in range(len(Pcol)) if S[x, j] != S[y, j]); sepP[(x, y)] = (sep, N//sep)
        gs = [two_sq(p) for _, p in phys]
        qs = []
        for q in primerange(3, X+1):
            A = 1
            while q**(A+1) <= X: A += 1
            qs.append((q, A))
        f0 = 0; fk = 0; viol0 = []
        for tr in range(ntr):
            g = [(u, v if rng.random() < 0.5 else -v) for (u, v) in gs]
            # m[eta][x][y]
            m = np.ones((4, M, M), dtype=object)
            for q, A in qs:
                Q = q**A
                z = []
                for x in range(M):
                    val = (1, 0)
                    for (j, p), (u, v) in zip(phys, g):
                        val = gmul(val, (u, v) if a[x, j] == 1 else (u, -v), Q)
                    z.append(val)
                for e, eta in enumerate(UNITS):
                    for x in range(M):
                        for y in range(M):
                            if x == y: continue
                            w = gmul(eta, z[y], Q); lev = 0
                            for L in range(1, A+1):
                                QL = q**L
                                if (z[x][0]-w[0]) % QL == 0 and (z[x][1]-w[1]) % QL == 0: lev = L
                                else: break
                            m[e, x, y] *= q**lev
            ok = np.zeros((4, M, M), dtype=bool)
            for (x, y), (sep, same) in sepP.items():
                for e in range(4):
                    ok[e, x, y] = sep >= int(m[e, x, y])**4*same   # z_x = i^e z_y
            nv0 = sum(not ok[0, x, y] for (x, y) in sepP)
            viol0.append(nv0)
            if nv0 == 0: f0 += 1
            # with windings: pair (x,y) unit-1 collision of zhat means z_x = i^{k_x-k_y} z_y
            allok = np.ones(K.shape[0], dtype=bool)
            for (x, y) in sepP:
                e = (K[:, x] - K[:, y]) % 4
                allok &= ok[e, x, y]
            if allok.any(): fk += 1
        print(f"X={X} t={R['t']} {key}: feasible k=0: {f0}/{ntr} (mean violated pairs {np.mean(viol0):.2f}); feasible for some k: {fk}/{ntr}", flush=True)
