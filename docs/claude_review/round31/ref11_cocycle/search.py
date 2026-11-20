"""Search for integral-cofactor configurations (one unit class) on q_P.

Criterion (proved in the report): a set S of distinct marks (one unit class,
rows pairwise distinct) is an integral-cofactor configuration for some
lambda = u/v, w iff  diam(S) <= C q^{1/4} * G(S),  G(S) = gcd_{pairs} T_ij,
and then the best choice is lambda = 1/G (u = 1, v = G).
|T_ij| = gcd(dx, dy)/2 exactly (content of 2 i eta H is 2).

Candidate generation: smooth part of G over odd moduli m (primes and
small prime powers coprime to q): A_i = A_j mod m  <=>  m | T_ij.
Every candidate is then checked exactly with Python integers.
"""
import math, pickle, sys, time
import numpy as np
from sympy import primerange

P = int(sys.argv[1]); CMAX = float(sys.argv[2]); i0 = int(sys.argv[3]); i1 = int(sys.argv[4])
ps, pis, q, marks = pickle.load(open(f"marks_{P}.pkl", "rb"))
N = len(marks)
q4 = q ** 0.25
mods, logs = [], []
for l in primerange(3, 600):
    if q % l == 0:
        continue
    e = 1
    while l ** e < 3000 and len(mods) < 200:
        mods.append(l ** e); logs.append(math.log(l)); e += 1
modl = [int(m) for m in mods]; mods = np.array(mods); logs = np.array(logs)
K = len(mods)
R = np.empty((N, K), dtype=np.int64)
for n, (x, y) in enumerate(marks):
    R[n] = [(x % m) * m + (y % m) for m in modl]
X = np.array([float(x) for x, y in marks]); Y = np.array([float(y) for x, y in marks])

def T(a, b):
    dx = marks[b][0] - marks[a][0]; dy = marks[b][1] - marks[a][1]
    g = math.gcd(dx, dy)
    assert g % 2 == 0
    return g // 2

def exact(S):
    G = 0; D2 = 0
    for a in range(len(S)):
        for b in range(a + 1, len(S)):
            G = math.gcd(G, T(S[a], S[b]))
            dx = marks[S[b]][0] - marks[S[a]][0]; dy = marks[S[b]][1] - marks[S[a]][1]
            D2 = max(D2, dx * dx + dy * dy)
    # C needed: diam / (G q^{1/4}); exact compare done by caller via squares
    return G, D2

out = open(f"found_{P}_{i0}_{i1}.txt", "w")
t0 = time.time()
nfound = 0
for i in range(i0, min(i1, N)):
    E = (R[i + 1:] == R[i])
    w = E @ logs
    d = np.hypot(X[i + 1:] - X[i], Y[i + 1:] - Y[i])
    need = np.log(np.maximum(d / (CMAX * q4), 1.0))
    J = np.nonzero(w >= need - 1e-9)[0]
    if len(J) < 2:
        continue
    EJ = E[J]; jj = J + i + 1
    for a in range(len(J)):
        common = EJ[a] & EJ[a + 1:]
        wc = common @ logs
        for bidx in np.nonzero(wc >= 0)[0]:
            b = a + 1 + bidx
            j, k = int(jj[a]), int(jj[b])
            dd = max(d[J[a]], d[J[b]], math.hypot(X[k] - X[j], Y[k] - Y[j]))
            if wc[bidx] < math.log(max(dd / (CMAX * q4), 1.0)) - 1e-9:
                continue
            S = [i, j, k]
            G, D2 = exact(S)
            if D2 * D2 <= int(CMAX) ** 4 * q * G ** 4:   # exact: diam <= CMAX q^(1/4) G
                cneed = math.sqrt(D2) / (G * q4)
                # try to extend to 4 points using remaining candidates
                ext = []
                for c in range(len(J)):
                    l4 = int(jj[c])
                    if l4 in (j, k):
                        continue
                    G4, D4 = exact(S + [l4])
                    c4 = math.sqrt(D4) / (G4 * q4)
                    if c4 <= CMAX:
                        ext.append((l4, G4, c4))
                out.write(f"{i} {j} {k} G={G} diam/q^(1/4)={math.sqrt(D2)/q4:.6g} Cneed={cneed:.6g} ext={ext}\n")
                out.flush(); nfound += 1
    if (i - i0) % 2000 == 0:
        print(i, nfound, round(time.time() - t0), flush=True)
print("done", nfound, round(time.time() - t0), flush=True)
