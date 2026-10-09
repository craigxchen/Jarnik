#!/usr/bin/env python3
"""Referee check R11: distribution of the Liouville ratio |lambda|/(sqrt2 n^(-1/2)) over the 45 pairs of the
10-point cluster, and the a-priori cap C sqrt(n/(2R)) of Prop 3.2(2).  Tests the Section 7 claim that a
cluster must consist of pair forms 'within bounded factors of the elementary bound'."""
import math, itertools, importlib.util, sys
spec = importlib.util.spec_from_file_location("rc", "ref_cluster.py")
src = open("ref_cluster.py").read().split("rows = []")[0]
ns = {}
exec(src, ns)
P, N, C = ns["P"], ns["N"], ns["C"]
ggcd, gmul, gconj, gnorm = ns["ggcd"], ns["gmul"], ns["gconj"], ns["gnorm"]
R = math.sqrt(N)
rat = []
for z, w in itertools.combinations(P, 2):
    g = ggcd(z, w); n = N // gnorm(g)
    q = gmul(z, gconj(w)); lam = abs(math.atan2(q[1], q[0]))
    rat.append((lam/(math.sqrt(2)/math.sqrt(n)), n))
rat.sort()
for thr in [1.01, 1.5, 2, 5, 10, 100]:
    print(f"pairs with ratio <= {thr:>6}: {sum(1 for r, n in rat if r <= thr)} of 45")
print("pairs with n >= R (n/R >= 1):", sum(1 for r, n in rat if n >= R), "; ratio cap C sqrt(n/2R) max =",
      round(max(C*math.sqrt(n/(2*R)) for r, n in rat), 1))
