"""CG d=6 (sum sigma=0, 20 points): K.Gamma~ on X_9 for a random sample of 9-subtuples
(the full enumeration of 167960 subtuples in ref_curve_intersections.py was cut by its timeout)."""
import sys, random
sys.path.insert(0, __file__.rsplit('/', 1)[0])
sys.argv = [sys.argv[0]]
import ref_curve_intersections as R
sig, pts = R.cg_family(6, 0)
rng = random.Random(7); vals = set(); n = 0
for _ in range(3000):
    S = sorted(rng.sample(range(20), 9))
    cg, (e, fin, inf) = R.curve_numbers([pts[i] for i in S], 6)
    assert cg[0] == 0
    vals.add((e, fin, inf, -12 + 7 * e)); n += 1
print(f"CG d=6 on X_9: {n} random 9-subtuples: (E.G, finite, at inf, K.G) values = {sorted(vals)}")
