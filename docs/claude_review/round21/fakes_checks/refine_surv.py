"""Re-test saved survivors at a larger truncation T2 and report proximity to transplant
vectors for those that remain (evidence for Conjecture R)."""
import sys, json, math
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from mc_gapset import enum_vectors, survivors
from mc_survivors import directions, nearest
cfg = json.loads(sys.argv[1])
primes = cfg["primes"]; pts = np.array(json.load(open(cfg["file"])))
logs = [math.log(p) for p in primes]
dirs = [directions(p) for p in primes]
spacing = [(math.pi/2)/len(d) for d in dirs]
res = []
for T2 in cfg["Ts"]:
    V = enum_vectors(logs, T2); w = np.abs(V) @ np.array(logs)
    m = survivors(pts, V, w)
    keep = pts[m]
    nd = np.array([[nearest(t, dirs[j])[0]/spacing[j] for j, t in enumerate(ph)] for ph in keep]) if len(keep) else np.zeros((0, len(primes)))
    res.append(dict(T=T2, nV=int(V.shape[0]), remaining=int(m.sum()),
                    median_normdist=float(np.median(nd)) if len(nd) else None,
                    max_coord_normdist_median=float(np.median(nd.max(axis=1))) if len(nd) else None,
                    all_within_0_1=int((nd.max(axis=1) < 0.1).sum()) if len(nd) else 0))
print(json.dumps(dict(primes=primes, start=len(pts), refine=res)))
