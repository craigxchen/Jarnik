"""Test of the rigidity picture (Conjecture R): are the Haar-random survivors of the
truncated gap test close to TRANSPLANT vectors (arg beta_j, beta_j primitive Gaussian
integer of odd norm <= p_j)?  For each survivor and coordinate j we compute the distance
to the nearest such direction (mod pi/2), normalised by the mean spacing of these
directions, and compare with uniform random points.  Also records the joint event
"every coordinate within 0.1 spacing of a transplant direction" and whether the norms of the
nearest directions are pairwise coprime (Lemma T).  Evidence only.
"""
import sys, json, math
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from mc_gapset import enum_vectors, survivors

def directions(p):
    out = []
    s = int(math.isqrt(p))
    for x in range(0, s + 1):
        for y in range(1, s + 1):
            n = x * x + y * y
            if n <= p and n % 2 == 1 and n > 1 and math.gcd(x, y) == 1:
                out.append((math.atan2(y, x) % (math.pi / 2), n))
    out.sort()
    return out

def nearest(t, dirs):
    q = math.pi / 2
    best = (9, None)
    for th, n in dirs:
        d = abs(t - th); d = min(d, q - d)
        if d < best[0]:
            best = (d, n)
    return best

def main():
    cfg = json.loads(sys.argv[1])
    primes = cfg["primes"]; T = cfg["T"]; n = cfg["n"]; seed = cfg.get("seed", 7)
    logs = [math.log(p) for p in primes]
    V = enum_vectors(logs, T); wts = np.abs(V) @ np.array(logs)
    rng = np.random.default_rng(seed)
    surv = []
    B = 5000; done = 0
    while done < n:
        b = min(B, n - done)
        ph = rng.uniform(0, math.pi / 2, size=(b, len(primes)))
        m = survivors(ph, V, wts)
        surv += ph[m].tolist(); done += b
    dirs = [directions(p) for p in primes]
    spacing = [(math.pi / 2) / len(d) for d in dirs]
    def stats(points):
        normd = []; close_all = 0; coprime_close = 0
        for ph in points:
            ds = []; ns = []
            for j, t in enumerate(ph):
                d, nn = nearest(t, dirs[j]); ds.append(d / spacing[j]); ns.append(nn)
            normd.append(ds)
            if max(ds) < 0.1:
                close_all += 1
                if all(math.gcd(ns[a], ns[b]) == 1 for a in range(len(ns)) for b in range(a)):
                    coprime_close += 1
        a = np.array(normd) if normd else np.zeros((0, len(primes)))
        return dict(count=len(points), median_normdist=float(np.median(a)) if len(a) else None,
                    frac_coord_within_0_05=float((a < 0.05).mean()) if len(a) else None,
                    all_coords_within_0_1=close_all, of_which_coprime_norms=coprime_close)
    rand = rng.uniform(0, math.pi / 2, size=(max(len(surv), 2000), len(primes))).tolist()
    out = dict(primes=primes, T=T, n=n, survivors=stats(surv), uniform=stats(rand),
               ndirs=[len(d) for d in dirs])
    print(json.dumps(out))
    with open(cfg.get("dump", "/dev/null"), "w") as f:
        json.dump(surv, f)

if __name__ == "__main__":
    main()
