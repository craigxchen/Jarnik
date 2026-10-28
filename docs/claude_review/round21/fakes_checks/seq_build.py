"""Sequential construction of gap-admissible angle vectors far beyond the union-bound
regime (evidence for 'Poisson positivity', §5.5 of fakes.md).

phi_1, ..., phi_r are chosen one at a time: phi_k is drawn uniformly from the admissible
set A_k(phi_{<k}) = { t : all v with last non-zero coordinate k and w(v) <= T pass },
computed exactly as a union of intervals (next_angle.allowed_set).  The product of the
measures mu_1(A_k) along the path is an unbiased estimate of the Haar measure of the
truncated set G_T (sequential importance sampling).  The final vector is re-verified
against ALL monomials with w(v) <= T2 (T2 > T) and its distance to transplant vectors is
recorded.  Evidence only (double precision).
"""
import sys, json, math
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from next_angle import allowed_set
from mc_gapset import enum_vectors, survivors
from mc_survivors import directions, nearest

def pick(comps, rng, shrink=0.1):
    lens = np.array([b - a for a, b in comps])
    i = rng.choice(len(comps), p=lens / lens.sum())
    a, b = comps[i]
    return a + (b - a) * (shrink + (1 - 2 * shrink) * rng.random())

def build(primes, T, rng):
    q = math.pi / 2
    angs = []; logmeas = 0.0; meas = []
    for k, p in enumerate(primes):
        comps, _ = allowed_set(p, primes[:k], angs, T)
        comps = [(a, b) for a, b in comps if b - a > 1e-12]
        if not comps:
            return None, logmeas, meas
        m = sum(b - a for a, b in comps) / q
        meas.append(m); logmeas += math.log(m)
        angs.append(pick(comps, rng))
    return angs, logmeas, meas

def main():
    cfg = json.loads(sys.argv[1])
    primes = cfg["primes"]; T = cfg["T"]; T2 = cfg["T2"]; tries = cfg.get("tries", 5)
    rng = np.random.default_rng(cfg.get("seed", 11))
    logs = [math.log(p) for p in primes]
    V2 = enum_vectors(logs, T2); w2 = np.abs(V2) @ np.array(logs)
    Lam = sum(p ** -0.5 for p in primes)
    Pi = 1.0
    for p in primes: Pi *= 1 + 2 / (math.sqrt(p) - 1)
    dirs = [directions(p) for p in primes]
    spacing = [(math.pi / 2) / len(d) for d in dirs]
    results = []
    for t in range(tries):
        angs, lm, meas = build(primes, T, rng)
        rec = dict(try_=t, reached=len(meas), log_path_measure=lm,
                   step_measures=[round(x, 6) for x in meas])
        if angs is not None:
            ok = bool(survivors(np.array([angs]), V2, w2, tol=1e-11)[0])
            nd = [nearest(a, dirs[j])[0] / spacing[j] for j, a in enumerate(angs)]
            rec.update(complete=True, verified_T2=ok, angles=angs,
                       median_normdist_to_transplant=float(np.median(nd)),
                       max_normdist=float(max(nd)))
        else:
            rec.update(complete=False)
        results.append(rec)
    print(json.dumps(dict(primes=primes, Lambda=Lam, Pi=Pi, T=T, T2=T2, nV2=int(V2.shape[0]),
                          results=results)))

if __name__ == "__main__":
    main()
