"""Sequential construction with one-step lookahead (evidence).  At step k, K candidate values
of phi_k are drawn uniformly from A_k(phi_{<k}); the one maximising the Haar measure of the
next admissible set A_{k+1} is kept.  Final vector re-verified at T2."""
import sys, json, math
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from next_angle import allowed_set
from mc_gapset import enum_vectors, survivors
from seq_build import pick

def main():
    cfg = json.loads(sys.argv[1])
    ps = cfg["primes"]; T = cfg["T"]; T2 = cfg["T2"]; K = cfg.get("K", 6)
    rng = np.random.default_rng(cfg.get("seed", 5))
    q = math.pi / 2
    angs = []; log = []
    comps, _ = allowed_set(ps[0], [], [], T)
    for k in range(len(ps)):
        comps = [(a, b) for a, b in comps if b - a > 1e-12]
        if not comps:
            print(json.dumps(dict(primes=ps, T=T, died_at=k, log=log))); return
        if k == len(ps) - 1:
            angs.append(pick(comps, rng)); log.append(dict(k=k, meas=sum(b - a for a, b in comps) / q)); break
        best = None
        for _ in range(K):
            t = pick(comps, rng)
            nxt, _ = allowed_set(ps[k + 1], ps[:k + 1], angs + [t], T)
            nxt = [(a, b) for a, b in nxt if b - a > 1e-12]
            m = sum(b - a for a, b in nxt) / q
            if best is None or m > best[0]:
                best = (m, t, nxt)
        log.append(dict(k=k, meas=sum(b - a for a, b in comps) / q, next_best=best[0]))
        angs.append(best[1]); comps = best[2]
    logs = [math.log(p) for p in ps]
    V2 = enum_vectors(logs, T2); w2 = np.abs(V2) @ np.array(logs)
    ok = bool(survivors(np.array([angs]), V2, w2, tol=1e-11)[0])
    print(json.dumps(dict(primes=ps, T=T, T2=T2, nV2=int(V2.shape[0]), complete=True, verified_T2=ok,
                          angles=angs, log=log)))

if __name__ == "__main__":
    main()
