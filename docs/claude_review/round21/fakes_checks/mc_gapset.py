"""Monte Carlo estimate of the Haar measure of the gap-admissible set G(p) in T^r.

G(p) = { phi : for all v != 0,
           dist(<v,phi>, (pi/2) Z)        >= arcsin(e^{-w(v)/2}),
           dist(<v,phi>, pi/4 + (pi/2) Z) >= arcsin(e^{-w(v)/2}/sqrt 2) },
w(v) = sum |v_j| log p_j  (overview Lemma 5.1).

We truncate to w(v) <= T; survivors of the truncated test form a superset of G(p),
so the estimate is an UPPER bound for mu(G(p)) up to sampling error, and the
tail beyond T can remove at most sum_{w(v)>T} (4/pi)(1+1/sqrt2) e^{-w/2} more.
Evidence only.
"""
import sys, math, json
import numpy as np

def split_primes(n_needed, start=5):
    out = []
    q = start
    while len(out) < n_needed:
        if q % 4 == 1 and all(q % d for d in range(2, int(q**0.5) + 1)):
            out.append(q)
        q += 1
    return out

def enum_vectors(logs, T):
    """all v in Z^r, v != 0, first nonzero coordinate > 0, w(v) <= T."""
    r = len(logs)
    res = []
    cur = [0] * r
    def rec(j, budget, started):
        if j == r:
            if started:
                res.append(list(cur))
            return
        kmax = int(budget / logs[j] + 1e-12)
        rng = range(0, kmax + 1) if not started else range(-kmax, kmax + 1)
        for k in rng:
            cur[j] = k
            rec(j + 1, budget - abs(k) * logs[j], started or k != 0)
        cur[j] = 0
    rec(0, T, False)
    return np.array(res, dtype=np.int64)

def tail_bound(logs, T, s=0.25):
    # sum_{w(v) > T} e^{-w/2} <= e^{-(1/2-s)T} prod (1+2/(p^s-1)) (Rankin)
    Z = 1.0
    for L in logs:
        ps = math.exp(s * L)
        Z *= 1 + 2 / (ps - 1)
    return (4 / math.pi) * (1 + 1 / math.sqrt(2)) * math.exp(-(0.5 - s) * T) * Z / 2

def survivors(phis, V, wts, tol=1e-12):
    """boolean mask of phis (B x r) surviving all v in V."""
    eps_ax = np.arcsin(np.exp(-wts / 2))
    eps_dg = np.arcsin(np.exp(-wts / 2) / math.sqrt(2))
    alive = np.ones(phis.shape[0], dtype=bool)
    q = math.pi / 2
    CH = 4000
    for i in range(0, V.shape[0], CH):
        Vc = V[i:i + CH].astype(np.float64)
        X = phis @ Vc.T  # B x ch
        Xm = np.mod(X, q)
        dax = np.minimum(Xm, q - Xm)
        ddg = np.abs(Xm - q / 2)
        bad = (dax < eps_ax[i:i + CH] - tol) | (ddg < eps_dg[i:i + CH] - tol)
        alive &= ~bad.any(axis=1)
        if not alive.any():
            break
    return alive

def main():
    cfg = json.loads(sys.argv[1])
    primes = cfg["primes"]
    T = cfg["T"]
    nsamp = cfg["n"]
    seed = cfg.get("seed", 1)
    logs = [math.log(p) for p in primes]
    V = enum_vectors(logs, T)
    wts = np.abs(V) @ np.array(logs)
    Lam = sum(p ** -0.5 for p in primes)
    Pi = 1.0
    for p in primes:
        Pi *= 1 + 2 / (math.sqrt(p) - 1)
    rng = np.random.default_rng(seed)
    alive_total = 0
    B = cfg.get("batch", 2000)
    done = 0
    while done < nsamp:
        b = min(B, nsamp - done)
        phis = rng.uniform(0, math.pi / 2, size=(b, len(primes)))
        alive_total += int(survivors(phis, V, wts).sum())
        done += b
    # genuine sanity check
    gen = []
    for p in primes:
        a = int(math.isqrt(p))
        found = None
        for x in range(1, a + 1):
            y2 = p - x * x
            y = math.isqrt(y2)
            if y * y == y2 and y > 0:
                found = (x, y); break
        gen.append(math.atan2(found[1], found[0]))
    gen_ok = bool(survivors(np.array([gen]), V, wts, tol=1e-9)[0])
    out = dict(primes=primes, T=T, nV=int(V.shape[0]), Lambda=Lam, Pi=Pi,
               n=nsamp, alive=alive_total, frac=alive_total / nsamp,
               tail=tail_bound(logs, T), genuine_survives=gen_ok)
    print(json.dumps(out))

if __name__ == "__main__":
    main()
