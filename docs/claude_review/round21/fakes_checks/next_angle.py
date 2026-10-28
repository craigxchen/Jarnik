"""Allowed set for one new angle t = phi_p, given a FIXED background angle vector
phi' for primes q_1..q_k (genuine: arg of a Gaussian prime above q_i, or Haar-random).

Constraint (overview Lemma 5.1) for v = (v', n), n >= 1 (n <= -1 is the same by v -> -v):
  dist(n t + <v',phi'>, (pi/2)Z)        >= arcsin(e^{-(w(v') + n log p)/2})
  dist(n t + <v',phi'>, pi/4 + (pi/2)Z) >= arcsin(e^{-(w(v') + n log p)/2}/sqrt2)
Truncated at w(v') + n log p <= T, so the computed allowed set is a SUPERSET of the
true allowed set.  Reports its Haar measure (fraction of [0, pi/2)), the number of
components, and which Gaussian directions arg(beta), Norm beta <= p, beta primitive,
lie in it.  Evidence only (floating point, double precision; windows below 1e-13
are dropped and reported).
"""
import sys, math, json
import numpy as np

def gaussian_rep(q):
    for x in range(1, int(math.isqrt(q)) + 1):
        y2 = q - x * x
        y = math.isqrt(y2)
        if y * y == y2 and y > 0:
            return x, y
    raise ValueError

def enum_vectors_all(logs, T):
    """all v' in Z^k (including 0) with w(v') <= T."""
    r = len(logs); res = []; cur = [0] * r
    def rec(j, budget):
        if j == r:
            res.append(list(cur)); return
        kmax = int(budget / logs[j] + 1e-12)
        for k in range(-kmax, kmax + 1):
            cur[j] = k; rec(j + 1, budget - abs(k) * logs[j])
        cur[j] = 0
    rec(0, T)
    return np.array(res, dtype=np.int64).reshape(-1, r)

def allowed_set(p, bg_primes, bg_angles, T):
    q = math.pi / 2
    logs = [math.log(x) for x in bg_primes]
    lp = math.log(p)
    if bg_primes:
        V = enum_vectors_all(logs, T - lp)
        wv = np.abs(V) @ np.array(logs)
        ang = V @ np.array(bg_angles)
    else:
        wv = np.zeros(1); ang = np.zeros(1)
    starts = []; ends = []
    dropped = 0
    n = 1
    while n * lp <= T + 1e-12:
        sel = wv + n * lp <= T + 1e-12
        w = wv[sel] + n * lp
        a = ang[sel]
        ea = np.arcsin(np.exp(-w / 2)); ed = np.arcsin(np.exp(-w / 2) / math.sqrt(2))
        for centre_off, eps in ((0.0, ea), (q / 2, ed)):
            keep = eps > 1e-13
            dropped += int((~keep).sum())
            # n t = -a + centre_off + k q  ->  t = (-a + centre_off + k q)/n, k = 0..n-1 mod q
            base = np.mod(-a[keep] + centre_off, q)
            e = eps[keep]
            for kk in range(n):
                c = (base + kk * q) / n
                starts.append(c - e / n); ends.append(c + e / n)
        n += 1
    s = np.concatenate(starts); e = np.concatenate(ends)
    # wrap into [0, q)
    s2 = []; e2 = []
    for (lo, hi) in ((s, e),):
        m1 = lo < 0
        s2 += [np.maximum(lo, 0), lo[m1] + q]; e2 += [np.minimum(hi, q), np.full(m1.sum(), q)]
        m2 = hi > q
        s2 += [np.zeros(m2.sum())]; e2 += [hi[m2] - q]
    s = np.concatenate(s2); e = np.concatenate(e2)
    order = np.argsort(s); s = s[order]; e = e[order]
    comps = []
    cur = 0.0
    run_end = 0.0
    for lo, hi in zip(s, e):
        if lo > run_end:
            comps.append((run_end, lo))
        run_end = max(run_end, hi)
    if run_end < q:
        comps.append((run_end, q))
    return comps, dropped

def main():
    cfg = json.loads(sys.argv[1])
    p = cfg["p"]; bg = cfg["bg"]; T = cfg["T"]; mode = cfg.get("mode", "genuine")
    rng = np.random.default_rng(cfg.get("seed", 0))
    if mode == "genuine":
        angs = [math.atan2(gaussian_rep(x)[1], gaussian_rep(x)[0]) for x in bg]
    else:
        angs = list(rng.uniform(0, math.pi / 2, size=len(bg)))
    comps, dropped = allowed_set(p, bg, angs, T)
    q = math.pi / 2
    meas = sum(b - a for a, b in comps) / q
    # Gaussian directions of norm <= p, primitive, modulo units: angles in [0, pi/2)
    dirs = []
    for x in range(0, int(math.isqrt(p)) + 1):
        for y in range(1, int(math.isqrt(p)) + 1):
            if x * x + y * y <= p and math.gcd(x, y) == 1:
                dirs.append((x * x + y * y, math.atan2(y, x) % q, (x, y)))
    hit = []
    for nrm, th, xy in dirs:
        for a, b in comps:
            if a - 1e-12 <= th <= b + 1e-12:
                hit.append((nrm, xy)); break
    big = sorted(comps, key=lambda c: c[0] - c[1])[:5]
    lam = sum(x ** -0.5 for x in bg)
    print(json.dumps(dict(p=p, bg=bg, mode=mode, T=T, Lambda_bg=lam, measure=meas,
                          components=len(comps), dropped_tiny=dropped,
                          largest=[(round(a, 8), round(b - a, 10)) for a, b in big],
                          gaussian_dirs_in_set=hit)))

if __name__ == "__main__":
    main()
