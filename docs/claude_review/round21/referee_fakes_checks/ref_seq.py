"""Independent re-verification of sequentially built vectors (referee).
Enumerates all v != 0 in Z^r (modulo sign) with w(v) <= T, written from scratch, and computes the
minimum ratio dist/threshold for both gap conditions. Doubles suffice: thresholds >= e^(-T/2) >~ 1e-8,
rounding error of <v,phi> <= sum|v_j| * 1e-16 * 2."""
import sys, json, math
import numpy as np

def enum(logs, T):
    r = len(logs); out = []
    cur = [0]*r
    # canonical sign: first nonzero coordinate positive
    def rec(i, rem, started):
        if i == r:
            if started: out.append(tuple(cur))
            return
        nmax = int(rem // logs[i] + 1e-12)
        rng = range(-nmax, nmax+1) if started else range(0, nmax+1)
        for n in rng:
            cur[i] = n
            rec(i+1, rem - abs(n)*logs[i], started or n != 0)
        cur[i] = 0
    rec(0, T, False)
    return np.array(out, dtype=np.int64)

def check(primes, phis, T, chunk=200000):
    logs = [math.log(p) for p in primes]
    V = enum(logs, T)
    w = np.abs(V) @ np.array(logs)
    res = []
    hp = math.pi/2
    for phi in phis:
        phi = np.array(phi)
        worst = np.inf; worst_v = None; nbad = 0
        for s in range(0, len(V), chunk):
            Vc = V[s:s+chunk]; wc = w[s:s+chunk]
            t = Vc @ phi
            ra = np.mod(t, hp); da = np.minimum(ra, hp - ra)
            rd = np.mod(t - math.pi/4, hp); dd = np.minimum(rd, hp - rd)
            ta = np.arcsin(np.exp(-wc/2)); td = np.arcsin(np.exp(-wc/2)/math.sqrt(2))
            ratio = np.minimum(da/ta, dd/td)
            k = int(np.argmin(ratio))
            nbad += int((ratio < 1).sum())
            if ratio[k] < worst: worst = float(ratio[k]); worst_v = Vc[k].tolist()
        res.append(dict(min_ratio=worst, worst_v=worst_v, nfail=nbad))
    return len(V), res

if __name__ == "__main__":
    cfg = json.loads(sys.argv[1])
    rows = [json.loads(l) for l in open(cfg["file"])]
    row = [r for r in rows if r["primes"] == cfg["primes"]][0]
    phis = [r["angles"] for r in row["results"] if r.get("complete")]
    flags = [r.get("verified_T2") for r in row["results"] if r.get("complete")]
    for T in cfg["Ts"]:
        n, res = check(row["primes"], phis, T)
        print(json.dumps(dict(primes=row["primes"], T=T, nV=n, flags=flags, res=res)), flush=True)
