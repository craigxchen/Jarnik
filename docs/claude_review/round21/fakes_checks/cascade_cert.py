"""Certification of the cascade (fakes.md §3.2, Proposition B').

Claim checked for one step: background angles a_i = arg(beta_i) (beta_i genuine Gaussian primes,
fixed), new prime p.  Let A(t) be the set of t in T = R/(pi/2)Z such that every monomial
v = (v', n), n >= 1, satisfies the two gap inequalities.  We certify A(t) = {arg pi_p, -arg pi_p}.

Method (rigorous up to the correctness of mpmath at 50 digits, error << 1e-40):
 (1) all windows with w(v) <= T are generated; each is SHRUNK by 1e-40 on both sides, so the
     union of shrunk windows is contained in the union of the true windows;
 (2) the complement of the union of shrunk windows is computed; certification requires every
     uncovered interval to lie within distance rho of one of the two candidates t*;
 (3) at each candidate t*, exact Gaussian-integer arithmetic finds TIGHT monomials (|Re A| = 1,
     |Im A| = 1 or |Re A +- Im A| = 1): their windows have an endpoint exactly at t* and length
     2 eps/n; we need one on each side, each of length > rho + 1e-30;
 (4) then A(t) minus {t*} is covered near t*, and A(t) contains t* because a genuine
     Gaussian-integer tuple satisfies Lemma 5.1.
Hence A(t) = {candidates}.  Output: a JSON certificate summary.
"""
import sys, json, itertools
from mpmath import mp, mpf, asin, atan2, pi, sqrt, exp, log, floor
mp.dps = 50

def gauss_rep(q):
    for x in range(1, int(q ** 0.5) + 2):
        y2 = q - x * x
        if y2 <= 0: break
        y = int(round(y2 ** 0.5))
        for yy in (y - 1, y, y + 1):
            if yy > 0 and yy * yy == y2:
                return (x, yy)
    raise ValueError(q)

def gmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])

def gpow(a, k):
    r = (1, 0)
    for _ in range(k): r = gmul(r, a)
    return r

def conj(a): return (a[0], -a[1])

def enum(logs, T):
    r = len(logs); out = []; cur = [0] * r
    def rec(j, budget):
        if j == r:
            out.append(tuple(cur)); return
        kmax = int(budget / logs[j] + 1e-9)
        for k in range(-kmax, kmax + 1):
            cur[j] = k; rec(j + 1, budget - abs(k) * logs[j])
        cur[j] = 0
    rec(0, T)
    return out

def certify(bg_primes, bg_beta, p, T):
    q = pi / 2
    SH = mpf(10) ** -40
    bg_ang = [atan2(b[1], b[0]) for b in bg_beta]
    logs = [log(x) for x in bg_primes]
    lp = log(p)
    Vp = enum([float(x) for x in logs], float(T - lp)) if bg_primes else [()]
    wins = []
    for vp in Vp:
        wv = sum(abs(k) * L for k, L in zip(vp, logs)) if vp else mpf(0)
        ang = sum(k * a for k, a in zip(vp, bg_ang)) if vp else mpf(0)
        n = 1
        while wv + n * lp <= T:
            w = wv + n * lp
            ea = asin(exp(-w / 2)); ed = asin(exp(-w / 2) / sqrt(2))
            for off, eps in ((mpf(0), ea), (q / 2, ed)):
                rad = eps / n - SH
                base = (off - ang) % q
                for kk in range(n):
                    c = (base + kk * q) / n
                    wins.append((c - rad, c + rad))
            n += 1
    # wrap into [0, q)
    W2 = []
    for a, b in wins:
        if a < 0: W2.append((a + q, q)); W2.append((mpf(0), b))
        elif b > q: W2.append((a, q)); W2.append((mpf(0), b - q))
        else: W2.append((a, b))
    W2.sort(key=lambda x: x[0])
    gaps = []; run = mpf(0)
    for a, b in W2:
        if a > run: gaps.append((run, a))
        if b > run: run = b
    if run < q: gaps.append((run, q))
    x, y = gauss_rep(p)
    cands = [atan2(y, x) % q, (-atan2(y, x)) % q]
    pis = [(x, y), (x, -y)]
    # distance of each gap to candidates
    def cdist(t, c):
        d = abs(t - c) % q
        return min(d, q - d)
    rho = mpf(0); bad_gaps = []
    for a, b in gaps:
        mid = (a + b) / 2
        dmin = min(cdist(mid, c) for c in cands)
        ext = dmin + (b - a) / 2
        if ext > mpf(10) ** -25:
            bad_gaps.append((float(a), float(b)))
        rho = max(rho, ext)
    # tight monomials at each candidate (exact arithmetic)
    tight_info = []
    bgb = list(bg_beta)
    for ci, (c, pp) in enumerate(zip(cands, pis)):
        left = mpf(0); right = mpf(0); found = []
        for vp in Vp:
            wv = sum(abs(k) * float(L) for k, L in zip(vp, logs)) if vp else 0.0
            n = 1
            while wv + n * float(lp) <= float(T) + 1e-9:
                A = gpow(pp, n)
                for k, b in zip(vp, bgb):
                    if k > 0: A = gmul(A, gpow(b, k))
                    elif k < 0: A = gmul(A, gpow(conj(b), -k))
                X, Y = A
                w = (sum(abs(k) * L for k, L in zip(vp, logs)) if vp else mpf(0)) + n * lp
                for kind, tight in (("ax", min(abs(X), abs(Y)) == 1), ("dg", min(abs(X + Y), abs(X - Y)) == 1)):
                    if not tight: continue
                    eps = asin(exp(-w / 2)) if kind == "ax" else asin(exp(-w / 2) / sqrt(2))
                    val = n * c + sum(k * a for k, a in zip(vp, bg_ang)) if vp else n * c
                    lat = q * floor(val / q + mpf(1) / 2) if kind == "ax" else q * floor((val - q / 2) / q + mpf(1) / 2) + q / 2
                    s = val - lat
                    length = 2 * eps / n
                    if s > 0: left = max(left, length)
                    else: right = max(right, length)
                    found.append((vp, n, kind, "L" if s > 0 else "R", float(length)))
                n += 1
        ok = left > rho + mpf(10) ** -30 and right > rho + mpf(10) ** -30
        tight_info.append(dict(candidate=float(c), left=float(left), right=float(right), ok=bool(ok),
                               examples=found[:4]))
    cert = (not bad_gaps) and all(t["ok"] for t in tight_info)
    return dict(bg=bg_primes, bg_beta=bg_beta, p=p, T=float(T), windows=len(wins),
                uncovered_intervals=len(gaps), rho=float(rho), bad_gaps=bad_gaps[:5],
                tight=tight_info, certified=bool(cert))

if __name__ == "__main__":
    cfg = json.loads(sys.argv[1])
    res = certify(cfg["bg"], [tuple(b) for b in cfg["bg_beta"]], cfg["p"], mpf(cfg["T"]))
    print(json.dumps(res))
