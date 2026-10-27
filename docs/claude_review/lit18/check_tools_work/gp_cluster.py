"""Gaussian-prime cluster test against the sparse-regime hypotheses of tools.md (QS^sparse, near-1 count, 2-function gcd).
For target s and log-height L0, find s Gaussian primes pi_j (first octant, pairwise non-associate, distinct prime norms)
inside one small square at modulus ~ exp(L0/2). Then S = places above the s split primes N(pi_j), |P| = s, and:
 (1) every ordered m-tuple of distinct pi_j solves the round-8 system (5.1) with delta up to delta_eff_min (computed as a
     rigorous worst case over ALL tuples), and H(x) >= exp(L0) >= exp(A s log s);
 (2) Chebyshev / determination: any family of proper subspaces containing all these tuples has >= s - m + 1 members;
 (3) all s(s-1) ratios pi_j/pi_k lie in {gamma: h(gamma) <= log R_gamma, |gamma-1| <= R_gamma^(-1/2)} (C = 1);
 (4) 2-function proximity ratio -log max(|u-1|,|v-1|)/h(1:u:v) for u = pi_j/pi_1, v = pi_k/pi_1.
Run: python3 -I gp_cluster.py"""
import math, itertools, random, sys, time
from sympy import isprime
random.seed(7)

def find(s, L0):
    # modulus ~ exp(L0/2) * 1.0001 ; point (3t, t) direction, arg = atan(1/3) ~ 0.32 rad, inside (0, pi/4)
    t = int(math.exp(L0 / 2 - math.log(math.sqrt(10)))) + 12345 if L0 < 1400 else None
    if t is None:
        # big: use integer exponent
        k = int((L0 / 2 - math.log(math.sqrt(10))) / math.log(10)) + 1
        t = 10 ** k + 12345
    a0, b0 = 3 * t, t
    found = []
    side = 8
    seen = set()
    while len(found) < s:
        for a in range(a0, a0 + side):
            for b in range(b0, b0 + side):
                if (a, b) in seen:
                    continue
                seen.add((a, b))
                n = a * a + b * b
                if isprime(n):
                    found.append((a, b, n))
        side *= 2 if len(found) < s else 1
    # keep the s points closest to the corner (a0,b0) in sup-norm, to keep the square small
    found.sort(key=lambda p: max(p[0] - a0, p[1] - b0))
    pts = found[:s]
    return pts, a0, b0

def run(s, A, m_list=(2, 3, 4, 5)):
    t0 = time.time()
    L0 = A * s * math.log(s) + 5.0          # target log H (H = norm), so H >= exp(A s log s)
    pts, a0, b0 = find(s, L0)
    norms = [p[2] for p in pts]
    assert len(set(norms)) == s, "norms not distinct"
    assert all(p[2] % 4 == 1 for p in pts)
    side = max(max(p[0] - a0, p[1] - b0) for p in pts) + 1
    H = max(norms); Nmin = min(norms)
    D2 = max((p[0]-q[0])**2 + (p[1]-q[1])**2 for p in pts for q in pts)
    logH = math.log(H)
    print(f"s={s} A={A}: log H = {logH:.1f} >= A s log s = {A*s*math.log(s):.1f}; square side used = {side}; "
          f"max|pi_i-pi_j|^2 = {D2}; H/Nmin - 1 = {H/Nmin-1:.2e}")
    for m in m_list:
        if m > s: continue
        # rigorous worst case over all ordered m-tuples: LHS <= H * D2^(m-1) / H^m / Nmin^m
        logLHS = logH + (m - 1) * math.log(D2) - m * logH - m * math.log(Nmin)
        dmin = -logLHS / logH - m
        # exact evaluation on a sample of tuples
        worst = float('inf')
        for _ in range(2000):
            tup = random.sample(pts, m)
            x1 = tup[0]
            lg = math.log(x1[2]) + sum(math.log((p[0]-x1[0])**2 + (p[1]-x1[1])**2) for p in tup[1:]) - m * logH \
                 - sum(math.log(p[2]) for p in tup)
            worst = min(worst, -lg / logH - m)
        print(f"   m={m}: delta_eff >= {dmin:.3f} for ALL ordered m-tuples (needed: any delta < (m-1)/2 = {(m-1)/2}); "
              f"sampled min {worst:.3f}; ordered tuples {math.perm(s, m)}; proper subspaces needed >= s-m+1 = {s-m+1}")
    # near-1 set (C = 1)
    ok = 0
    for p in pts:
        for q in pts:
            if p is q: continue
            hg = 0.5 * math.log(max(p[2], q[2]))         # h(gamma) = log R_gamma
            dist = math.sqrt((p[0]-q[0])**2 + (p[1]-q[1])**2) / math.sqrt(q[2])
            if math.log(dist) <= -0.5 * hg:
                ok += 1
    print(f"   near-1 set: {ok} of s(s-1) = {s*(s-1)} ratios gamma = pi_j/pi_k satisfy |gamma-1| <= R_gamma^(-1/2), h(gamma) = log R_gamma")
    # two-function proximity ratio
    p1 = pts[0]; rmin = 9
    for p, q in itertools.permutations(pts[1:], 2):
        du = math.sqrt((p[0]-p1[0])**2 + (p[1]-p1[1])**2) / math.sqrt(p1[2])
        dv = math.sqrt((q[0]-p1[0])**2 + (q[1]-p1[1])**2) / math.sqrt(p1[2])
        h = 0.5 * math.log(max(p1[2], p[2], q[2]))
        rmin = min(rmin, -math.log(max(du, dv)) / h)
    print(f"   2-function proximity ratio lambda_inf/h for (pi_j/pi_1, pi_k/pi_1): min {rmin:.3f} (cluster value 1/2); time {time.time()-t0:.0f}s")
    return pts

if __name__ == "__main__":
    for s, A in [(12, 0.5), (40, 0.5), (40, 1.0), (100, 0.5)]:
        run(s, A)
        sys.stdout.flush()
