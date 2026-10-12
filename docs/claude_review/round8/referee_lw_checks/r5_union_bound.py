"""Referee check R5: arithmetic of lw.md Theorem 6.1, Corollary 3.3, and the size of the representation count.

1. Theorem 6.1(a): 4.3 * sum_{n>=2} (n + pi/C0) (C0/20)^n < 1/16 for C0 in (0,1].
2. Theorem 6.1(c), character part, evaluated with the proof's own bounds at P' = P_1 (log-space):
   (4M/C0)(1+1/eps) e^{-eps W/8} * 2Z/(1-2Z)^2, Z = Z(eps/2) <= 4r/(P'^{eps/2}-1), W >= 2q log P'.
3. Alternative split (1+eps)h' >= W/4 + x0 + (1-theta) eps W/4 + theta eps h': the smallest log P' that makes the
   character part < 1/8, as a multiple of log M (shows the proof's 4/eps can be lowered to about 2/eps;
   the 1/eps itself is intrinsic to summing over representations).
4. Corollary 3.3: the proof's claim  E <= W/(2 log M)  versus the actual size of the error term.
5. Bell numbers: log of the number of groupings of a support of size M/2+1, versus (M/2) log M.
Usage: python3 r5_union_bound.py
"""
import math
out = []
def log(s): print(s); out.append(s)

# 1
for C0 in (1.0, 0.5, 0.1, 0.01):
    S = 4.3 * sum((n + math.pi / C0) * (C0 / 20) ** n for n in range(2, 200))
    log(f"(a) C0={C0}: 4.3 sum (n+pi/C0)(C0/20)^n = {S:.5f}  (< 1/16 = 0.0625: {S < 0.0625})")

def logZbound(s, r, logP):
    # log of 4r/(P^s - 1), valid when 2r <= P^s - 1
    Ps = math.exp(min(s * logP, 700))
    if 2 * r > Ps - 1: return None
    return math.log(4 * r / (Ps - 1))

# 2
for eps in (1.0, 0.5, 0.25, 0.1, 0.05):
    for C0 in (1.0,):
        Mmin = None
        for q in [7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83, 103, 107, 127, 131, 139, 151, 163, 167, 179, 191, 199, 211, 223, 227, 239, 251]:
            M = q + 1; r = 2 * q
            logP = max(4 * math.log(40 * M / C0), (2 / eps) * math.log(64 * M))
            W = 2 * q * logP
            lz = logZbound(eps / 2, r, logP)
            if lz is None: continue
            twoZ = 2 * math.exp(lz)
            if twoZ > 0.25: continue
            val = math.log(4 * M / C0) + math.log(1 + 1 / eps) - eps * W / 8 + math.log(twoZ / (1 - twoZ) ** 2)
            if val < math.log(1 / 8):
                Mmin = M; break
        log(f"(c) eps={eps}, C0={C0}: character part < 1/8 at P'=P_1 already for M = {Mmin}; "
            f"W/(M log M) at that M ~ {2*q*logP/(M*math.log(M)):.1f} (asymptotic max(8,4/eps) = {max(8,4/eps)})")

# 3: alternative split.  Character part <= (4M/C0)(1+1/eps) e^{-x0} e^{-(1-theta) eps W/4} sum_K K (2 Z(theta eps))^K
def min_logP(M, eps, C0, theta):
    q = M - 1; r = 2 * q
    lo, hi = 1.0, 1e5
    for _ in range(200):
        mid = (lo + hi) / 2
        lz = logZbound(theta * eps, r, mid)
        ok = False
        if lz is not None:
            twoZ = 2 * math.exp(lz)
            if twoZ < 1:
                W = 2 * q * mid
                val = (math.log(4 * M / C0) + math.log(1 + 1 / eps) - (mid / 2 - 1 / 16)
                       - (1 - theta) * eps * W / 4 + math.log(twoZ / (1 - twoZ) ** 2))
                ok = val < math.log(1 / 8)
        if ok: hi = mid
        else: lo = mid
    return hi
for eps in (0.25, 0.1, 0.05):
    for M in (10**3, 10**5, 10**7):
        best = min((min_logP(M, eps, 1.0, th), th) for th in (0.5, 0.8, 0.9, 0.95, 0.99))
        log(f"(3) eps={eps}, M={M}: smallest log P' for the character part (best theta={best[1]}) = "
            f"{best[0]/math.log(M):.3f} log M  -> W/(M log M) ~ {2*best[0]/math.log(M):.2f} "
            f"(proof's choice: 4/eps = {4/eps:.0f}; 2/eps = {2/eps:.0f})")

# 4: Corollary 3.3
for M in (10**4, 10**6, 10**9):
    L = math.log(M); W = 3 * M * L
    Kp1 = W / L**3
    err = Kp1 * (L**2 + 3 * L + 9) + L**2
    log(f"(Cor 3.3) M={M}: error term with log(1/C(eps)) = 0 is {err/(W/L):.3f} * W/log M; "
        f"the proof claims <= 0.5 * W/log M: {'TRUE' if err <= W/(2*L) else 'FALSE'}; "
        f"the conclusion needs <= (1+eps) W/log M")

# 5: Bell numbers
def logbell(n):
    # Bell numbers via the Bell triangle (exact big integers)
    row = [1]
    for _ in range(n - 1):
        nr = [row[-1]]
        for x in row: nr.append(nr[-1] + x)
        row = nr
    return math.log(row[-1])
for M in (100, 400, 1000):
    m = M // 2 + 1
    log(f"(Bell) M={M}: log Bell(M/2+1) = {logbell(m):.1f}; (M/2) log M = {(M/2)*math.log(M):.1f}; "
        f"(M/2)(log M - loglog M) = {(M/2)*(math.log(M)-math.log(math.log(M))):.1f}")
open('r5_union_bound_out.txt', 'w').write('\n'.join(out) + '\n')
