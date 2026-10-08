# aux_paley_design.py -- abstract feasibility of the complete per-pair larger-sieve system at W = O(M log M).
# Exponent code: b copies of the q non-constant columns of the normalized Paley-Hadamard matrix of order
# M = q+1 (q = 3 mod 4 prime), placed on r = b q DISTINCT ACTUAL split primes just above M^4.
# Residue design: point x (labelled 1..M) gets, for every prime power l^a, the class
#     inert l, l=2 :  x mod l^a                                 (<= nu(l^a) classes)
#     split l      :  x mod l^a if l does not divide x, else (x - 1/2) mod l^a   (lands in the
#                     (l-1) l^(a-1) residues prime to l; refines across a)
# Every pair (x,y) then collides exactly at levels l^a dividing n, 2n-1 or 2n+1 (split, mixed type),
# n = |x-y|, so its sieve demand  sum_{colliding levels} 2 log l  is at most 2log n + 2log(2n+1) + 2log(2n-1).
# We verify EXACTLY (rational primes, exact class maps) that demand_xy <= d_xy - W/2 + 2 log C for C = 1,
# where d_xy is the actual weighted conductor distance.  This is an abstract model (no lattice points are
# claimed to lie on a short arc); it shows that no inequality derived only from per-pair collision data,
# cut structure and class-count constraints can force W >> M log M.
import math
from aux_gauss import is_prime

def legendre(a, q):
    a %= q
    return 0 if a == 0 else (1 if pow(a, (q - 1) // 2, q) == 1 else -1)

def paley(q):
    M = q + 1
    H = [[1] * M]                 # row 'infinity'
    for x in range(q):
        row = [1]
        for j in range(q):
            row.append(-1 if x == j else -legendre(x - j, q))
        H.append(row)
    for i in range(M):
        for k in range(i + 1, M):
            assert sum(a * c for a, c in zip(H[i], H[k])) == 0
    return H

def split_primes_from(X, r):
    out = []; p = X
    while len(out) < r:
        if p % 4 == 1 and is_prime(p): out.append(p)
        p += 1
    return out

def class_of(x, l, a):
    m = l ** a
    if l == 2 or l % 4 == 3: return x % m
    if x % l: return x % m
    inv2 = pow(2, -1, m)
    return (x - inv2) % m

def nu(l, a):
    if l == 2: return 2 if a == 1 else 2 ** (a + 1)
    return (l + 1) * l ** (a - 1) if l % 4 == 3 else (l - 1) * l ** (a - 1)

def run(q, b=5, C=1.0):
    H = paley(q); M = q + 1
    cols = [[H[i][j] for i in range(M)] for j in range(1, M)] * b
    r = len(cols)
    primes = split_primes_from(M ** 4, r)
    w = [math.log(p) for p in primes]
    W = sum(w)
    # prime powers that can matter: l^a <= 4M+2 (beyond that no pair difference/2n+-1 is divisible)
    levels = [(l, a) for l in range(2, 4 * M + 3) if is_prime(l) for a in range(1, 40) if l ** a <= 4 * M + 2]
    # class-count and refinement checks
    for (l, a) in levels:
        cls = {class_of(x, l, a) for x in range(1, M + 1)}
        assert len(cls) <= nu(l, a)
        if l % 4 == 1:
            assert all(c % l for c in cls)
        if (l, a + 1) in levels:
            for x in range(1, M + 1):
                assert class_of(x, l, a + 1) % (l ** a) == class_of(x, l, a)
    worst = 1e9; worst_line = 1e9
    for x in range(M):
        for y in range(x + 1, M):
            d = sum(wk for wk, col in zip(w, cols) if col[x] != col[y])
            excess = d - W / 2 + 2 * math.log(C)
            demand = sum(2 * math.log(l) for (l, a) in levels if class_of(x + 1, l, a) == class_of(y + 1, l, a))
            n = y - x
            assert demand <= 2 * math.log(n) + 2 * math.log(2 * n + 1) + 2 * math.log(2 * n - 1) + 1e-9
            worst = min(worst, excess - demand)
            # line-metric (ordered arc) model: positions x*C/M on an arc of normalized length C;
            # gamma := log of (chord^2 / gcd norm) = (d - W/2) + 2 log(C n / M) must lie in [demand, excess]
            gamma = (d - W / 2) + 2 * math.log(C * n / M)
            assert gamma <= excess + 1e-12
            worst_line = min(worst_line, gamma - demand)
    return M, r, W, W / (M * math.log(M)), worst, worst_line

if __name__ == "__main__":
    for q in (11, 19, 23, 31, 43):
        M, r, W, ratio, margin, mline = run(q)
        print("q=%d M=%d primes=%d W=%.1f W/(M log M)=%.2f  min(excess-demand)=%.3f  min(gamma_line-demand)=%.3f %s"
              % (q, M, r, W, ratio, margin, mline, "OK" if min(margin, mline) >= 0 else "FAIL"))
