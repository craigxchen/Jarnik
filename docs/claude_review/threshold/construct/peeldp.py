"""Rigorous upper bound for reg of symmetric cut polytopes by sorted coordinate peeling.

Pi(M, U, s) = { k in Z^M : sum k = s, k(B) <= U[|B|-1] for all proper nonempty B }.
Slicing k_1 = c gives Pi(M-1, U', s-c) with U'[b-1] = min(U[b-1], U[b] - c), b = 1..M-2.
Peeling lemma (upper.md Lemma 2.1, re-verified in construct.md):
    reg(A) <= max_i ( reg(A_(i)) + i - 1 ),  slices sorted by decreasing (bound on) reg.
Integer feasibility: Pi(M,U,s) has an integer point iff the most balanced integer vector with sum s
(entries floor/ceil of s/M) satisfies all top-b constraints (it is majorization-minimal)."""
import sys
from functools import lru_cache

def feasible(M, U, s):
    q, r = divmod(s, M)
    # balanced vector: r entries q+1, M-r entries q ; top-b sum
    for b in range(1, M):
        top = b * q + min(b, r)
        if top > U[b - 1]:
            return False
    return True

def tighten(M, U, s):
    """replace U by the true maxima of top-b sums over the real polytope's integer points, using
    only valid implications: top_b <= U_b ; top_b <= top_a + top_{b-a} ; top_b <= s - (M-b)*min_entry ...
    keep it simple: subadditive closure."""
    U = list(U)
    changed = True
    while changed:
        changed = False
        for b in range(2, M):
            for a in range(1, b):
                v = U[a - 1] + U[b - a - 1]
                if v < U[b - 1]:
                    U[b - 1] = v; changed = True
    return tuple(U)

@lru_cache(maxsize=None)
def bound(M, U, s):
    if M == 1:
        return 0
    U = tighten(M, U, s)
    lo = s - U[M - 2]; hi = U[0]
    rs = []
    for c in range(lo, hi + 1):
        U2 = tuple(min(U[b - 1], U[b] - c) for b in range(1, M - 1))
        s2 = s - c
        if M - 1 == 1:
            # single coordinate k_2 = s2 ; constraint set empty (no proper subsets)
            rs.append(0)
            continue
        if not feasible(M - 1, U2, s2):
            continue
        rs.append(bound(M - 1, U2, s2))
    rs.sort(reverse=True)
    return max(r + i for i, r in enumerate(rs))

if __name__ == "__main__":
    M = int(sys.argv[1])
    prof = [int(x) for x in sys.argv[2].split(",")]   # u_1..u_{M-1}
    for n in range(1, int(sys.argv[3]) + 1):
        U = tuple(n * u for u in prof)
        b = bound(M, U, 0)
        print(f"M={M} profile={prof} n={n}: peeling bound {b}  per n {b/n:.3f}", flush=True)
