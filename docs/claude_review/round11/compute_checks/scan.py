"""Screen full-rank sign profiles for the supercritical window.
rho(A) computed by rank mod p = 2^61-1 (rank_p <= rank_Q, so rho_Q <= rho_p; exact confirmation later).
rmax = max_p range_A(k(S_p))  (so h_A(beta_p) = 2 range).  (ii) needs 2*rmax < m <= rho; (ii') rmax < m <= rho.
Also prints seg = longest lattice segment length (N-1) in A over primitive directions with |v|_inf <= 2,
an upper bound proxy for 1-dim sections."""
import itertools, sys, random
import flint
from core import build_A, monomials_lt, mono_val, width
from models import fullrank, canon, random_profile
P = 2**61 - 1

def rho_p(A):
    U = len(A); j = 0
    while True:
        rows = monomials_lt(len(A[0]), j + 1)
        M = flint.nmod_mat([[mono_val(k, a) % P for k in A] for a in rows], P)
        if M.rank() == U:
            return j
        j += 1

def pats(n):
    out = []
    for s in range(1, n // 2 + 1):
        for S in itertools.combinations(range(n), s):
            out.append(canon(n, S))
    return sorted(set(out))

if __name__ == "__main__":
    n = int(sys.argv[1]); nprof = int(sys.argv[2]); Tmax = int(sys.argv[3])
    rng = random.Random(3)
    allp = pats(n)
    profs = []
    if n == 4:
        profs = [p for r in range(3, len(allp) + 1) for p in itertools.combinations(allp, r) if fullrank(n, p)]
    else:
        while len(profs) < nprof:
            p = random_profile(n, rng)
            if fullrank(n, p) and p not in profs:
                profs.append(p)
    for prof in profs:
        for L in (1, 2, 3):
            for T in range(1, Tmax + 1):
                if (2 * T + 1) ** (n - 1) > 3500:
                    continue
                A = build_A(n, T, [(S, L) for S in prof])
                if len(A) == (2 * T + 1) ** (n - 1):
                    continue
                rmax = max(width(A, S) for S in prof)
                r = rho_p(A)
                tag = "II" if r > 2 * rmax else ("II'" if r > rmax else "-")
                if tag != "-":
                    print(n, prof, "T=%d L=%d |A|=%d rmax=%d rho=%d %s" % (T, L, len(A), rmax, r, tag), flush=True)
