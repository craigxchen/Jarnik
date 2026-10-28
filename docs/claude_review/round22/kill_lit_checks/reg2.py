# Independent recomputation of the maximal vanishing order at e of a section of m(-K) on Lbar_k^tw.
# Sections over Qbar: Laurent polys with support in S_m = mP cap M, P = {u in R^k: sum u = 0, <u,1_A> >= -1 all A}.
# ord_e(sum c_u x^u) >= r  iff  sum_u c_u q(u) = 0 for every polynomial q of degree < r  (x = exp(t)).
# max ord = largest r such that the evaluation map Poly_{<r} -> F^{S} is NOT surjective  (i.e. rank < |S|).
# Differences from lit_checks/reg.py: lattice points enumerated from the full facet list (not ||u||_1 <= 2m),
# coordinates = LAST k-1 entries, basis = binomial polynomials prod binom(u_i + m, a_i) (integer valued), ranks by
# python-flint nmod_mat with two different primes, plus exact fmpz rank where cheap.
import itertools, sys
from math import comb
import flint

def lattice_points(k, m):
    A_list = [A for r in range(1, k) for A in itertools.combinations(range(k), r)]
    pts = []
    for u in itertools.product(range(-m, m + 1), repeat=k - 1):
        last = -sum(u)
        v = u + (last,)
        if abs(last) > m:
            continue
        if all(sum(v[i] for i in A) >= -m for A in A_list):
            pts.append(v)
    return pts

def exps(n, d):  # all exponent vectors of total degree <= d in n variables
    out = []
    for tot in range(d + 1):
        for c in itertools.combinations(range(tot + n - 1), n - 1):
            prev = -1; e = []
            for x in c:
                e.append(x - prev - 1); prev = x
            e.append(tot + n - 1 - prev - 1)
            out.append(tuple(e))
    return out

def rank_mod(S, m, d, p, coords):
    n = len(coords)
    E = exps(n, d)
    rows = []
    for u in S:
        x = [u[c] + m for c in coords]  # shift to nonnegative: binomial basis is still a basis of Poly_{<=d}
        rows.append([ (lambda e: _prod([comb(x[i], e[i]) for i in range(n)]) % p)(e) for e in E])
    M = flint.nmod_mat(rows, p)
    return M.rank(), len(E)

def _prod(l):
    r = 1
    for a in l: r *= a
    return r

def max_ord(k, m, p, coords, start=0):
    S = lattice_points(k, m)
    d = start
    while True:
        rk, nE = rank_mod(S, m, d, p, coords)
        if rk == len(S):
            return len(S), d   # Poly_{<=d} surjective => no section with ord >= d+1; Poly_{<=d-1} not (checked by loop)
        d += 1

if __name__ == '__main__':
    k = int(sys.argv[1]); mmax = int(sys.argv[2])
    for p in (1000003, 2305843009213693951 % (2**61-1) or 998244353):
        pass
    primes = [1000003, 998244353]
    for m in range(1, mmax + 1):
        res = []
        for p in primes:
            nS, d = max_ord(k, m, p, coords=list(range(1, k)))
            res.append(d)
        print(f"k={k} m={m} |S|={nS} max_ord(mod p, p in {primes}) = {res}  ratio {res[0]/m:.4f}", flush=True)
