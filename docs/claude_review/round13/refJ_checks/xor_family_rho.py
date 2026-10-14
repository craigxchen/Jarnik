"""rho of the image of a Sylvester box under an XOR-closed family of columns (all nonzero a in a t-dim subspace V):
Y_t(j) = {kappa in Z^{2^t} : sum kappa = 0, |kappa(S_a)| <= j for all a in F_2^t \\ 0} (cell sums; kappa(S_a) = sum over
cells c with a.c = 0).  For large n every such kappa is attained (cells have n/2^t points, entries of k in {-1,0,1}),
so Y_t(j) is the image.  Coordinates: drop kappa_0 (= -sum others).  Hilbert function via mod-p ranks (EVIDENCE) and
exact fmpz ranks where small."""
import itertools, sys
from flint import nmod_mat, fmpz_mat
PR = 2**31 - 1
def dot(a, x): return bin(a & x).count('1') & 1
def Y(t, j):
    N = 2**t; B = 7 * j  # crude bound on |kappa_x|
    hyper = [[c for c in range(N) if dot(a, c) == 0] for a in range(1, N)]
    # enumerate via Fourier: kappa_x = (1/N) sum_a khat(a) (-1)^{a.x}, khat(a) = 2 kappa(S_a) in [-2j, 2j] even
    pts = []
    for ks in itertools.product(range(-j, j + 1), repeat=N - 1):
        kh = [0] + [2 * s for s in ks]
        kap = []
        ok = True
        for x in range(N):
            v = sum(kh[a] * (-1) ** dot(a, x) for a in range(N))
            if v % N: ok = False; break
            kap.append(v // N)
        if ok: pts.append(tuple(kap[1:]))
    return pts
def monos(nv, i):
    if nv == 0: return [()] if i == 0 else []
    return [(a,) + m for a in range(i + 1) for m in monos(nv - 1, i - a)]
t = int(sys.argv[1]); out = []
for j in range(1, int(sys.argv[2]) + 1):
    P = Y(t, j); nv = 2**t - 1
    rows = []; ranks = []; i = 0
    while True:
        for e in monos(nv, i):
            r = []
            for p in P:
                v = 1
                for z in range(nv): v *= p[z] ** e[z]
                r.append(v % PR)
            rows.append(r)
        rk = nmod_mat(len(rows), len(P), [x for r in rows for x in r], PR).rank()
        ranks.append(rk)
        if rk == len(P): break
        i += 1
    s = f"t={t} j={j}: |Y|={len(P)} rho_p={len(ranks)-1} (rho_p <= rho_Q only if ranks agree; evidence)  ratio rho/j={(len(ranks)-1)/j:.3f}"
    print(s, flush=True); out.append(s)
open(f'xor_family_rho_t{t}_out.txt', 'w').write("\n".join(out) + "\n")
