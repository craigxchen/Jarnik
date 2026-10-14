r"""Complete-block images Y_t(j) (loctheory §5) and which products of linear forms lie in J(Y_t(j))^perp.

Y_t(j) = { y = M c : c in Z^{2^t}, sum c = 0, |y_a| <= j for a in F_2^t \ 0 },  y_a = sum_{z : a.z = 1} c_z.
By loctheory Lemma 5.1 this is the image pi_F(A) of the complete t-block F = U \ 0 of Sylvester columns for every large n.

A restricted generator is a_i = M c_i, c_i the class-sum vector of z_i (cosets of U^perp).  Types (t = 3):
  light : c = e_z - e_z'                      (z_i a pair, or a coset difference with H subset U^perp)
  heavy2: c = 1_{x'+H'} - 1_{y'+H'}, |H'| = 2   (coset difference with |H : H cap U^perp| = 2)
  heavy4: c = 1_{x'+H'} - 1_{y'+H'}, |H'| = 4
usage: python3 yblock.py t j emax
Prints H_Y, rho, alpha, and for each degree e <= emax the maximal number of heavy2 generators (and of light
generators) whose product lies in J(Y)^perp_e, by exhaustive search over multisets subject to the (necessary)
per-axis grid count <= 2j.  Membership: mod-p rank test against T_e = {mu(c) : c in V(Y,e)} (a negative is rigorous:
T_e(mod p) contains the reduction of T_e(Q)); positives are re-checked exactly.
"""
import sys, itertools, time
import numpy as np
from math import factorial
from flint import nmod_mat, fmpz_mat, fmpq_mat, fmpq

PR = 2147475331


def build_Y(t, j):
    Z = list(range(2 ** t))
    axes = list(range(1, 2 ** t))
    M = np.array([[bin(a & z).count('1') % 2 for z in Z] for a in axes], dtype=np.int64)
    R = 2 * j - 1
    pts = set()
    # enumerate c with sum 0, |c_z| <= R (|c_z| < 2j by loctheory Lemma 5.1 proof)
    rng = np.arange(-R, R + 1)
    grids = np.array(np.meshgrid(*([rng] * (2 ** t - 1)), indexing='ij')).reshape(2 ** t - 1, -1).T
    last = -grids.sum(axis=1)
    ok = np.abs(last) <= R
    C = np.hstack([grids[ok], last[ok][:, None]])
    Y = C @ M.T
    good = (np.abs(Y) <= j).all(axis=1)
    Y = np.unique(Y[good], axis=0)
    return M, Y


def monos(b, e):
    out = []
    for c in itertools.combinations_with_replacement(range(b), e):
        a = [0] * b
        for x in c:
            a[x] += 1
        out.append(tuple(a))
    return out


def evalm(Y, mons):
    E = np.ones((len(mons), len(Y)), dtype=object)
    Yo = Y.astype(object)
    for r, a in enumerate(mons):
        col = np.ones(len(Y), dtype=object)
        for x, ax in enumerate(a):
            if ax:
                col = col * Yo[:, x] ** ax
        E[r] = col
    return E


def to_nmod(E, p=PR):
    r, c = E.shape
    return nmod_mat(r, c, [int(x) % p for x in E.flatten()], p)


def main():
    t, j, emax = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    M, Y = build_Y(t, j)
    b = Y.shape[1]
    print(f't={t} j={j} |Y|={len(Y)} b={b}', flush=True)
    # Hilbert function mod p (upper bound for rho is rigorous; H itself mod p)
    cum = 0
    H = []
    low = []
    for e in range(0, 4 * b * j + 1):
        low_e = low + monos(b, e)
        r = to_nmod(evalm(Y, low_e)).rank()
        H.append(r - cum)
        cum = r
        low = low_e
        if cum == len(Y):
            break
    print('H_Y (mod p) =', H, ' rho <=', len(H) - 1, flush=True)
    grid = [1]
    for _ in range(b):
        new = [0] * (len(grid) + 2 * j)
        for i, g in enumerate(grid):
            for k in range(2 * j + 1):
                new[i + k] += g
        grid = new
    alpha = next((e for e in range(len(H)) if H[e] != grid[e]), None)
    print('H_grid =', grid[:len(H) + 2], ' alpha (first e with H_Y < H_grid) =', alpha, flush=True)
    # generator types
    Zs = list(range(2 ** t))
    gens = {'light': [], 'heavy2': [], 'heavy4': []}
    for z1, z2 in itertools.combinations(Zs, 2):
        c = np.zeros(2 ** t, dtype=np.int64)
        c[z1], c[z2] = 1, -1
        gens['light'].append(M @ c)
    seen = set()
    for h in range(1, 2 ** t):
        Hs = [0, h]
        cosets = sorted({tuple(sorted({x ^ y for y in Hs})) for x in Zs})
        for C1, C2 in itertools.combinations(cosets, 2):
            c = np.zeros(2 ** t, dtype=np.int64)
            for x in C1:
                c[x] += 1
            for x in C2:
                c[x] -= 1
            a = tuple(M @ c)
            if a not in seen and tuple(-x for x in a) not in seen:
                seen.add(a)
                gens['heavy2'].append(np.array(a))
    print('#light', len(gens['light']), '#heavy2', len(gens['heavy2']), flush=True)
    print('heavy2 example restricted vector a =', gens['heavy2'][0], flush=True)
    rng = np.random.default_rng(1)
    nsamp = int(sys.argv[4]) if len(sys.argv) > 4 else 3000
    for e in range(max(2, (alpha or 2) - 1), emax + 1):
        lowm = []
        for f in range(e):
            lowm += monos(b, f)
        N = to_nmod(evalm(Y, lowm))
        X, nul = N.nullspace()
        me = monos(b, e)
        Em = to_nmod(evalm(Y, me).T)  # |Y| x N_e
        if nul == 0:
            print(f'e={e}: V = 0', flush=True)
            continue
        Vm = nmod_mat(nul, X.nrows(), [int(X[i, k]) for k in range(nul) for i in range(X.nrows())], PR)
        T = (Vm * Em).rref()[0]
        rT = T.rank()
        Tl = [[int(T[i, c]) for c in range(T.ncols())] for i in range(rT)]
        piv = [next(c for c in range(len(me)) if Tl[i][c]) for i in range(rT)]
        Ta = np.array(Tl, dtype=np.int64)
        idx = {a: i for i, a in enumerate(me)}
        mult = [factorial(e) // int(np.prod([factorial(int(x)) for x in a])) for a in me]
        minv = [pow(int(q) % PR, -1, PR) for q in mult]

        def member(G):
            poly = {tuple([0] * b): 1}
            for g in G:
                new = {}
                for ex, cf in poly.items():
                    for x in range(b):
                        if g[x]:
                            e2 = list(ex)
                            e2[x] += 1
                            e2 = tuple(e2)
                            new[e2] = new.get(e2, 0) + cf * int(g[x])
                poly = {k: v for k, v in new.items() if v}
            vec = np.zeros(len(me), dtype=np.int64)
            for ex, cf in poly.items():
                vec[idx[ex]] = (cf % PR) * minv[idx[ex]] % PR
            for i in range(rT):
                c = int(vec[piv[i]])
                if c:
                    vec = (vec - (c * Ta[i]) % PR) % PR
            return not vec.any()

        lines = [(a, c, a ^ c) for a in range(1, 2 ** t) for c in range(a + 1, 2 ** t) if (a ^ c) > c]
        for typ in ('light', 'heavy2'):
            G = gens[typ]
            nz = [np.nonzero(g)[0] for g in G]
            stats = {}
            for _ in range(nsamp):
                # random multiset with axis counts <= 2j, built greedily
                cnt = np.zeros(b, dtype=np.int64)
                ch = []
                tries = 0
                while len(ch) < e and tries < 200:
                    i = int(rng.integers(len(G)))
                    tries += 1
                    if (cnt[nz[i]] + 1 <= 2 * j).all():
                        cnt[nz[i]] += 1
                        ch.append(i)
                if len(ch) < e:
                    continue
                mx = int(cnt.max())
                # weighted axis cost: sum over generators of |entry| / (gcd of entries) on each axis
                wcost = np.zeros(b, dtype=np.int64)
                for i in ch:
                    g = G[i]
                    gg = int(np.gcd.reduce(np.abs(g[g != 0])))
                    wcost += np.abs(g) // gg if typ == 'light' else np.abs(g)
                wmx = int(wcost.max())
                nline = max(sum(1 for i in ch if any(G[i][q - 1] for q in L)) for L in lines) if t >= 2 else 0
                parity_ok = nline <= 3 * j - 1
                mem = member([G[i] for i in ch])
                key = (mx, wmx, parity_ok)
                a, t_ = stats.get(key, (0, 0))
                stats[key] = (a + int(mem), t_ + 1)
            print(f'e={e} dimT={rT} {typ}: (max count, max weighted cost, lines<=3j-1) -> (in J^perp, sampled):', flush=True)
            for k, v in sorted(stats.items()):
                print(f'      {k}: {v}', flush=True)


if __name__ == '__main__':
    main()
