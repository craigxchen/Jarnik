# Pushforward of sub-cluster systems to the coordinates t of the global column space (Section 3 of subcluster.md).
# Exact dimensions (mod p ranks) + numerical base locus in t-space (EVIDENCE).
# Global profile: points x with exponent vectors b_x in {0..e}^rho; columns = coordinates (r = rho, e_p = e);
# delta_x = sum_i t_i b_{x,i}; pi_J(k) = (sum_{x in J} k_x b_{x,i})_i.
import sys, itertools, math, json, numpy as np
from common import *

def model_systems(pts, e, n, T, w, m, use_all_bases=True, max_J=None, rng=None, verbose=False):
    M = len(pts); rho = len(pts[0])
    cols = [tuple(pts[x][i] for x in range(M)) for i in range(rho)]
    Js = list(itertools.combinations(range(M), n))
    if max_J and len(Js) > max_J:
        idx = rng.choice(len(Js), max_J, replace=False); Js = [Js[i] for i in sorted(idx)]
    out = []
    for J in Js:
        for b in (J if use_all_bases else [J[-1]]):
            bi = J.index(b); others = [i for i in range(n) if i != bi]
            A = []
            for ks in itertools.product(range(-T, T+1), repeat=n-1):
                k = [0]*n
                for i, v in zip(others, ks): k[i] = v
                k[bi] = -sum(ks)
                y = tuple(sum(k[j]*pts[J[j]][i] for j in range(n)) for i in range(rho))
                # slab |<k,beta_p>| <= w with <k,beta_p> = (2/e) y_p
                if all(2*abs(yy) <= w*e for yy in y): A.append((tuple(k), y))
            if not A: continue
            hs = []
            for i in range(rho):
                v = [a[1][i] for a in A]; hs.append((max(v)-min(v))*2/e)
            LP = sum(hs)/rho
            if not m > LP:   # model admissibility: LP supercritical at leading order
                continue
            out.append(dict(J=J, base=b, A=A, LP=LP))
    return out

def poly_rows_mod_p(Ak, bi, m):
    n = len(Ak[0]); others = [i for i in range(n) if i != bi]
    mons = monomials(n-1, m)
    return np.array([[math.prod(k[o]**a for o, a in zip(others, al)) % P31 for k in Ak] for al in mons], dtype=np.int64)

def analyse(pts, e, n, T, w, m, label, max_J=None, seed=0, grid=0, alphas=True):
    rng = np.random.default_rng(seed)
    rho = len(pts[0]); M = len(pts)
    systems = model_systems(pts, e, n, T, w, m, max_J=max_J, rng=rng)
    s = (e*w)//2
    grid_pts = list(itertools.product(range(-s, s+1), repeat=rho)); gidx = {g: i for i, g in enumerate(grid_pts)}
    G = len(grid_pts)
    print(f"\n[{label}] M={M} rho={rho} e={e} n={n} T={T} w={w} m={m}: admissible (J,base) systems: {len(systems)}; y-box {2*s+1}^{rho}", flush=True)
    Xs = []; dims = []
    for sy in systems:
        Ak = [a[0] for a in sy['A']]; Ay = [a[1] for a in sy['A']]
        bi = sy['J'].index(sy['base'])
        Pm = poly_rows_mod_p(Ak, bi, m)
        V = nullspace_mod_p(Pm)               # rows: basis of V(A,m) mod p
        Pi = np.zeros((len(Ak), G), dtype=np.int64)
        for j, y in enumerate(Ay): Pi[j, gidx[y]] = 1
        X = (V @ Pi) % P31 if len(V) else np.zeros((0, G), dtype=np.int64)
        rX = rank_rows_mod_p(X) if len(X) else 0
        S = sorted(set(Ay))
        sy.update(dimV=len(V), dimX=rX, S=S, X=X)
        dims.append((len(V), rX, len(S), len(Ak)))
    Xall = np.concatenate([sy['X'] for sy in systems if len(sy['X'])], axis=0) if systems else np.zeros((0, G))
    R, piv = rref_mod_p(Xall) if len(Xall) else (np.zeros((0, G), dtype=np.int64), [])
    dim_comb = len(piv)
    S_tot = sorted(set(y for sy in systems for y in sy['S']))
    # full order-m space on the full y-box, and on S_tot
    Pbox = np.array([[math.prod(g[i]**al[i] for i in range(rho)) % P31 for g in grid_pts] for al in monomials(rho, m)], dtype=np.int64)
    dim_Vbox = G - rank_rows_mod_p(Pbox)
    single = [d[1] for d in dims]
    print(f"  single systems: dim V(A,m) range {min(d[0] for d in dims)}..{max(d[0] for d in dims)}, "
          f"|A| range {min(d[3] for d in dims)}..{max(d[3] for d in dims)}, |S| range {min(d[2] for d in dims)}..{max(d[2] for d in dims)}")
    print(f"  dim pi_*V per system: min {min(single)} max {max(single)}; number with pi_*V = 0: {sum(1 for x in single if x == 0)}")
    print(f"  combined dim sum_J pi_*V = {dim_comb};  |S_tot| = {len(S_tot)};  dim V(y-box, m) = {dim_Vbox}", flush=True)
    res = dict(label=label, M=M, rho=rho, n=n, T=T, w=w, m=m, nsys=len(systems), single_max=max(single), dim_comb=dim_comb, dim_Vbox=dim_Vbox)
    if alphas:
        # product differences Delta^alpha (anchored at the lowest corner) with |alpha| = m inside the y-box
        found_single = 0; found_comb = 0; tested = 0
        Rset = R
        for al in itertools.product(range(0, 2*s+1), repeat=rho):
            if sum(al) != m: continue
            tested += 1
            vec = np.zeros(G, dtype=np.int64)
            for be in itertools.product(*[range(a+1) for a in al]):
                y = tuple(-s + bb for bb in be)
                coef = (-1)**(sum(al)-sum(be)) * math.prod(math.comb(a, bb) for a, bb in zip(al, be))
                vec[gidx[y]] = coef % P31
            in_comb = rank_rows_mod_p(np.vstack([Rset, vec[None, :]])) == dim_comb
            found_comb += in_comb
            fs = 0
            for sy in systems:
                if sy['dimX'] and rank_rows_mod_p(np.vstack([sy['X'], vec[None, :]])) == sy['dimX']: fs = 1; break
            found_single += fs
        print(f"  product differences Delta^alpha, |alpha|=m, anchored at corner: tested {tested}; in some single pi_*V: {found_single}; in combined: {found_comb}", flush=True)
        res.update(alpha_tested=tested, alpha_single=found_single, alpha_comb=found_comb)
    return res, systems, grid_pts

if __name__ == "__main__":
    cube3 = [tuple(b) for b in itertools.product([0, 1], repeat=3)]
    for (n, T, w, m) in [(5, 2, 2, 5), (5, 2, 2, 6), (6, 2, 2, 5), (6, 2, 2, 6), (5, 3, 2, 5), (5, 3, 2, 6)]:
        analyse(cube3, 1, n, T, w, m, "cube rho=3")
