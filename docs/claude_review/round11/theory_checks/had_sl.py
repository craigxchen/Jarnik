"""SL-full test on the widened slab box of a Sylvester-Hadamard profile (evidence, far below the regime).

Points x in F_2^k (n = 2^k); columns beta_q(x) = (-1)^{q.x}, q != 0 (balanced, orthogonal, Plotkin tight).
Slice coordinates k' in [-T,T]^{n-1}, k_n = -sum k'  (base point = last point).
A_p = {k : |k_x| <= T (x<n), |<k,beta_q>| <= w (q != p), |<k,beta_p>| <= 2*m}   (e_p = 1, l_p = <k,beta_p>/2).
For each candidate k* (sampled), S = {k*} cup {k in A_p : l(k) > l(k*)}:
  char0: is e_{k*} hit, i.e. rank Ev_S == rank Ev_{S \\ k*} + 0   (some c in V(S,m) with c_{k*} != 0)
  pint : col_{k*} in Z_(p)-span of the other columns (c_{k*} = 1, others p-integral  =>  UM_p)
Ranks mod 2^31-1 (lower bounds for rational ranks).  Polynomials: binomial basis prod binom(k_x+T, a_x).
"""
import itertools, random, sys
from math import comb
import numpy as np
from bool_sl import rank_modP, member_Zp

def profile(kdim):
    pts = list(itertools.product((0, 1), repeat=kdim))
    cols = []
    for q in pts[1:]:
        cols.append([(-1) ** (sum(a * b for a, b in zip(q, x)) % 2) for x in pts])
    return pts, cols

def build(kdim, T, w, m, p_idx):
    pts, cols = profile(kdim)
    n = len(pts)
    d = n - 1
    B = np.array(cols)            # r x n
    Bp = B[:, :d] - B[:, [d]]     # slice coefficients
    box = np.array(list(itertools.product(range(-T, T + 1), repeat=d)), dtype=np.int64)
    vals = box @ Bp.T             # |box| x r : <k, beta_q>
    ok = np.ones(len(box), dtype=bool)
    for q in range(B.shape[0]):
        if q == p_idx:
            ok &= np.abs(vals[:, q]) <= 2 * m
        else:
            ok &= np.abs(vals[:, q]) <= w
    A = box[ok]
    lev = vals[ok][:, p_idx] // 2
    return A, lev, d

def monos(d, m):
    out = []
    for deg in range(m):
        for a in itertools.product(range(deg + 1), repeat=d):
            if sum(a) == deg:
                out.append(a)
    return out

def evmat(S, mon, T):
    # rows = monomials, cols = points
    G = np.zeros((len(mon), len(S)), dtype=np.int64)
    P = 2147483647
    for j, a in enumerate(mon):
        v = np.ones(len(S), dtype=object)
        for x, ax in enumerate(a):
            if ax:
                v = v * np.array([comb(int(t) + T, ax) for t in S[:, x]], dtype=object)
        G[j] = np.array([int(t) % P for t in v], dtype=np.int64)
    return G

if __name__ == '__main__':
    kdim, T, w, m, p, nsamp = (int(a) for a in sys.argv[1:7])
    rng = random.Random(1)
    A, lev, d = build(kdim, T, w, m, 0)
    mon = monos(d, m)
    print(f"n={2**kdim} d={d} T={T} w={w} m={m}: |A_p|={len(A)} E={len(mon)} levels {lev.min()}..{lev.max()}")
    Gall = evmat(A, mon, T)
    rA = rank_modP(Gall)
    print(f"  rank Ev_A = {rA}, dim V(A_p,m) >= {len(A) - rA}")
    idx = list(range(len(A)))
    rng.shuffle(idx)
    res = []
    for i in idx[:nsamp]:
        lstar = lev[i]
        Sidx = [i] + [j for j in range(len(A)) if lev[j] > lstar]
        if lev[Sidx].max() - lstar < m:
            continue
        G = Gall[:, Sidx]
        r1 = rank_modP(G); r0 = rank_modP(G[:, 1:])
        char0 = (r1 == r0)
        pint = None
        if char0 and p:
            Gi = evmat(A[Sidx], mon, T)  # small ints; recompute without reduction issues
            pint, mv, _ = member_Zp(Gi[:, 1:], Gi[:, 0], p, 19 if p == 3 else 13)
        res.append((tuple(A[i]), int(lstar), len(Sidx), len(Sidx) - r1, char0, pint))
        print("   k*=", tuple(int(t) for t in A[i]), "l*=", int(lstar), "|S|=", len(Sidx), "dimV(S)>=", len(Sidx) - r1,
              "char0 UM:", char0, "p-integral UM:", pint); sys.stdout.flush()
