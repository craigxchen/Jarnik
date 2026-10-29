"""Round 28 theory lane: own library (written fresh; conventions = round-27 killer's klib.build_box).
Slab box: k' in {-1,0,1}^d, k_n = -sum k', slabs (S, L): |sum_{x in S} k_x| <= L (S subset {0..n-1}).
On subsets of the ternary cube k_x^3 = k_x, so everything is computed in STANDARD monomials (exponents <= 2):
  M_e(S) = { (sum_k c_k k^alpha)_{alpha std, |alpha|=e} : c in V(S,e) },  J(S)_e (mod J(C)) = M_e(S)^perp (dot product).
Mod-p routines use flint.nmod_mat; exact routines use fmpz_mat / fmpq_mat.
Central symmetry A = -A of every slab box gives M_e(A_{-l}) = M_e(A_l) (sign (-1)^e), so L_pair is halved."""
import itertools, math
import numpy as np
import flint

P1 = 2147483629
P2 = 1000000007

# ---------------- boxes ----------------
def canon(n, S):
    Sc = tuple(x for x in range(n) if x not in S)
    return min(tuple(S), Sc, key=lambda s: (len(s), s))

def syl(k):
    H = [[1]]
    for _ in range(k):
        H = [r + r for r in H] + [r + [-x for x in r] for r in H]
    n = 2 ** k
    return n, tuple(sorted(set(canon(n, tuple(x for x in range(n) if r[x] == 1)) for r in H[1:])))

def legendre(a, q):
    a %= q
    return 0 if a == 0 else (1 if pow(a, (q - 1) // 2, q) == 1 else -1)

def sdc(q):
    idx = ['inf'] + list(range(q)); M = q + 1
    K = np.zeros((M, M), dtype=int)
    for a, x in enumerate(idx):
        for b, y in enumerate(idx):
            if x == 'inf': K[a, b] = 1
            elif y == 'inf': K[a, b] = -1
            else: K[a, b] = (1 if x == y else 0) + legendre(y - x, q)
    if not ((K + K.T) == 2 * np.eye(M, dtype=int)).all(): K = K.T
    Dg = np.diag(K[:, 0]); K = Dg @ K @ Dg
    pats = set()
    for a in range(1, M):
        S = [x for x in range(M) if K[x, a] == -1]
        pats.add(canon(M, tuple(S))); pats.add(canon(M, tuple(sorted(set(S) ^ {0, a}))))
    return M, tuple(sorted(pats))

def kfull(kp):
    return tuple(kp) + (-sum(kp),)

def build_box(n, slabs):
    d = n - 1
    out = []
    for kp in itertools.product((-1, 0, 1), repeat=d):
        kf = kp + (-sum(kp),)
        if all(abs(sum(kf[x] for x in S)) <= L for S, L in slabs):
            out.append(kp)
    return out

def Lpair_half(d):
    L = []
    for x in range(d):
        e = [0] * d; e[x] = 1; L.append(tuple(e))
    for x in range(d):
        for y in range(x + 1, d):
            e = [0] * d; e[x] = 1; e[y] = -1; L.append(tuple(e))
    return L

def Lcube_half(d):
    """one representative of each +-pair of nonzero {-1,0,1}^d shifts"""
    out = []
    for l in itertools.product((-1, 0, 1), repeat=d):
        if any(l) and next(v for v in l if v) > 0: out.append(l)
    return out

def shiftset(S, l, Sset=None):
    Sset = Sset if Sset is not None else set(S)
    return [k for k in S if tuple(a + b for a, b in zip(k, l)) in Sset]

# ---------------- standard monomials ----------------
_sm = {}
def smonos(d, e):
    key = (d, e)
    if key not in _sm:
        _sm[key] = [a for a in itertools.product(range(3), repeat=d) if sum(a) == e]
    return _sm[key]

def smonos_lt(d, e):
    out = []
    for t in range(e): out += smonos(d, t)
    return out

def T(d, e):
    return len(smonos(d, e))

def evalmat(S, mons):
    """int64 matrix |mons| x |S| of values k^alpha (entries in {-1,0,1})"""
    P = np.array(S, dtype=np.int64)
    if len(mons) == 0: return np.zeros((0, len(S)), dtype=np.int64)
    Ex = np.array(mons, dtype=np.int64)
    out = np.ones((len(mons), len(S)), dtype=np.int64)
    for x in range(P.shape[1]):
        col = P[:, x][None, :] ** Ex[:, x][:, None]
        out = out * col
    return out

def nmod(A, p):
    return flint.nmod_mat(A.shape[0], A.shape[1], [int(v) % p for v in A.flatten()], p)

def Me_modp(S, e, p, d):
    """nmod_mat (rows = std monomials of degree e) whose columns span M_e(S) mod p; None if 0."""
    if len(S) == 0: return None
    lt = smonos_lt(d, e)
    if lt:
        Mlt = nmod(evalmat(S, lt), p)
        K, nul = Mlt.nullspace()
        if nul == 0: return None
        K = flint.nmod_mat(len(S), nul, [int(K[i, j]) for i in range(len(S)) for j in range(nul)], p)
    else:
        K = flint.nmod_mat(len(S), len(S), [1 if i == j else 0 for i in range(len(S)) for j in range(len(S))], p)
    E = nmod(evalmat(S, smonos(d, e)), p)
    return E * K

def colbasis(Mm, p):
    """column basis (as nmod_mat) of the column space of Mm"""
    if Mm is None: return None
    nr = Mm.nrows()
    R, r = Mm.transpose().rref()
    if r == 0: return None
    return flint.nmod_mat(nr, r, [int(R[j, i]) for i in range(nr) for j in range(r)], p)

def hcat(mats, p):
    mats = [m for m in mats if m is not None]
    if not mats: return None
    nr = mats[0].nrows()
    cols = sum(m.ncols() for m in mats)
    ent = []
    for i in range(nr):
        for m in mats:
            ent += [int(m[i, j]) for j in range(m.ncols())]
    return flint.nmod_mat(nr, cols, ent, p)

def H_modp(S, e, p, d):
    M = Me_modp(S, e, p, d)
    return 0 if M is None else M.rank()

def shift_sum_modp(S, Ls, e, p, d, stop_at=None):
    """column basis of sum_{l in Ls} M_e(S_l) (mod p)"""
    Sset = set(S); acc = None; cur = 0
    for l in Ls:
        Sl = shiftset(S, l, Sset)
        if len(Sl) <= 0: continue
        Ml = Me_modp(Sl, e, p, d)
        if Ml is None: continue
        acc = colbasis(hcat([acc, Ml], p), p)
        cur = 0 if acc is None else acc.ncols()
        if stop_at is not None and cur >= stop_at: break
    return acc

def shift_dims(S, Ls, e, p, d):
    hs = H_modp(S, e, p, d)
    acc = shift_sum_modp(S, Ls, e, p, d, stop_at=hs)
    return hs, (0 if acc is None else acc.ncols())

# ---------------- exact membership in J(S)_e ----------------
def J_member_space(S, e, d, forms):
    """forms: list of dicts {std alpha (deg e): int coeff}. Returns an fmpq basis (list of lists) of the space of
    coefficient vectors a with sum_i a_i F_i in J(S)_e, i.e. (sum a_i F_i)|_S in span(std monomials of deg < e)|_S. Exact."""
    lt = smonos_lt(d, e)
    E = evalmat(S, lt).T if lt else np.zeros((len(S), 0), dtype=np.int64)     # |S| x |lt|
    Fv = np.zeros((len(S), len(forms)), dtype=object)
    P = np.array(S, dtype=np.int64)
    for i, F in enumerate(forms):
        col = np.zeros(len(S), dtype=object)
        for a, c in F.items():
            v = np.ones(len(S), dtype=np.int64)
            for x, ax in enumerate(a):
                if ax: v = v * P[:, x] ** ax
            col = col + c * v.astype(object)
        Fv[:, i] = col
    Mfull = flint.fmpz_mat([[int(E[r, j]) for j in range(E.shape[1])] + [int(-Fv[r, i]) for i in range(len(forms))] for r in range(len(S))])
    X, nul = Mfull.nullspace()
    nE = E.shape[1]
    vecs = [[int(X[nE + i, j]) for i in range(len(forms))] for j in range(nul)]
    if not vecs: return []
    Mv = flint.fmpq_mat(vecs)
    R, r = Mv.rref()
    return [[R[i, j] for j in range(len(forms))] for i in range(r)]

# ---------------- power-of-linear-form obstructions ----------------
def power_obstructions(S, e, d, Ls, p=P1, vals=(-1, 0, 1)):
    """phi in vals^d (up to sign) with phi^e in cap_{l} J(S_l)_e but not in J(S)_e (mod p).
    Returns (H_S(e), dim Sigma, list of phi, dim span of their classes mod J(S)_e)."""
    from math import factorial
    MA = Me_modp(S, e, p, d)
    HA = 0 if MA is None else MA.rank()
    Sg = shift_sum_modp(S, Ls, e, p, d)
    dS = 0 if Sg is None else Sg.ncols()
    if HA == 0: return HA, dS, [], 0
    mons = smonos(d, e)
    mult = [factorial(e) // math.prod(factorial(int(a)) for a in al) for al in mons]
    phis = [ph for ph in itertools.product(vals, repeat=d) if any(ph) and next(v for v in ph if v) > 0]
    rows = []
    for ph in phis:
        for i, al in enumerate(mons):
            v = mult[i]
            for x in range(d):
                if al[x]: v *= ph[x] ** al[x]
            rows.append(v % p)
    PH = flint.nmod_mat(len(phis), len(mons), rows, p)
    PS = None if Sg is None else PH * Sg
    PM = PH * MA
    obs = []
    for i, ph in enumerate(phis):
        if PS is not None and any(int(PS[i, j]) for j in range(PS.ncols())): continue
        if any(int(PM[i, j]) for j in range(PM.ncols())): obs.append(i)
    dim = 0
    if obs:
        dim = flint.nmod_mat(len(obs), PM.ncols(), [int(PM[i, j]) for i in obs for j in range(PM.ncols())], p).rank()
    return HA, dS, [phis[i] for i in obs], dim
