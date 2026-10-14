"""Core library for the locality computations (round 13).

Everything lives in the reduced ring R = Q[x_1..x_d]/(x_1^3,...,x_d^3)  (monomials with exponents in {0,1,2}).
x_x^3 lies in J(A) for every A in {-1,0,1}^d and in every local ideal that contains the coordinate width powers,
so all ideals considered contain (x^3) and are compared inside R.

dim J(A)_i (in R) = N_i - H_A(i).

Soundness conventions (mod-p sandwich):
  * every generator is an exact integer form lying in J(A) over Q;
  * dim_p of the ideal generated mod p by the reductions <= dim_Q of the rational ideal  (any subset of rows too);
  * N_{<=i} - rank_p(eval_{<=i}) >= N_{<=i} - rank_Q(eval_{<=i}) = dim_Q J(A)_{<=i}.
  So cumulative equality of the two mod-p numbers proves J_loc = J(A) over Q in all degrees <= i.
  Non-equality is certified only by EXACT ranks over Q (fmpz).
"""
import itertools, math
import numpy as np
from flint import nmod_mat, fmpz_mat

PR = 2147475331          # prime, = 1 mod 2310 (so mod 11, 15, 22 roots of unity exist)


class Ring:
    def __init__(self, d, maxdeg):
        self.d = d; self.maxdeg = maxdeg
        self.mons = []      # mons[t] = int array (N_t, d) of exponents
        self.code2idx = {}
        pw = 3 ** np.arange(d, dtype=np.int64)
        self.pw = pw
        allm = {t: [] for t in range(maxdeg + 2)}
        for e in itertools.product((0, 1, 2), repeat=d):
            t = sum(e)
            if t <= maxdeg + 1:
                allm[t].append(e)
        for t in range(maxdeg + 2):
            arr = np.array(allm[t], dtype=np.int64).reshape(-1, d)
            # sort lexicographically for reproducibility
            self.mons.append(arr)
        self.N = [len(m) for m in self.mons]
        self.idx = []
        for t in range(maxdeg + 2):
            codes = self.mons[t] @ pw
            self.idx.append({int(c): i for i, c in enumerate(codes)})
        # multiplication tables: mul[t][l] = index in degree t+1 of x_l * mon, or -1
        self.mul = []
        for t in range(maxdeg + 1):
            tab = np.full((d, self.N[t]), -1, dtype=np.int64)
            codes = self.mons[t] @ pw
            for l in range(d):
                ok = self.mons[t][:, l] < 2
                nxt = self.idx[t + 1]
                tab[l, ok] = [nxt[int(c + pw[l])] for c in codes[ok]]
            self.mul.append(tab)

    def mul_var(self, vec, t, l, mod=PR):
        """x_l * vec (vec dense over R_t) -> dense over R_{t+1}"""
        out = np.zeros(self.N[t + 1], dtype=object if mod is None else np.int64)
        tab = self.mul[t][l]
        ok = tab >= 0
        if mod is None:
            for a, b in zip(np.nonzero(ok)[0], tab[ok]):
                out[b] += vec[a]
        else:
            np.add.at(out, tab[ok], vec[ok])
            out %= mod
        return out

    def mul_lin(self, vec, t, v, mod=PR):
        """<v,x> * vec"""
        if mod is None:
            out = np.zeros(self.N[t + 1], dtype=object)
            for l in range(self.d):
                if v[l]:
                    tab = self.mul[t][l]; ok = tab >= 0
                    for a, b in zip(np.nonzero(ok)[0], tab[ok]):
                        out[b] += v[l] * vec[a]
            return out
        out = np.zeros(self.N[t + 1], dtype=np.int64)
        for l in range(self.d):
            if v[l] % mod:
                tab = self.mul[t][l]; ok = tab >= 0
                np.add.at(out, tab[ok], (vec[ok] * (v[l] % mod)) % mod)
        return out % mod

    def one(self, mod=PR):
        return np.ones(1, dtype=object if mod is None else np.int64)


def hilbert_modp(A, ring, maxdeg, mod=PR):
    """H_A(i), i=0..maxdeg, via mod-p cumulative rank profile of reduced monomials evaluated on A.
    (cumulative rank_p <= rank_Q, hence N_{<=i} - cum rank_p is an UPPER bound for dim_Q J(A)_{<=i}.)"""
    A = np.array(A, dtype=np.int64)
    rows = []
    ranks = []
    for t in range(maxdeg + 1):
        M = ring.mons[t]
        # value of x^e at k: prod k_x^{e_x}; with k in {-1,0,1}
        vals = np.ones((len(M), len(A)), dtype=np.int64)
        for x in range(ring.d):
            e = M[:, x][:, None]
            kx = A[:, x][None, :]
            vals *= np.where(e == 0, 1, np.where(e == 1, kx, kx * kx))
        rows.append(vals % mod)
        allrows = np.vstack(rows)
        Mt = nmod_mat(allrows.shape[0], allrows.shape[1], [int(z) for z in allrows.ravel()], mod)
        ranks.append(Mt.rank())
    H = [ranks[0]] + [ranks[i] - ranks[i - 1] for i in range(1, len(ranks))]
    return H


def hilbert_exact(A, ring, maxdeg):
    A = np.array(A, dtype=np.int64)
    rows = []; ranks = []
    for t in range(maxdeg + 1):
        M = ring.mons[t]
        vals = np.ones((len(M), len(A)), dtype=np.int64)
        for x in range(ring.d):
            e = M[:, x][:, None]; kx = A[:, x][None, :]
            vals *= np.where(e == 0, 1, np.where(e == 1, kx, kx * kx))
        rows.append(vals)
        allrows = np.vstack(rows)
        Mt = fmpz_mat(allrows.shape[0], allrows.shape[1], [int(z) for z in allrows.ravel()])
        ranks.append(Mt.rank())
    return [ranks[0]] + [ranks[i] - ranks[i - 1] for i in range(1, len(ranks))]


# ---------------- projections ----------------

def bmonos(b, i):
    return [tuple(sum(1 for z in c if z == q) for q in range(b)) for c in itertools.combinations_with_replacement(range(b), i)]


def image(A, F):
    Av = np.array(A, dtype=np.int64)
    Fm = np.array(F, dtype=np.int64).reshape(len(F), -1)
    return sorted(set(map(tuple, (Av @ Fm.T).tolist())))


def jproj_gens(P, b, maxdeg):
    """Exact (over Z) spanning set of NEW generators of J(P) in each degree i<=maxdeg:
    a basis of J(P)_i modulo y*J(P)_{i-1} (selected by exact rank).  Returns {i: [dict gamma->int]}, rho(P)."""
    out = {}
    prev_span = []          # list of integer vectors spanning J(P)_{i-1} (in degree i-1 monomial basis)
    rho = None
    for i in range(1, maxdeg + 1):
        mons = [e for t in range(i + 1) for e in bmonos(b, t)]
        topm = bmonos(b, i); topidx = {e: a for a, e in enumerate(topm)}
        Mt = fmpz_mat(len(P), len(mons), [int(np.prod([p[z] ** e[z] for z in range(b)])) for p in P for e in mons])
        X, nul = Mt.nullspace()
        forms = []
        off = len(mons) - len(topm)
        for c in range(nul):
            v = [int(X[off + r, c]) for r in range(len(topm))]
            if any(v):
                forms.append(v)
        # y * J(P)_{i-1}
        shifted = []
        prevm = bmonos(b, i - 1)
        for v in prev_span:
            for q in range(b):
                w = [0] * len(topm)
                for a, e in enumerate(prevm):
                    if v[a]:
                        e2 = list(e); e2[q] += 1
                        w[topidx[tuple(e2)]] += v[a]
                shifted.append(w)
        base_rank = fmpz_mat(shifted).rank() if shifted else 0
        cur = list(shifted); new = []
        r = base_rank
        for v in forms:
            r2 = fmpz_mat(cur + [v]).rank()
            if r2 > r:
                cur.append(v); new.append(v); r = r2
        out[i] = [{e: v[a] for a, e in enumerate(topm) if v[a]} for v in new]
        # full span of J(P)_i
        prev_span = cur if cur else []
        if r == len(topm) and rho is None:
            rho = i - 1   # J(P)_i = everything  =>  H_P(i) = 0
    return out, rho


def pullback(ring, F, form, i, mod=PR, cache=None):
    """pull back a degree-i form in b variables (dict gamma->int) along y_q = <F_q, x>, reduced in R."""
    if cache is None:
        cache = {}
    b = len(F)

    def prod(gamma):
        key = (tuple(map(tuple, F)), gamma)
        if key in cache:
            return cache[key]
        if sum(gamma) == 0:
            res = ring.one(mod)
        else:
            q = max(z for z in range(b) if gamma[z] > 0)
            g2 = list(gamma); g2[q] -= 1
            res = ring.mul_lin(prod(tuple(g2)), sum(gamma) - 1, F[q], mod)
        cache[key] = res
        return res

    out = np.zeros(ring.N[i], dtype=object if mod is None else np.int64)
    for gamma, c in form.items():
        if mod is None:
            out = out + c * prod(gamma)
        else:
            out = (out + (c % mod) * prod(gamma)) % mod
    return out


# ---------------- ideal dimensions ----------------

def ideal_dims_modp(ring, gens, maxdeg, mod=PR, return_basis=False):
    """gens: {i: [dense vectors over R_i]} (mod p).  Returns dims of the ideal generated, degree 0..maxdeg (mod p)."""
    dims = []; basis = []
    for i in range(maxdeg + 1):
        rows = []
        if i > 0:
            for r in basis:
                for l in range(ring.d):
                    rows.append(ring.mul_var(r, i - 1, l, mod))
        rows += [g % mod for g in gens.get(i, [])]
        if not rows:
            dims.append(0); basis = []; continue
        M = np.array(rows, dtype=np.int64)
        Mt = nmod_mat(M.shape[0], M.shape[1], [int(z) for z in M.ravel()], mod)
        Rr, rk = Mt.rref()
        dims.append(rk)
        ent = np.array([int(x) for x in Rr.entries()], dtype=np.int64).reshape(M.shape[0], M.shape[1])
        basis = list(ent[:rk])
    if return_basis:
        return dims, basis
    return dims


def ideal_dims_exact(ring, gens, maxdeg):
    """exact over Q: gens {i: [object vectors (ints)]}.  Uses fraction-free rref and keeps integer bases."""
    dims = []; basis = []
    for i in range(maxdeg + 1):
        rows = []
        if i > 0:
            for r in basis:
                for l in range(ring.d):
                    rows.append(ring.mul_var(r, i - 1, l, None))
        rows += list(gens.get(i, []))
        if not rows:
            dims.append(0); basis = []; continue
        Mt = fmpz_mat([[int(z) for z in r] for r in rows])
        Rr, den = Mt.rref()
        rk = Mt.rank()
        dims.append(rk)
        basis = []
        for a in range(rk):
            v = np.array([int(Rr[a, c]) for c in range(Mt.ncols())], dtype=object)
            g = 0
            for z in v:
                g = math.gcd(g, int(z))
            if g > 1:
                v = v // g
            basis.append(v)
    return dims
