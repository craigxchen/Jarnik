"""Isotypic (Specht-module) computation of the regularity / maximal vanishing order of measures
on a finite S_M-stable point set A in Z^M with constant coordinate sum s (a shell S^M(p,n) or a
slice Q^M(p,n)).

Notation:  reg(A) = min{d : polynomials of degree <= d interpolate A}
                  = max{ord(c) : c != 0 a measure on A}       (duality, see compute.md, Lemma 1)

For a partition lam of M, with t the row-reading tableau, a_t = row symmetriser, b_t = signed
column antisymmetriser:
  * the lam-isotypic measures are detected inside E_lam = a_t b_t C^A, of dimension
    N_lam = sum_orbits K_{lam, mu(orbit)}; we use the spanning vectors c_T = a_t b_t 1_{k_T},
    T a semistandard tableau of shape lam whose content is the orbit's value multiset
    (k_T = the point carrying T's entries at the positions of t);
  * <c_T, P> = (b_t a_t P)(k_T), and b_t a_t Poly = Sym . span{F_t^S : S standard}
    (Ariki-Terasoma-Yamada higher Specht polynomials F_t^S = b_t a_t x_t^{i(S)});
  * on A the power sum p_1 = s is constant, so Sym acts through C[p_2..p_M].
Rows (generators) are g(v_O) F_t^S(k_T) with g a monomial in p_2..p_M; the generator has degree
wdeg(g) + charge(S).  h_lam(d) = rank of the generators of degree <= d (mod p) and
reg_lam(A) = min{d : h_lam(d) = N_lam} = max ord over nonzero lam-isotypic measures.

Rigour: rank mod p <= rank over Q, and the generators are genuine pairings with polynomials of
degree <= d, so h_lam(d) mod p is a lower bound for the true h_lam(d); hence the computed
reg_lam is a rigorous UPPER bound for the true reg_lam (whatever the completeness of the
generators).  Equality is supported by (i) two primes, (ii) agreement of
sum_lam f^lam h_lam(d) with the direct Hilbert function where that is computable, and
(iii) exact certificates (cert.py) where stated.
"""
import itertools
import time
import numpy as np
from math import factorial
from modla import Echelon, mm, P1, P2
from combi import (partitions, conjugate, hook_dim, row_positions, col_positions,
                   standard_tableaux, aty_index, ssyt_fillings, content_of, pos_of)


def group_perms(blocks, M, signed):
    """all permutations preserving each block (list of position lists); returns (perms, signs)
    perms[i] is an index array with (sigma k)[pos] = k[perms[i][pos]]."""
    per_block = []
    for blk in blocks:
        opts = []
        for pr in itertools.permutations(range(len(blk))):
            # sign of pr
            sg = 1
            seen = [False] * len(pr)
            for a in range(len(pr)):
                if not seen[a]:
                    l = 0
                    b = a
                    while not seen[b]:
                        seen[b] = True
                        b = pr[b]
                        l += 1
                    if l % 2 == 0:
                        sg = -sg
            opts.append(([blk[x] for x in pr], sg))
        per_block.append(opts)
    perms, signs = [], []
    for combo in itertools.product(*per_block):
        perm = list(range(M))
        sg = 1
        for blk, (img, s) in zip(blocks, combo):
            for a, b in zip(blk, img):
                perm[a] = b
            sg *= s
        perms.append(perm)
        signs.append(sg if signed else 1)
    return np.array(perms, dtype=np.int64), np.array(signs, dtype=np.int64)


def _powers(Y, emax, p):
    """Y[..., q] -> list pw[e] (e = 0..emax) of arrays Y^e mod p."""
    pw = [np.ones_like(Y)]
    for e in range(1, emax + 1):
        pw.append((pw[-1] * Y) % p)
    return pw


def subset_dp(pw, exps, L, p, signed):
    """permanent (signed=False) or determinant (signed=True) of the L x L matrix
    A[i][q] = y_q^{exps[i]}, vectorised; pw[e][..., q] = y_q^e."""
    shape = pw[0].shape[:-1]
    dp = {0: np.ones(shape, dtype=np.int64)}
    for mask in range(1, 1 << L):
        i = bin(mask).count('1') - 1
        e = exps[i]
        acc = np.zeros(shape, dtype=np.int64)
        for q in range(L):
            if mask >> q & 1:
                term = (dp[mask ^ (1 << q)] * pw[e][..., q]) % p
                if signed and (bin(mask >> (q + 1)).count('1') % 2):
                    acc = (acc - term) % p
                else:
                    acc = (acc + term) % p
        dp[mask] = acc
    return dp[(1 << L) - 1]


class LamData:
    """per-partition static data: tableaux, indices, group elements."""

    def __init__(self, lam):
        self.lam = lam
        self.M = M = sum(lam)
        self.rows = row_positions(lam)
        self.cols = col_positions(lam)
        P = pos_of(lam)
        self.SYT = standard_tableaux(lam)
        self.alpha = []      # exponent vector by position
        self.charge = []
        for S in self.SYT:
            I, ch = aty_index(lam, S)
            a = [0] * M
            for b, v in I.items():
                a[P[b]] = v
            self.alpha.append(a)
            self.charge.append(ch)
        self.f = len(self.SYT)
        nR = 1
        for r in self.rows:
            nR *= factorial(len(r))
        nC = 1
        for c in self.cols:
            nC *= factorial(len(c))
        costC = nC * sum(len(r) * 2 ** len(r) for r in self.rows)
        costR = nR * sum(len(c) * 2 ** len(c) for c in self.cols)
        self.method = 'C' if costC <= costR else 'R'
        self.nR, self.nC = nR, nC
        if self.method == 'C':
            self.perms, self.signs = group_perms(self.cols, M, True)
        else:
            self.perms, _ = group_perms(self.rows, M, False)
            # inverse permutations: beta = alpha[tau_inv]
            inv = np.empty_like(self.perms)
            for i, pr in enumerate(self.perms):
                inv[i, pr] = np.arange(M)
            self.perms_inv = inv

    def F_values(self, K, p, chunk_elems=4_000_000):
        """F_t^S(k) mod p for all S (rows) and all points k in K (n x M int array)."""
        n = K.shape[0]
        out = np.zeros((self.f, n), dtype=np.int64)
        if n == 0:
            return out
        Kp = K % p
        emax = max(max(a) for a in self.alpha)
        if self.method == 'C':
            nC = len(self.perms)
            step = max(1, chunk_elems // (nC * self.M))
            # row patterns
            pats = {}
            for si, a in enumerate(self.alpha):
                for r, rp in enumerate(self.rows):
                    pats.setdefault((r, tuple(sorted(a[x] for x in rp))), None)
            for s0 in range(0, n, step):
                Kc = Kp[s0:s0 + step]
                KS = Kc[:, self.perms]                     # (m, nC, M)
                vals = {}
                for (r, e) in pats:
                    rp = self.rows[r]
                    Y = KS[:, :, rp]
                    pw = _powers(Y, max(e) if e else 0, p)
                    vals[(r, e)] = subset_dp(pw, e, len(rp), p, False)
                for si, a in enumerate(self.alpha):
                    prod = None
                    for r, rp in enumerate(self.rows):
                        v = vals[(r, tuple(sorted(a[x] for x in rp)))]
                        prod = v if prod is None else (prod * v) % p
                    tot = (prod * self.signs[None, :]).sum(axis=1) % p
                    out[si, s0:s0 + step] = tot
        else:
            step = max(1, chunk_elems // self.M)
            for s0 in range(0, n, step):
                Kc = Kp[s0:s0 + step]
                colpw = []
                for c, cp in enumerate(self.cols):
                    colpw.append(_powers(Kc[:, cp], emax, p))
                cache = {}
                for si, a in enumerate(self.alpha):
                    a = np.array(a, dtype=np.int64)
                    betas = {}
                    for inv in self.perms_inv:
                        b = tuple(a[inv])
                        betas[b] = betas.get(b, 0) + 1
                    tot = np.zeros(Kc.shape[0], dtype=np.int64)
                    for b, mult in betas.items():
                        prod = None
                        for c, cp in enumerate(self.cols):
                            key = (c, tuple(b[x] for x in cp))
                            if key not in cache:
                                cache[key] = subset_dp(colpw[c], key[1], len(cp), p, True)
                            v = cache[key]
                            prod = v if prod is None else (prod * v) % p
                        tot = (tot + (prod * (mult % p)) % p) % p
                    out[si, s0:s0 + Kc.shape[0]] = tot
        return out


def power_sums(V, M, p):
    """V: (norb, M) orbit representatives -> dict i -> p_i values mod p (i = 2..M)."""
    Vp = V % p
    out = {}
    pw = Vp.copy()
    for i in range(1, M + 1):
        if i > 1:
            pw = (pw * Vp) % p
        out[i] = pw.sum(axis=1) % p
    return out


def sym_filtration(V, M, p, emax=10 ** 9, parity=None):
    """greedy basis of Sym'_{<=e} = C[p_2..p_M]_{<=e} restricted to the orbit set V
    (if parity is 0/1: only monomials of weighted degree = parity mod 2).
    Returns list of (e_j, values_j, key) in birth order, and the profile u(e)."""
    norb = V.shape[0]
    ps = power_sums(V, M, p)
    E = Echelon(norb, p)
    basis = []
    prof = []
    # monomials of weighted degree e: generated as partitions of e into parts 2..M
    level = {0: {(): np.ones(norb, dtype=np.int64)}}
    e = 0
    while E.rank < norb and e <= emax:
        if e > 0:
            cur = {}
            for i in range(2, M + 1):
                if e - i < 0:
                    continue
                for key, val in level[e - i].items():
                    if key and key[-1] > i:
                        continue   # parts weakly increasing -> unique representation
                    cur[key + (i,)] = (val * ps[i]) % p
            level[e] = cur
        keys = list(level[e].keys())
        if keys and (parity is None or e % 2 == parity):
            B = np.array([level[e][k] for k in keys], dtype=np.int64)
            ind = E.add(B)
            for i in ind:
                basis.append((e, B[i].copy(), keys[i]))
        prof.append(E.rank)
        # free memory of old levels no longer needed
        for old in list(level.keys()):
            if old < e - M:
                del level[old]
        e += 1
    return basis, prof


def columns_for(lam, orbits):
    """columns (k_T, orbit index) for the semistandard tableaux of all orbits."""
    P = pos_of(lam)
    M = sum(lam)
    pts, oidx, relevant = [], [], []
    for o, v in enumerate(orbits):
        fills = ssyt_fillings(lam, content_of(v))
        if fills:
            relevant.append(o)
        for T in fills:
            k = [0] * M
            for b, val in T.items():
                k[P[b]] = val
            pts.append(k)
            oidx.append(o)
    return np.array(pts, dtype=np.int64).reshape(-1, M), np.array(oidx, dtype=np.int64), relevant


def lam_profile(lam, orbits, p, ld=None, verbose=False, max_rows=6000, dstop=None, dguess=None):
    """rank profile h_lam(d) (mod p) and reg_lam for the point set with the given orbits.
    Returns dict with N, reg, profile, timings.  If dstop is given, stop after degree dstop
    (then reg may be None = 'greater than dstop')."""
    t0 = time.time()
    M = sum(lam)
    if ld is None:
        ld = LamData(lam)
    K, oidx, relevant = columns_for(lam, orbits)
    N = K.shape[0]
    res = {'lam': lam, 'N': N, 'f': ld.f, 'method': ld.method}
    if N == 0:
        res.update(reg=None, profile=[], time=0.0)
        return res
    V = np.array([orbits[o] for o in relevant], dtype=np.int64)
    remap = {o: i for i, o in enumerate(relevant)}
    oloc = np.array([remap[o] for o in oidx], dtype=np.int64)
    basis, sprof = sym_filtration(V, M, p)
    t1 = time.time()
    Fv = ld.F_values(K, p)
    t2 = time.time()
    G = np.array([b[1] for b in basis], dtype=np.int64)          # (nbasis, nrel)
    ge = np.array([b[0] for b in basis], dtype=np.int64)
    ch = np.array(ld.charge, dtype=np.int64)
    dmax_gen = int(ge.max() + ch.max())
    if dstop is not None:
        dmax_gen = min(dmax_gen, dstop)
    D0 = dmax_gen if dguess is None else min(dmax_gen, dguess)
    E = Echelon(N, p)
    indep_deg = []
    ngen = 0

    def run_degrees(dlo, dhi):
        nonlocal ngen
        pairs = []
        for d in range(dlo, dhi + 1):
            for s in range(ld.f):
                for j in np.flatnonzero(ge == d - ch[s]):
                    pairs.append((int(j), s, d))
        mr = max(500, min(max_rows, 20_000_000 // max(1, N)))
        for b0 in range(0, len(pairs), mr):
            if E.rank == N:
                break
            sub = pairs[b0:b0 + mr]
            js = np.array([x[0] for x in sub], dtype=np.int64)
            ss = np.array([x[1] for x in sub], dtype=np.int64)
            B = (G[js][:, oloc] * Fv[ss]) % p
            ind = E.add(B)
            ngen += len(sub)
            indep_deg.extend(sub[i][2] for i in ind)

    run_degrees(0, D0)
    d = D0
    while E.rank < N and d < dmax_gen:
        d += 1
        run_degrees(d, d)
    if E.rank < N:
        if dstop is not None:
            prof = [sum(1 for x in indep_deg if x <= dd) for dd in range(dmax_gen + 1)]
            res.update(reg=None, profile=prof, time=time.time() - t0, partial=True)
            return res
        raise RuntimeError(f'columns dependent? lam={lam} rank {E.rank} < N {N}')
    reg = max(indep_deg)
    profile = [sum(1 for x in indep_deg if x <= dd) for dd in range(reg + 1)]
    res['ngen'] = ngen
    res['beyond_guess'] = (dguess is not None and reg > dguess)
    res.update(reg=len(profile) - 1, profile=profile, time=time.time() - t0,
               t_sym=t1 - t0, t_F=t2 - t1, sym_reg=len(sprof) - 1, nrel=len(relevant))
    return res


def neg_orbit(v):
    return tuple(sorted((-x for x in v), reverse=True))


def lam_profile_parity(lam, orbits, p, ld=None, max_rows=6000, dguess=None, extra=6):
    """same output as lam_profile, for an orbit set closed under k -> -k (s = 0), using the
    decomposition E_lam = E^+ (+) E^- under negation.  Columns: one orbit of each pair {O,-O}
    and all self-negating orbits; generators of degree parity eps act on E^eps only."""
    t0 = time.time()
    M = sum(lam)
    if ld is None:
        ld = LamData(lam)
    oset = set(orbits)
    assert all(neg_orbit(v) in oset for v in orbits), 'orbit set not closed under negation'
    K, oidx, relevant = columns_for(lam, orbits)
    N = K.shape[0]
    res = {'lam': lam, 'N': N, 'f': ld.f, 'method': ld.method, 'parity': True}
    if N == 0:
        res.update(reg=None, profile=[], time=0.0)
        return res
    kept_orb, selfneg = [], set()
    for o in relevant:
        v = orbits[o]
        w = neg_orbit(v)
        if w == v:
            kept_orb.append(o)
            selfneg.add(o)
        elif v > w:
            kept_orb.append(o)
    kset = set(kept_orb)
    kcols = np.array([i for i in range(N) if int(oidx[i]) in kset], dtype=np.int64)
    Kk = K[kcols]
    ok_idx = oidx[kcols]
    Fv = ld.F_values(Kk, p)
    ch = np.array(ld.charge, dtype=np.int64)
    # dimensions of E^eps
    dim = {0: 0, 1: 0}            # key = parity of the generator degree (0: E^+, 1: E^-)
    for o in kept_orb:
        cols = np.flatnonzero(ok_idx == o)
        KO = len(cols)
        if o not in selfneg:
            dim[0] += KO
            dim[1] += KO
            continue
        rk = {}
        for eps in (0, 1):
            rows = Fv[np.flatnonzero(ch % 2 == eps)][:, cols]
            Ee = Echelon(KO, p)
            if rows.shape[0]:
                Ee.add(rows)
            rk[eps] = Ee.rank
        Ea = Echelon(KO, p)
        Ea.add(Fv[:, cols])
        if not (Ea.rank == KO and rk[0] + rk[1] == KO):
            raise RuntimeError(f'parity split inconsistent at orbit {orbits[o]} lam={lam}')
        dim[0] += rk[0]
        dim[1] += rk[1]
    V = np.array([orbits[o] for o in kept_orb], dtype=np.int64)
    remap = {o: i for i, o in enumerate(kept_orb)}
    oloc = np.array([remap[int(o)] for o in ok_idx], dtype=np.int64)
    D0 = dguess if dguess is not None else 4 * M * M
    emax = D0 + extra
    filt = {}
    for par in (0, 1):
        basis, _ = sym_filtration(V, M, p, emax=emax, parity=par)
        filt[par] = (np.array([b[1] for b in basis], dtype=np.int64).reshape(len(basis), V.shape[0]),
                     np.array([b[0] for b in basis], dtype=np.int64))
    out = {}
    ngen = 0
    for eps in (0, 1):
        target = dim[eps]
        if target == 0:
            out[eps] = []
            continue
        E = Echelon(len(kcols), p)
        indep_deg = []
        dmax = emax
        def run(dlo, dhi):
            nonlocal ngen
            pairs = []
            for d in range(dlo, dhi + 1):
                if d % 2 != eps:
                    continue
                for s in range(ld.f):
                    par = (d - int(ch[s])) % 2
                    G, ge = filt[par]
                    for j in np.flatnonzero(ge == d - ch[s]):
                        pairs.append((int(j), s, d, par))
            mr = max(500, min(max_rows, 20_000_000 // max(1, N)))
            for b0 in range(0, len(pairs), mr):
                if E.rank == target:
                    break
                sub = pairs[b0:b0 + mr]
                B = np.empty((len(sub), len(kcols)), dtype=np.int64)
                for par in (0, 1):
                    idx = [i for i, x in enumerate(sub) if x[3] == par]
                    if not idx:
                        continue
                    js = np.array([sub[i][0] for i in idx], dtype=np.int64)
                    ss = np.array([sub[i][1] for i in idx], dtype=np.int64)
                    B[idx] = (filt[par][0][js][:, oloc] * Fv[ss]) % p
                ind = E.add(B)
                ngen += len(sub)
                indep_deg.extend(sub[i][2] for i in ind)
        run(0, D0)
        d = D0
        while E.rank < target and d < dmax:
            d += 1
            run(d, d)
        if E.rank < target:
            raise RuntimeError(f'parity {eps}: rank {E.rank} < target {target} by degree {dmax} lam={lam}')
        if E.rank > target:
            raise RuntimeError('rank exceeds target')
        out[eps] = indep_deg
    allind = out[0] + out[1]
    reg = max(allind)
    profile = [sum(1 for x in allind if x <= dd) for dd in range(reg + 1)]
    if profile[-1] != N:
        raise RuntimeError('parity profile does not reach N')
    res.update(reg=reg, profile=profile, time=time.time() - t0, ngen=ngen,
               reg_even=max(out[0]) if out[0] else None, reg_odd=max(out[1]) if out[1] else None,
               dim_even=dim[0], dim_odd=dim[1], beyond_guess=(dguess is not None and reg > dguess))
    return res
