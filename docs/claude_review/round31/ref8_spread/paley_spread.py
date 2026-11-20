"""Referee's own checks of the VIII.3 inputs and the VIII.4 spread background.

Usage: python3 paley_spread.py q [d_num d_den] [seed] [mode]
  q      prime = 3 mod 4
  d      exponent of P = q^d (default 199/100, as in Astra's computed paragraph)
  mode   'astra' -> S = first (q-1)/2 - 2 positive rows of column 0 (Astra's choice)
         'random' -> S = Bernoulli(rho=0.99) subset of the positive rows (the proof's choice)
Exact integer / Fraction arithmetic throughout; numpy int64 only for reductions mod small moduli
with entries far below 2^63.
"""
import sys, random, math
from fractions import Fraction
from itertools import combinations
import numpy as np


def legendre_matrix(q):
    sq = {x * x % q for x in range(1, q)}
    chi = [0] + [1 if x in sq else -1 for x in range(1, q)]
    K = np.array([[chi[(j - i) % q] for j in range(q)] for i in range(q)], dtype=np.int64)
    return K, chi


def rank_mod_p(M, p):
    """rank of integer matrix M over F_p (p prime), plain Gaussian elimination."""
    A = (np.array(M, dtype=np.int64) % p).copy()
    r = 0
    rows, cols = A.shape
    for c in range(cols):
        piv = None
        for k in range(r, rows):
            if A[k, c] % p:
                piv = k
                break
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), -1, p)
        A[r] = (A[r] * inv) % p
        for k in range(rows):
            if k != r and A[k, c]:
                A[k] = (A[k] - A[k, c] * A[r]) % p
        r += 1
        if r == rows:
            break
    return r


def solve_mod_p(M, R, p):
    """Solve M x = R (mod p) for x, R may have several columns. Returns x (free vars = 0) or None."""
    M = np.array(M, dtype=np.int64) % p
    R = np.array(R, dtype=np.int64) % p
    if R.ndim == 1:
        R = R[:, None]
    A = np.concatenate([M, R], axis=1)
    rows, cols = M.shape
    r = 0
    pivcols = []
    for c in range(cols):
        piv = None
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0:
            continue
        piv = r + nz[0]
        A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), -1, p)
        A[r] = (A[r] * inv) % p
        col = A[:, c].copy()
        col[r] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr):
            A[nzr] = (A[nzr] - np.outer(col[nzr], A[r])) % p
        pivcols.append(c)
        r += 1
        if r == rows:
            break
    if np.any(A[r:, cols:] % p):
        return None
    x = np.zeros((cols, R.shape[1]), dtype=np.int64)
    for k, c in enumerate(pivcols):
        x[c] = A[k, cols:]
    return x


def main():
    q = int(sys.argv[1])
    dn, dd = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (199, 100)
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    mode = sys.argv[5] if len(sys.argv) > 5 else 'astra'
    rng = random.Random(seed)
    assert q % 4 == 3 and all(q % k for k in range(2, int(q ** .5) + 1))
    K, chi = legendre_matrix(q)
    I = list(range(1, q))
    # basic identities
    assert (K.T == -K).all()
    assert (K @ K.T == q * np.eye(q, dtype=np.int64) - 1).all()
    assert (K.sum(0) == 0).all() and (K.sum(1) == 0).all()
    A = [i for i in range(q) if K[i, 0] == 1]
    print(f"q={q}: positive rows of column 0: {len(A)} (=(q-1)/2: {len(A) == (q - 1) // 2})")
    # pair distances after deleting column 0, over all of A
    KI = K[:, 1:]
    dists = set()
    npairs = 0
    for a, b in combinations(A, 2):
        dists.add(int(np.abs(KI[a] - KI[b]).sum()))
        npairs += 1
    print(f"  pairs in A: {npairs}; set of l1 distances on I: {sorted(dists)}")
    if mode == 'astra':
        S = A[:-2]
    else:
        S = [i for i in A if rng.random() < 0.99]
    Sset = set(S)
    Sc = [i for i in range(q) if i not in Sset]
    print(f"  mode={mode}: |S|={len(S)}, |S^c|={len(Sc)}, pairs in S: {len(S)*(len(S)-1)//2}")

    # ---------- moment lifts, proof's construction (random rounding + completion mod q) ----------
    J = q // 4
    KS = K[S]
    # completion system: K_{S^c}^T x = rhs (mod q); matrix indexed (j, i in S^c)
    Msys = K[Sc].T  # shape q x |Sc|
    rhs_all = []
    taus = []
    ds = []
    for s in range(1, J + 1):
        c = [pow(j, s, q) for j in range(q)]
        c = [x - q if x > q // 2 else x for x in c]
        # mean-zero random rounding between adjacent representatives
        d = []
        for x in c:
            if x == 0:
                d.append(0)
            else:
                alt = x - q if x > 0 else x + q
                d.append(alt if rng.random() < abs(x) / q else x)
        d = np.array(d, dtype=np.int64)
        Kd = K @ d
        assert (Kd % q == 0).all(), "tau not integral"
        tau = Kd // q
        taus.append(tau)
        ds.append(d)
        rhs_all.append((-(K[S].T @ tau[S])) % q)
    X = solve_mod_p(Msys, np.array(rhs_all).T, q)
    assert X is not None, "completion mod q not solvable"
    maxsup = 0
    maxtau = 0
    for s in range(1, J + 1):
        tau = taus[s - 1]
        d = ds[s - 1]
        maxtau = max(maxtau, int(np.abs(tau).max()))
        tp = np.zeros(q, dtype=np.int64)
        tp[S] = tau[S]
        xs = X[:, s - 1] % q
        # random zero-mean rounding of residues on S^c, magnitude < q
        for k, i in enumerate(Sc):
            x = int(xs[k])
            tp[i] = x - q if rng.random() < x / q else x
        assert ((K.T @ tp) % q == 0).all()
        lift = d - K.T @ tp
        assert ((lift - np.array([pow(j, s, q) for j in range(q)])) % q == 0).all()
        ev = KS[:, 1:] @ lift[1:]
        assert len(set(ev.tolist())) == 1, "row evaluations not constant"
        ev0 = KS @ lift
        assert len(set(ev0.tolist())) == 1
        # randomisation q K_{S^c}^T xi keeps congruence and constancy
        xi = np.array([int(round(rng.gauss(0, q ** 0.12))) for _ in Sc], dtype=np.int64)
        D = lift + q * (K[Sc].T @ xi)
        assert ((D - lift) % q == 0).all()
        assert len(set((KS[:, 1:] @ D[1:]).tolist())) == 1
        maxsup = max(maxsup, int(np.abs(lift).max()))
    print(f"  moment lifts s=1..{J}: congruences and constant K_S-evaluations (with and without column 0) OK;"
          f" max|tau|={maxtau} (sqrt(q log q)={math.sqrt(q*math.log(q)):.1f}); max sup of lift={maxsup}"
          f" (q^1.5 sqrt(log q)={q**1.5*math.sqrt(math.log(q)):.0f})")

    # ---------- spread background ----------
    # b_j = 1/12 + (j/q + z_j)/3 = (q + 4j + 4 q z_j)/(12 q); T = labels with some z placing b_j within 10 P^{-1/2} of Z/2.
    # dist(b, Z/2) = dist(n, 6qZ)/(12q); compare r/(12q) <= 10 q^{-d/2}  <=>  r^{2dd} q^{dn} <= (120 q)^{2dd}
    def near(n):
        r = n % (6 * q)
        r = min(r, 6 * q - r)
        return r ** (2 * dd) * q ** dn <= (120 * q) ** (2 * dd)
    T = [j for j in range(q) if any(near(q + 4 * j + 4 * q * z) for z in range(3))]
    rkT = rank_mod_p(K[np.ix_(Sc, T)], 3)
    print(f"  spread: |T|={len(T)} (sqrt q={q**.5:.1f}); 0 in T: {0 in T}; rank_F3 K_(S^c,T) = {rkT}")
    if rkT < len(T):
        print("  FULL COLUMN RANK FAILS")
        return
    J1 = np.arange(q, dtype=np.int64)
    KJ = K @ J1
    assert (KJ % q == 0).all()
    tau1 = KJ // q
    u = np.zeros(q, dtype=np.int64)
    u[S] = tau1[S] % 3
    # (K^T u)_T = 0 mod 3:  K^T_{T,S^c} u_{S^c} = -K^T_{T,S} u_S
    Mt = K.T[np.ix_(T, Sc)]
    rhs = (-(K.T[np.ix_(T, S)] @ u[S])) % 3
    x = solve_mod_p(Mt, rhs, 3)
    assert x is not None
    u[Sc] = x[:, 0] % 3
    qinv3 = pow(q, -1, 3)
    z = (-qinv3 * (K.T @ u)) % 3
    assert all(z[j] == 0 for j in T)
    w = (tau1 + K @ z) % 3
    assert len(set(w[S].tolist())) == 1, "tau1 + K z not constant mod 3 on S"
    b = [Fraction(q + 4 * j + 4 * q * int(z[j]), 12 * q) for j in range(q)]
    # exact windings
    for cols in (range(q), I):
        vals = [2 * sum(int(K[i, j]) * b[j] for j in cols) for i in S]
        assert all((v - vals[0]).denominator == 1 for v in vals), "windings not integral"
    vals = [2 * sum(int(K[i, j]) * b[j] for j in I) for i in S]
    kwind = [int(v - vals[0]) for v in vals]
    print(f"  integral row windings on S (with and without column 0): OK; winding range [{min(kwind)},{max(kwind)}]")
    # clearance
    for j in range(q):
        n = q + 4 * j + 4 * q * int(z[j])
        if j in T:
            assert Fraction(1, 12) <= b[j] <= Fraction(5, 12)
        else:
            assert not near(n)
    # separation over I, in quarter-turn units mod 1, and up to conjugation
    def dist1(f):
        f = f - math.floor(f)
        return min(f, 1 - f)
    mind = min(dist1(b[j] - b[k]) for j, k in combinations(I, 2))
    minc = min(dist1(b[j] + b[k]) for j in I for k in I if k >= j)
    print(f"  separation: min |b_j-b_k| = {mind} (1/(3q)={Fraction(1,3*q)}), min |b_j+b_k| = {minc} (1/(6q)={Fraction(1,6*q)})")
    assert mind >= Fraction(1, 3 * q) and minc >= Fraction(1, 6 * q)
    # distance of every single b_j (j in I) to Z/2
    dmin1 = min(min(dist1(2 * b[j]) / 2, 1) for j in I)
    print(f"  min dist(b_j, Z/2) over I = {float(dmin1):.3e}; 10 P^(-1/2) = {10*q**(-dn/(2*dd)):.3e}")
    # short-word centre count: balanced words of length 4 at a centre (v.b in Z/2) -- first power sum must vanish
    cnt = 0; bad = 0
    lab = I[:60]
    for a, bb in combinations(lab, 2):
        for c_, e_ in combinations(lab, 2):
            if len({a, bb, c_, e_}) < 4:
                continue
            val = b[a] + b[bb] - b[c_] - b[e_]
            if (2 * val).denominator == 1:
                cnt += 1
                if (a + bb - c_ - e_) % q:
                    bad += 1
    print(f"  length-4 balanced words (distinct labels <60) at a centre: {cnt}; with nonvanishing first power sum: {bad}")


if __name__ == '__main__':
    main()
