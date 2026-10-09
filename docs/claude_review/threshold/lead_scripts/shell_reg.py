# For S = {k in Z^M : sum k = s, ||k||_1 = t}, compute reg(S) = max m such that some nonzero c on S
# is orthogonal to all polynomials of degree < m (in the M-1 coordinates k_2..k_M).
# Rank computed modulo a large prime (rank mod p <= rank over Q, so reg_mod >= reg_Q is an upper bound
# check; we also re-check with a second prime).
import itertools, sys
P1 = 2**61 - 1
def shell(M, s, t):
    pts = []
    def rec(prefix, remaining_l1):
        j = len(prefix)
        if j == M - 1:
            last = s - sum(prefix)
            if abs(last) == remaining_l1:
                pts.append(tuple(prefix) + (last,))
            return
        for a in range(-remaining_l1, remaining_l1 + 1):
            rec(prefix + [a], remaining_l1 - abs(a))
    rec([], t)
    return pts
def monomials(n, d):  # exponent tuples in n vars of total degree <= d
    out = []
    for deg in range(d + 1):
        for c in itertools.combinations_with_replacement(range(n), deg):
            e = [0] * n
            for i in c: e[i] += 1
            out.append(tuple(e))
    return out
def rank_mod(rows, p):
    rows = [r[:] for r in rows]
    rank = 0; ncol = len(rows[0]) if rows else 0
    for col in range(ncol):
        piv = None
        for i in range(rank, len(rows)):
            if rows[i][col] % p: piv = i; break
        if piv is None: continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        inv = pow(rows[rank][col], p - 2, p)
        rows[rank] = [x * inv % p for x in rows[rank]]
        for i in range(len(rows)):
            if i != rank and rows[i][col] % p:
                f = rows[i][col]
                rows[i] = [(a - f * b) % p for a, b in zip(rows[i], rows[rank])]
        rank += 1
    return rank
def reg(M, s, t, p=P1):
    S = shell(M, s, t)
    n = len(S)
    if n == 0: return None, 0
    coords = [k[:-1] for k in S]  # drop last coord (determined)
    m = 0
    while True:
        mons = monomials(M - 1, m)  # degree <= m
        # matrix: rows = monomials, cols = points
        rows = []
        for e in mons:
            rows.append([ (lambda k: eval_mon(k, e, p))(k) for k in coords])
        r = rank_mod(rows, p)
        if r == n:  # polys of degree <= m interpolate S => reg <= m
            return m, n
        m += 1
def eval_mon(k, e, p):
    v = 1
    for ki, ei in zip(k, e):
        if ei: v = v * pow(ki % p, ei, p) % p
    return v
if __name__ == '__main__':
    Mmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    tmax = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    best = {}
    for M in range(2, Mmax + 1):
        for t in range(1, tmax + 1):
            for s in range(0, t + 1):
                if (t - s) % 2: continue
                r, n = reg(M, s, t)
                if r is None: continue
                ratio = r / t
                print(f"M={M} t={t} s={s} |S|={n} reg={r} ratio={ratio:.4f}", flush=True)
                if ratio > best.get(M, (0,))[0]: best[M] = (ratio, t, s, r)
        print("BEST", M, best.get(M), flush=True)
