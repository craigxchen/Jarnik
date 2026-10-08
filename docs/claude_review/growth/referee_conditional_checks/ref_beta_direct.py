"""Referee check R1: independent computation of h^0(Bl_Y X_M, NH - mE) for small (M,N).

Independent of the author's A_s (torus-weight) reduction.  We work directly in the REAL circle
model X_M = {x_1^2+y_1^2 = ... = x_M^2+y_M^2} subset P^(2M-1) and parametrize a neighbourhood of
the diagonal line Y by
    z_1 = (1, y)           (affine chart x_1 = 1 of the direction of z_1),
    z_j = z_1 * c(s_j),    c(s) = ((1-s^2) + 2 s i)/(1+s^2)  (unit circle, c(0)=1),   j = 2..M,
so Y = {s_2 = ... = s_M = 0} and (y, s_2..s_M) are local coordinates (dc/ds(0) = 2i != 0).
A degree-N form F restricts to  G = F(z_1, z_1 c_2, ...) * prod_j (1+s_j^2)^N, a polynomial in
(y, s).  Because the parametrization is dominant onto the chart, R_N = H^0(X_M, O(N)) (X_M is
projectively normal) injects into these polynomials.  ord_Y F = least total s-degree occurring in G.
Hence  h^0(NH - mE) = rank{G} - rank{G truncated to s-degree < m}.
Compared against beta_vanishing.beta (author's code) for the totals h^0(NH) and sum_m h^0(NH-mE).
Arithmetic mod a large prime (two primes for safety).
"""
import sys, itertools
sys.path.insert(0, __file__.rsplit('/', 1)[0] + '/rerun')

def polymul(P, Q, p):
    out = {}
    for a, c in P.items():
        for b, d in Q.items():
            k = tuple(x + y for x, y in zip(a, b))
            out[k] = (out.get(k, 0) + c * d) % p
    return {k: v for k, v in out.items() if v}

def polypow(P, e, nv, p):
    R = {tuple([0] * nv): 1}
    for _ in range(e):
        R = polymul(R, P, p)
    return R

def build(M, N, p):
    nv = M  # variables: y (index 0), s_2..s_M (indices 1..M-1)
    def var(i, coef=1, power=1):
        e = [0] * nv; e[i] = power
        return {tuple(e): coef % p}
    def add(*Ps):
        out = {}
        for P in Ps:
            for k, v in P.items():
                out[k] = (out.get(k, 0) + v) % p
        return {k: v for k, v in out.items() if v}
    one = {tuple([0] * nv): 1}
    y = var(0)
    # coordinates as numerators; pair j>=2 carries denominator (1+s_j^2)
    X = [one]; Yc = [y]
    for j in range(1, M):
        s = var(j); s2 = var(j, 1, 2)
        ys = polymul(y, s, p); ys2 = polymul(y, s2, p)
        X.append(add(one, {k: (-v) % p for k, v in s2.items()}, {k: (-2 * v) % p for k, v in ys.items()}))
        Yc.append(add({k: (2 * v) % p for k, v in s.items()}, y, {k: (-v) % p for k, v in ys2.items()}))
    onep = [None] + [add(one, var(j, 1, 2)) for j in range(1, M)]
    coords = []
    for j in range(M):
        coords.append((j, X[j])); coords.append((j, Yc[j]))
    Gs = []
    for mono in itertools.combinations_with_replacement(range(2 * M), N):
        G = one
        deg_pair = [0] * M
        for c in mono:
            j, P = coords[c]
            G = polymul(G, P, p); deg_pair[j] += 1
        for j in range(1, M):
            G = polymul(G, polypow(onep[j], N - deg_pair[j], nv, p), p)
        Gs.append(G)
    return Gs

def rank_mod(vecs, p):
    # vecs: list of dicts key->val ; gaussian elimination with dict rows
    basis = {}  # pivot key -> row (dict), row normalized pivot=1
    r = 0
    for v in vecs:
        v = dict(v)
        while v:
            k = max(v)
            if k in basis:
                c = v[k]; row = basis[k]
                for kk, vv in row.items():
                    v[kk] = (v.get(kk, 0) - c * vv) % p
                    if v[kk] == 0: del v[kk]
            else:
                inv = pow(v[k], p - 2, p)
                basis[k] = {kk: vv * inv % p for kk, vv in v.items()}
                r += 1
                break
    return r

def direct(M, N, p):
    Gs = build(M, N, p)
    h0 = rank_mod(Gs, p)
    tot = 0; m = 1; hist = []
    while True:
        tr = [{k: v for k, v in G.items() if sum(k[1:]) < m} for G in Gs]
        hm = h0 - rank_mod(tr, p)
        hist.append(hm)
        if hm == 0: break
        tot += hm; m += 1
    return h0, tot, hist

if __name__ == '__main__':
    from beta_vanishing import beta
    for arg in sys.argv[1:]:
        M, N = map(int, arg.split(','))
        res = [direct(M, N, p) for p in (1000003, 999983)]
        assert res[0] == res[1], res
        h0, tot, hist = res[0]
        T, A, b = beta(M - 1, N)
        print(f"M={M} N={N}: direct h0={h0} sum_m h0(NH-mE)={tot}  | author A_s: h0={A} sumT={T}  "
              f"agree={h0 == A and tot == T}  beta_N={tot/(N*h0):.4f}  h0(NH-mE), m=1..: {hist[:-1]}", flush=True)
