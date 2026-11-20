#!/usr/bin/env python3
"""Judge's own check of X.3 (mod-four lower bound on owner twist weights) on the IX.3 Legendre design.
For q = 11, 19, 23, full positive class S. Lower bound: parity count (proved argument).
Upper bound: explicit integral V with D V = (H - D)/2 and owner weight equal to the parity bound,
found by search with an own Hermite-form integer solver. Unit weights (log p_k = 1)."""
import itertools, math, sys

def chi(a, q):
    a %= q
    if a == 0: return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1

def egcd(a, b):
    if b == 0: return (a, 1, 0) if a >= 0 else (-a, -1, 0)
    g, x, y = egcd(b, a % b)
    return g, y, x - (a // b) * y

def solve_int(A, c):
    """Integer solution y of A y = c, or None. Column Hermite reduction."""
    r = len(A); n = len(A[0]) if r else 0
    M = [row[:] for row in A]
    U = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    def colop(j1, j2, a, b, cc, d):  # new col j1 = a*c1 + b*c2, col j2 = cc*c1 + d*c2
        for Mat in (M, U):
            for row in Mat:
                x1, x2 = row[j1], row[j2]
                row[j1] = a * x1 + b * x2; row[j2] = cc * x1 + d * x2
    piv = 0; pivots = []
    for i in range(r):
        if piv >= n: pivots.append(None); continue
        for j in range(piv + 1, n):
            if M[i][j] != 0:
                a1, a2 = M[i][piv], M[i][j]
                g, s, t = egcd(a1, a2)
                # [s t; -a2/g a1/g] unimodular
                colop(piv, j, s, t, -a2 // g, a1 // g)
        if M[i][piv] != 0:
            pivots.append(piv); piv += 1
        else:
            pivots.append(None)
    z = [0] * n
    for i in range(r):
        acc = sum(M[i][j] * z[j] for j in range(n))
        if pivots[i] is None:
            if acc != c[i]: return None
        else:
            j = pivots[i]
            acc -= M[i][j] * z[j]
            rem = c[i] - acc
            if rem % M[i][j] != 0: return None
            z[j] = rem // M[i][j]
    y = [sum(U[a][b] * z[b] for b in range(n)) for a in range(n)]
    assert all(sum(A[i][k] * y[k] for k in range(n)) == c[i] for i in range(r))
    return y

def run(q):
    labels = list(range(1, q))           # column 0 deleted
    S = [i for i in labels if chi(0 - i, q) == 1]   # K_i0 = chi(0 - i) = 1
    m = len(S); i0 = S[0]
    a = {i: [1 + chi(k - i, q) for k in labels] for i in labels}
    idx = {k: t for t, k in enumerate(labels)}
    Drows = []; Hrows = []
    for i in S[1:]:
        d = [a[i][t] - a[i0][t] for t in range(q - 1)]
        h = [0] * (q - 1); h[idx[i]] += 1; h[idx[i0]] -= 1
        assert sum(abs(x) for x in d) == q + 1 and max(abs(x) for x in d) <= 2
        assert all((d[t] - h[t]) % 2 == 0 for t in range(q - 1))
        Drows.append(d); Hrows.append(h)
    Bm = [[(Hrows[r][t] - Drows[r][t]) // 2 for t in range(q - 1)] for r in range(m - 1)]
    owners = [idx[i] for i in S]; nonown = [t for t in range(q - 1) if t not in owners]
    DS = [[row[t] for t in owners] for row in Drows]
    DN = [[row[t] for t in nonown] for row in Drows]
    total_parity = 0; total_found = 0
    for k in range(q - 1):
        # parity: V_ik - V_i0k = B_rk mod 2, i = S[r+1]
        par = [0] + [Bm[r][k] % 2 for r in range(m - 1)]
        best = None
        for c in (0, 1):
            pat = [(c + x) % 2 for x in par]
            w = sum(pat)
            if best is None or w < best[0]: best = (w, [pat])
            elif w == best[0]: best[1].append(pat)
        total_parity += best[0]
        found = None
        for pat in best[1]:
            odd = [t for t in range(m) if pat[t] == 1]
            for signs in itertools.product((1, -1), repeat=len(odd)):
                xS = [0] * m
                for t, sg in zip(odd, signs): xS[t] = sg
                rhs = [Bm[r][k] - sum(DS[r][u] * xS[u] for u in range(m)) for r in range(m - 1)]
                y = solve_int(DN, rhs)
                if y is not None:
                    found = xS; break
            if found: break
        if found is None:
            print(f"  q={q} column {k}: parity bound NOT attained with entries in {{0,+-1}}")
        else:
            total_found += sum(abs(x) for x in found)
    T = q - 1; B = q - 1
    cs = 0.5 * (m * T - math.sqrt(m * q * B) - m)
    print(f"q={q} m={m}: parity lower bound {total_parity} (= m(q-3)/2 = {m*(q-3)//2}); "
          f"explicit integral V with owner weight {total_found}; Astra CS bound {cs:.1f}")

if __name__ == "__main__":
    for q in (11, 19, 23):
        run(q)
