"""Exact checks for Section 6.4 of t1rigour.md (single-cut ternary cubes; projection-cheap components of TC).

(1) rho(Lambda_t(s)) = s - |t|  for the ternary slice Lambda_t(s) = {k in {-1,0,1}^s : sum k = t}.
(2) rho(A2(s,j)) = s + j for A2(s,j) = {k in {-1,0,1}^s : |sum k| <= j}   (upper bound: Lemma 6.7; lower: zonotope).
    rho(L) = least e such that every function on L is a polynomial of degree <= e
           = least e with rank[ k^alpha ]_{k in L, |alpha|<=e, alpha_x<=2} = |L|.
    Ranks over Q are computed exactly with flint fmpz_mat.
(3) Counterexample model: d = 6, S' = {0,1,2,3}, j = 1, A = {k in {-1,0,1}^6 : |k_0+k_1+k_2+k_3| <= 1}.
    For m = 6, 7, 8: exact basis of V(A,m); check Phi_c(u) = sum_k c_k <k,u>^m == 0 for u supported on S'
    (generic rational), and != 0 for some c at a generic u of full support.
"""
import itertools, sys
from flint import fmpz_mat, fmpq

def monos(s, e):
    return [a for a in itertools.product(range(3), repeat=s) if sum(a) <= e]

def evalmat(L, s, e):
    ms = monos(s, e)
    rows = []
    for k in L:
        row = []
        for a in ms:
            v = 1
            for x in range(s):
                if a[x]:
                    v *= k[x] ** a[x]
            row.append(v)
        rows.append(row)
    return fmpz_mat(rows)

def rho(L, s):
    for e in range(0, 2 * s + 1):
        if evalmat(L, s, e).rank() == len(L):
            return e
    raise RuntimeError

out = []
ok = True
for s in range(1, 7):
    cube = list(itertools.product((-1, 0, 1), repeat=s))
    for t in range(0, s + 1):
        Lt = [k for k in cube if sum(k) == t]
        r = rho(Lt, s)
        good = (r == s - t)
        ok &= good
        out.append(f"slice s={s} t={t}: |L|={len(Lt)} rho={r} (s-|t|={s-t}) {'OK' if good else 'MISMATCH'}")
    for j in range(0, s + 1):
        A2 = [k for k in cube if abs(sum(k)) <= j]
        r = rho(A2, s)
        good = (r == min(s + j, 2 * s))
        ok &= good
        out.append(f"A2 s={s} j={j}: |A2|={len(A2)} rho={r} (s+j={s+j}) {'OK' if good else 'MISMATCH'}")
print("\n".join(out))
print("ALL rho formulas:", "OK" if ok else "FAIL")

# (3) counterexample model
d, Sp, j = 6, [0, 1, 2, 3], 1
A = [k for k in itertools.product((-1, 0, 1), repeat=d) if abs(sum(k[x] for x in Sp)) <= j]
print(f"\nmodel: d={d}, S'={Sp}, j={j}, |A|={len(A)}, rho(A) = {rho(A, d)}")

def all_monos(dim, e):
    res = []
    for tot in range(e + 1):
        for a in itertools.product(range(tot + 1), repeat=dim):
            if sum(a) == tot:
                res.append(a)
    return res

def V_basis(A, d, m):
    ms = all_monos(d, m - 1)
    M = []
    for a in ms:
        row = []
        for k in A:
            v = 1
            for x in range(d):
                if a[x]:
                    v *= k[x] ** a[x]
            row.append(v)
        M.append(row)
    M = fmpz_mat(M)
    X, nul = M.nullspace()
    cols = [[int(X[i, c]) for i in range(X.nrows())] for c in range(nul)]
    return cols

def Phi(c, A, u, m):
    tot = 0
    for ck, k in zip(c, A):
        if ck:
            tot += ck * sum(fmpq(k[x]) * u[x] for x in range(len(u))) ** m
    return tot

u_S = [fmpq(1), fmpq(3), fmpq(7), fmpq(19), fmpq(0), fmpq(0)]
u_S1 = [fmpq(1), fmpq(3), fmpq(7), fmpq(19), fmpq(5), fmpq(0)]
u_gen = [fmpq(1), fmpq(3), fmpq(7), fmpq(19), fmpq(43), fmpq(101)]
for m in (6, 7, 8, 9):
    B = V_basis(A, d, m)
    zS = all(Phi(c, A, u_S, m) == 0 for c in B)
    zS1 = all(Phi(c, A, u_S1, m) == 0 for c in B)
    zG = all(Phi(c, A, u_gen, m) == 0 for c in B)
    print(f"m={m}: dim V(A,m)={len(B)}; Phi_c(u_S')==0 for all c: {zS}; "
          f"u on S' + one coordinate: {zS1}; generic full-support u: {zG}")
