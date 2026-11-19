"""Which binary quadratic forms Q = alpha X^2 + beta XY + gamma Y^2 over F_q[t] (q = 2^d)
satisfy Q(P_a) = c_a * n_a with c_a in F_q for every column P_a = (f_a, 1) of the trace construction?

Linear system over F_q in the coefficients of alpha, beta, gamma (degree <= B) and the c_a.
Reports the dimension of the solution space and, for every solution with beta != 0 (separable form),
whether some mark is isotropic (c_a = 0), which exhibits an F_q(t)-rational isotropic vector.
"""
import sys
from itertools import product
from gf2d import GF


def solve_nullspace(F, rows, ncols):
    """Row-reduce over GF(q); return basis of nullspace (list of vectors)."""
    M = [r[:] for r in rows]
    piv_cols = []
    r = 0
    for c in range(ncols):
        p = None
        for i in range(r, len(M)):
            if M[i][c]:
                p = i
                break
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        iv = F.inv(M[r][c])
        M[r] = [F.mul(iv, x) for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [x ^ F.mul(f, y) for x, y in zip(M[i], M[r])]
        piv_cols.append(c)
        r += 1
        if r == len(M):
            break
    free = [c for c in range(ncols) if c not in piv_cols]
    basis = []
    for fc in free:
        v = [0] * ncols
        v[fc] = 1
        for i, pc in enumerate(piv_cols):
            v[pc] = M[i][fc]  # char 2: minus = plus
        basis.append(v)
    return basis


def run(d, B):
    F = GF(d)
    q = F.q
    nz = list(range(1, q))
    cols = {}
    for a in nz:
        marked = [x for x in nz if F.tr(F.mul(a, x)) == 0]
        na = F.from_roots(marked)
        fa = F.add([0] * (q // 2 - 1) + [1], na)
        cols[a] = (fa, na)
    nvar = 3 * (B + 1) + len(nz)
    maxdeg = 2 * (q // 4 - 1) + B + 1
    rows = []
    for idx, a in enumerate(nz):
        fa, na = cols[a]
        fa2 = F.pmul(fa, fa)
        # coefficient of t^e in alpha*fa2 + beta*fa + gamma + c_a*na
        for e in range(maxdeg + 1):
            row = [0] * nvar
            for i in range(B + 1):
                j = e - i
                if 0 <= j < len(fa2):
                    row[i] = fa2[j]
                if 0 <= j < len(fa):
                    row[B + 1 + i] = fa[j]
            if e <= B:
                row[2 * (B + 1) + e] = 1
            if e < len(na):
                row[3 * (B + 1) + idx] = na[e]
            if any(row):
                rows.append(row)
    basis = solve_nullspace(F, rows, nvar)
    print(f"q={q}, B={B}: solution space dimension {len(basis)}")
    for v in basis:
        al = F.strip(v[:B + 1]); be = F.strip(v[B + 1:2 * B + 2]); ga = F.strip(v[2 * B + 2:3 * B + 3])
        print("   alpha deg", F.deg(al), " beta", be, " gamma deg", F.deg(ga))
    # all separable solutions (beta != 0): check one c_a vanishes
    if len(basis) <= 3:
        bad = 0
        for coeffs in product(range(q), repeat=len(basis)):
            if not any(coeffs):
                continue
            v = [0] * nvar
            for c, b in zip(coeffs, basis):
                if c:
                    v = [x ^ F.mul(c, y) for x, y in zip(v, b)]
            be = F.strip(v[B + 1:2 * B + 2])
            cs = v[3 * (B + 1):]
            al = F.strip(v[:B + 1])
            if be:
                # isotropic witnesses: some column P_a with Q(P_a)=0, or (1,0) with Q(1,0)=alpha=0
                if 0 not in cs and al:
                    bad += 1
        print("   separable solutions with no rational isotropic vector among {P_a,(1,0)}:", bad)


if __name__ == '__main__':
    for d, B in [(3, 4), (3, 8), (3, 16), (4, 8), (4, 16), (4, 24)]:
        run(d, B)
