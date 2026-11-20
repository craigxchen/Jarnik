"""Checks of the mod-3 support statement used for full column rank of K_{S^c,T} (VIII.4).
(a) q = 11, 23 (3 splits in Q(sqrt(-q)) since q = 2 mod 3): enumerate the F_3 eigenspaces of the Paley
    conference matrix C (order q+1) for eigenvalues g, -g with g^2 = -q mod 3; minimum weight vs sqrt(q).
(b) q = 127, 383: for all nonzero v in F_3^q with |supp v| <= 3 (affine symmetry: translate and dilate so that
    supp contains 0 and, if |supp|>=2, 1), the support of Kv mod 3 vs the bound sqrt(q) - |supp v| - 1.
"""
import itertools, sys
import numpy as np


def paley(q):
    sq = {x * x % q for x in range(1, q)}
    chi = [0] + [1 if x in sq else -1 for x in range(1, q)]
    K = np.array([[chi[(j - i) % q] for j in range(q)] for i in range(q)], dtype=np.int64)
    C = np.zeros((q + 1, q + 1), dtype=np.int64)
    C[:q, :q] = K
    C[:q, q] = 1
    C[q, :q] = -1
    return K, C


def nullspace_mod3(M):
    M = M.copy() % 3
    rows, cols = M.shape
    piv = []
    r = 0
    for c in range(cols):
        nz = [k for k in range(r, rows) if M[k, c]]
        if not nz:
            continue
        M[[r, nz[0]]] = M[[nz[0], r]]
        M[r] = (M[r] * M[r, c]) % 3  # inverse of 1 is 1, of 2 is 2
        for k in range(rows):
            if k != r and M[k, c]:
                M[k] = (M[k] - M[k, c] * M[r]) % 3
        piv.append(c)
        r += 1
    free = [c for c in range(cols) if c not in piv]
    basis = []
    for f in free:
        v = np.zeros(cols, dtype=np.int64)
        v[f] = 1
        for k, c in enumerate(piv):
            v[c] = (-M[k, f]) % 3
        basis.append(v)
    return np.array(basis)


def part_a(q):
    K, C = paley(q)
    assert ((C @ C + q * np.eye(q + 1, dtype=np.int64)) == 0).all()
    gs = [g for g in (1, 2) if (g * g + q) % 3 == 0]
    for g in gs:
        Bs = nullspace_mod3(C - g * np.eye(q + 1, dtype=np.int64))
        dim = len(Bs)
        best = q + 2
        for coeffs in itertools.product(range(3), repeat=dim):
            if not any(coeffs):
                continue
            w = np.count_nonzero((np.array(coeffs) @ Bs) % 3)
            best = min(best, w)
        print(f"q={q} eigenvalue g={g} mod 3: eigenspace dim {dim}, min weight {best}, sqrt(q)={q**.5:.2f}")


def part_b(q):
    K, _ = paley(q)
    worst = {}
    for k in (1, 2, 3):
        best = q
        if k == 1:
            shapes = [((0,), c) for c in [(1,)]]  # +-e_0 (sign irrelevant)
        elif k == 2:
            shapes = [((0, 1), c) for c in [(1, 1), (1, 2)]]
        else:
            shapes = [((0, 1, x), c) for x in range(2, q) for c in itertools.product((1, 2), repeat=3) if c[0] == 1]
        for supp, c in shapes:
            v = np.zeros(q, dtype=np.int64)
            for s_, cc in zip(supp, c):
                v[s_] = cc
            w = np.count_nonzero((K @ v) % 3)
            best = min(best, w)
        worst[k] = best
        print(f"q={q} |supp v|={k}: min |supp(Kv mod 3)| = {best}; bound sqrt(q)-k-1 = {q**.5-k-1:.2f}")


if __name__ == '__main__':
    for q in (11, 23):
        part_a(q)
    for q in (127, 383):
        part_b(q)
