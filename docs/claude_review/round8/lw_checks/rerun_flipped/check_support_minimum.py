"""Reach of the bounded-log (Lang-Waldschmidt) route: minimal support weight of rational row
combinations on fully flipped Hadamard-core profiles (Section 5 of flipped.md).

Every zero-sum rational row vector is lambda = H n / M with n rational on the nonconstant labels.
Its column image v = lambda^T S is (up to the orientation sign)
   unflipped copy of label a      : n_a
   flipped copy at row x          : n_{a(x)} - 2 H(x,a(x)) (Hn)_x / M.
With weight u per unflipped prime and f per flipped prime the support weight is u*A(n) + f*B(n),
A(n) = sum_{a in supp n} (#unflipped copies of a),  B(n) = #{x : M n_{a(x)} != 2 H(x,a(x)) (Hn)_x}.
We enumerate integer n (ranges below), record all achieved (A,B) pairs, and report for several
f/u the minimum of (support weight)/W, W = u*(#unflipped) + f*M, together with the best label pair.
The bounded-log route needs a value < 1/4 (times 1/kappa).  Exact integer arithmetic (numpy int64).
"""
import itertools, sys
import numpy as np
from fcommon import paley_I, sylvester


def frontier(H, b, asg, amp, maxsupp=None, chunk=1 << 20):
    H = np.array(H, dtype=np.int64)
    M = H.shape[0]
    Hn_mat = H[:, 1:]                       # M x (M-1)
    flips = np.zeros(M, dtype=np.int64)     # flips per label (index = label)
    for x in range(M):
        flips[asg[x]] += 1
    unf = (b - flips)[1:]                   # unflipped copies per nonconstant label
    assert (unf >= 1).all()
    ax = np.array(asg) - 1                  # column index (0-based among nonconstant) of a(x)
    sig = H[np.arange(M), np.array(asg)]    # H(x, a(x))
    seen = {}
    vals = list(range(-amp, amp + 1))
    L = M - 1

    def process(N):
        N = N[np.any(N != 0, axis=1)]       # n = 0 is not a character
        if N.shape[0] == 0:
            return
        Hn = N @ Hn_mat.T                   # (#vec, M)
        A = (N != 0).astype(np.int64) @ unf
        coef = M * N[:, ax] - 2 * sig[None, :] * Hn
        B = np.count_nonzero(coef, axis=1)
        key = A * 1000 + B
        uk, idx = np.unique(key, return_index=True)
        for k, i in zip(uk.tolist(), idx.tolist()):
            if k not in seen:
                seen[k] = N[i].copy()

    if maxsupp is None:
        # full box {-amp..amp}^L, enumerated in chunks via mixed radix
        base = 2 * amp + 1
        total = base ** L
        start = 1
        while start < total:
            stop = min(total, start + chunk)
            idx = np.arange(start, stop, dtype=np.int64)
            N = np.zeros((stop - start, L), dtype=np.int64)
            t = idx.copy()
            for k in range(L):
                N[:, k] = t % base - amp
                t //= base
            process(N)
            start = stop
    else:
        nz = [v for v in vals if v != 0]
        buf = []
        for s in range(1, maxsupp + 1):
            for supp in itertools.combinations(range(L), s):
                for co in itertools.product(nz, repeat=s):
                    v = np.zeros(L, dtype=np.int64)
                    v[list(supp)] = co
                    buf.append(v)
                    if len(buf) >= chunk:
                        process(np.array(buf))
                        buf = []
        if buf:
            process(np.array(buf))
    pts = sorted((k // 1000, k % 1000) for k in seen)
    return pts, seen, int(unf.sum())


def report(name, H, b, asg, amp, maxsupp=None):
    M = len(H)
    pts, seen, nunf = frontier(H, b, asg, amp, maxsupp)
    print(f"--- {name} (M={M}, b={b}, entries |n|<={amp}" + (f", supp<={maxsupp}" if maxsupp else ", full box") + ")")
    # Pareto frontier
    par = []
    for A, B in pts:
        if not any(A2 <= A and B2 <= B and (A2, B2) != (A, B) for A2, B2 in pts):
            par.append((A, B))
    print("    Pareto frontier (A = #unflipped primes, B = #flipped primes in support):", par)
    for rho in [0.25, 0.5, 1, 2, 3, 4, 6, 10]:
        W = nunf + rho * M
        best = min(pts, key=lambda p: p[0] + rho * p[1])
        n = seen[best[0] * 1000 + best[1]]
        s = int(np.count_nonzero(n))
        print(f"    f/u={rho:5.2f}: W_F/W={rho*M/W:5.3f}  min support/W={(best[0]+rho*best[1])/W:6.3f}"
              f"  at (A,B)={best}, |supp n|={s}, n-entries={sorted(set(abs(int(t)) for t in n if t))}")


def main():
    b = 5
    H = paley_I(11)
    report("Paley-12 canonical", H, b, [1] + list(range(1, 12)), amp=2)
    asg = [1, 3, 5, 7, 9, 11, 2, 4, 6, 8, 10, 3]
    report("Paley-12 other cap-2", H, b, asg, amp=2)
    H = sylvester(3)
    report("Walsh-8 canonical", H, b, [1] + list(range(1, 8)), amp=3)
    H = sylvester(4)
    report("Walsh-16 canonical", H, b, [1] + list(range(1, 16)), amp=1)
    H = paley_I(19)
    report("Paley-20 canonical", H, b, [1] + list(range(1, 20)), amp=1, maxsupp=6)
    report("Paley-20 canonical", H, b, [1] + list(range(1, 20)), amp=2, maxsupp=3)
    print("DONE")


if __name__ == "__main__":
    main()
