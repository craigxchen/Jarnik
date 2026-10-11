"""Concrete skew double-core fake circles (Theorem 5.1 of exact.md).

usage: python3 instance_skew.py spec C Q0mode [refine_nmax]
  spec  : q (Paley, M = q+1) or s<t> (doubling, M = 2^t)
  C     : arc constant (rational, <= 1)
  Q0mode: 'A0' (no residue data, L = 1) or an integer f (residue moduli Q0 = f*M; L = M-th prime
          = 1 mod 4 exceeding Q0)

EXACT (integers, 50-digit Decimal logarithms):
  * P_* = max(160000 M^4 L^2 / C^2, 21 r^2), r = 2(M-1); the r smallest primes = 1 mod 4 above P_*,
    the log-width of their window (<= 1/(4r) required), W = sum log p, W/(M log M), 1.33(Pi - 1) < 1;
  * the theorem's union bound  (4(pi+1)e^(1/16)/C) * sum_{n>=2} n x^n,  x = 2 M L^(1/2) e^(-tau/4),
    evaluated exactly (must be <= 1/2);
  * optional refinement: the union bound with the exact slacks s_c = sum_S w_S(|c(S)| - 1/2) of every
    integral character with ||c||_1 <= refine_nmax (one of each +-pair), plus the theorem's tail for
    larger n; then the smallest P (primes re-chosen) for which the refined bound is <= 1/2 is searched
    by bisection on the exponent, giving a smaller W for the same profile.
"""
import sys, itertools, math
from fractions import Fraction
from decimal import Decimal, getcontext
import numpy as np
from skew_double import paley_skew, kron_skew, cuts

getcontext().prec = 50
PI = Decimal('3.14159265358979323846264338327950288419716939937510')


def is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d = n - 1; s = 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True      # deterministic for n < 3.3e24


def primes_1mod4_from(start, count):
    out = []
    n = start + (1 - start) % 4
    while len(out) < count:
        if is_prime(n):
            out.append(n)
        n += 4
    return out


def build(spec):
    K = kron_skew(int(spec[1:])) if spec.startswith('s') else paley_skew(int(spec))
    M = K.shape[0]
    S, S2 = cuts(K)
    P = np.array([[1 if x in s else 0 for x in range(M)] for s in S + S2], dtype=np.int64)
    return M, P


def theorem_bound(M, L, C, tau):
    x = 2 * M * Decimal(L).sqrt() * (-tau / 4).exp()
    pre = 4 * (PI + 1) * Decimal('0.0625').exp() / Decimal(C.numerator) * Decimal(C.denominator)
    if x >= 1:
        return Decimal('Infinity'), x
    s = x * x * (2 - x) / (1 - x) ** 2          # sum_{n>=2} n x^n
    return pre * s, x


def chars_upto(M, nmax):
    """All nonzero zero-sum integer c with ||c||_1 <= nmax, one of each +-pair, as an array."""
    out = []
    def rec(i, cur, l1, s):
        if i == M - 1:
            v = -s
            if l1 + abs(v) <= nmax:
                c = cur + [v]
                if any(c):
                    first = next(t for t in c if t)
                    if first > 0:
                        out.append(c)
            return
        for t in range(-(nmax - l1), nmax - l1 + 1):
            rec(i + 1, cur + [t], l1 + abs(t), s + t)
    rec(0, [], 0, 0)
    return np.array(out, dtype=np.int64)


def refined_bound(M, P, logs, L, C, chars, tau, nmax):
    """sum over characters with ||c||_1 <= nmax of the per-character bound with exact slack (both signs:
    factor 2 since the list holds one of each +-pair), plus the theorem's tail over n > nmax."""
    CS = np.abs(chars @ P.T)                                 # |c(S)|
    n = np.abs(chars).sum(axis=1)
    W = sum(logs)
    tot = Decimal(0)
    pre = 4 * (PI + 1) * Decimal('0.0625').exp() / Decimal(C.numerator) * Decimal(C.denominator)
    pre0 = 4 * (PI + 1) / Decimal(C.numerator) * Decimal(C.denominator)
    # s_c exact:  sum_S w_S |c(S)| - W/2
    for row, nn in zip(CS, n):
        s = sum(Decimal(int(a)) * w for a, w in zip(row, logs)) - W / 2
        tot += 2 * pre0 * int(nn) * Decimal(L) ** (Decimal(int(nn)) / 2) * (-s / 2).exp()
    x = 2 * M * Decimal(L).sqrt() * (-tau / 4).exp()
    tail = Decimal(0)
    if x < 1:
        # sum_{n > nmax} n x^n
        k = nmax + 1
        tail = pre * (x ** k) * (k - (k - 1) * x) / (1 - x) ** 2
    else:
        tail = Decimal('Infinity')
    return tot + tail


def main():
    spec = sys.argv[1]; C = Fraction(sys.argv[2]); mode = sys.argv[3]
    nref = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    M, P = build(spec)
    r = 2 * (M - 1)
    if mode == 'A0':
        L = 1; Q0 = None
    else:
        Q0 = int(mode) * M
        ells = primes_1mod4_from(Q0 + 1, M)
        L = ells[-1]
    Pst = max(math.ceil(Fraction(160000 * M ** 4 * L * L) / (C * C)), 21 * r * r)
    ps = primes_1mod4_from(Pst, r)
    logs = [Decimal(p).ln() for p in ps]
    tau = logs[0]; eta = logs[-1] - logs[0]
    W = sum(logs)
    Pi = 1.0
    for p in ps:
        Pi *= 1 + 2 / (p ** 0.5 - 1)
    ub, x = theorem_bound(M, L, C, tau)
    print(f'[{spec}] M={M} r={r} C={C} Q0={Q0} L={L}: P_*={Pst:.4e}, primes {ps[0]}..{ps[-1]}, '
          f'eta={float(eta):.2e} (<= 1/(4r)={1/(4*r):.2e}: {eta <= Decimal(1)/(4*r)}), tau={float(tau):.4f}')
    print(f'   W = {float(W):.2f},  W/(M log M) = {float(W)/(M*math.log(M)):.2f},  1.33(Pi-1) = {1.33*(Pi-1):.2e},  '
          f'theorem union bound = {float(ub):.4f} (x = {float(x):.4f})')
    if nref:
        chars = chars_upto(M, nref)
        print(f'   refinement: {len(chars)} characters with ||c||_1 <= {nref} (one per +-pair)')
        lo, hi = 0.0, math.log10(Pst)
        best = None
        for it in range(30):
            mid = (lo + hi) / 2
            P0 = max(int(10 ** mid), 21 * r * r)
            ps2 = primes_1mod4_from(P0, r)
            lg = [Decimal(p).ln() for p in ps2]
            if lg[-1] - lg[0] > Decimal(1) / (4 * r):
                lo = mid; continue
            b = refined_bound(M, P, lg, L, C, chars, lg[0], nref)
            if b <= Decimal('0.5'):
                hi = mid; best = (ps2, lg, b)
            else:
                lo = mid
        if best:
            ps2, lg, b = best
            W2 = sum(lg)
            print(f'   refined: primes {ps2[0]}..{ps2[-1]}, W = {float(W2):.2f}, W/(M log M) = {float(W2)/(M*math.log(M)):.2f}, '
                  f'refined union bound = {float(b):.4f}')


if __name__ == '__main__':
    main()
