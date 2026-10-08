#!/usr/bin/env python3
"""Exact exhaustive checker for the chain lemma (proof.md, Lemma 4.1).

For a distribution n_0..n_e (n_t >= 0, sum = M) of M points into exponent levels
0..e, put S_tau = n_0 + ... + n_(tau-1)  (tau = 1..e), lambda = 1/(p-1), and

    Phi = sum_(tau=1)^e (S_tau - M/2)^2 + lambda * sum_(t=0)^e n_t^2 .

Chain lemma:  Phi >= c(lambda) M^2/2,   c(lambda) = (sqrt(1+4 lambda)-1)/2,
i.e. c is the positive root of c^2 + c = lambda.  All comparisons below are
exact (fractions.Fraction); "Phi >= c M^2/2" is tested via the equivalent exact
condition  y^2 + y >= lambda  with  y = 2 Phi / M^2 >= 0.

Checks performed
  (A) brute force over ALL compositions (empty levels allowed) for small M, e,
      several p (and some lambda not of the form 1/(p-1));
  (B) an exact dynamic programme over (tau, S_tau) for larger M, e, cross-checked
      against (A) where both run;
  (C) the exact continuous optimum Phi_cont(M, e, lambda) (real S_tau), obtained
      by solving the tridiagonal stationarity system exactly; checks
      integer min >= Phi_cont(e) >= c M^2/2 and that inf_e Phi_cont(e) -> c M^2/2
      (sharpness of the constant c);
  (D) the recursion gamma_0 = lambda, gamma_n = lambda(1+gamma_(n-1))/(1+gamma_(n-1)+lambda)
      from the proof: gamma_n decreases and stays >= c;
  (E) the integer statement actually used in the proof:
          min over distributions of  sum_t E(n_t, p-1) + (1/2) sum_tau (S_tau - M/2)^2
            >= c M^2/4 - M/2,
      where E(n, q) is the exact minimal number of colliding pairs (formula (3.5));
  (F) the restricted minimum Psi_p(M) used in the configuration-free finite
      inequality (levels 0 and e nonempty, all e >= 1), via a DP over strictly
      increasing S-sequences, including prime-power depths; checks Psi_p(M) >= c M^2/4 - M/2.
"""
from fractions import Fraction as Fr
from itertools import combinations
from math import comb, sqrt
import sys


def c_of(lam):
    return (sqrt(1 + 4 * float(lam)) - 1) / 2


def ge_cM2over2(Phi, M, lam):
    """exact test Phi >= c(lam) M^2 / 2."""
    y = Fr(2) * Phi / (M * M)
    return y * y + y >= lam


def E(n, q):
    """minimal number of equal-class unordered pairs for n objects in <= q classes."""
    if n <= 0:
        return 0
    b, r = divmod(n, q)
    return q * comb(b, 2) + r * b


def compositions(M, parts):
    """all (n_0..n_{parts-1}) of nonnegative integers summing to M (stars and bars)."""
    for bars in combinations(range(M + parts - 1), parts - 1):
        prev = -1
        out = []
        for b in bars:
            out.append(b - prev - 1)
            prev = b
        out.append(M + parts - 1 - prev - 1)
        yield out


def phi_of(ns, M, lam):
    half = Fr(M, 2)
    S = 0
    tot = Fr(0)
    for t in range(len(ns) - 1):  # tau = t+1 = 1..e
        S += ns[t]
        tot += (S - half) ** 2
    tot += lam * sum(n * n for n in ns)
    return tot


def bracket_of(ns, M, q):
    """sum_t E(n_t, q) + (1/2) sum_tau (S_tau - M/2)^2   (depth one, q = p-1)."""
    half = Fr(M, 2)
    S = 0
    tot = Fr(0)
    for t in range(len(ns) - 1):
        S += ns[t]
        tot += (S - half) ** 2
    return tot / 2 + sum(E(n, q) for n in ns)


def brute_min(M, e, lam, q=None):
    best_phi = None
    best_br = None
    arg = None
    for ns in compositions(M, e + 1):
        ph = phi_of(ns, M, lam)
        if best_phi is None or ph < best_phi:
            best_phi, arg = ph, ns
        if q is not None:
            br = bracket_of(ns, M, q)
            if best_br is None or br < best_br:
                best_br = br
    return best_phi, arg, best_br


def dp_min(M, e, lam):
    """exact min of Phi over integer nondecreasing S_0=0 <= S_1 <= ... <= S_e <= S_{e+1}=M."""
    half = Fr(M, 2)
    INF = None
    f = [INF] * (M + 1)
    f[0] = Fr(0)
    for tau in range(1, e + 1):
        g = [INF] * (M + 1)
        for S in range(M + 1):
            best = INF
            for Sp in range(S + 1):
                if f[Sp] is None:
                    continue
                v = f[Sp] + lam * (S - Sp) ** 2
                if best is None or v < best:
                    best = v
            if best is not None:
                g[S] = best + (S - half) ** 2
        f = g
    best = None
    for Sp in range(M + 1):
        if f[Sp] is None:
            continue
        v = f[Sp] + lam * (M - Sp) ** 2
        if best is None or v < best:
            best = v
    return best


def cont_opt(M, e, lam):
    """exact minimum over REAL x_1..x_e of sum x_tau^2 + lam sum (x_{t+1}-x_t)^2,
    x_0 = -M/2, x_{e+1} = M/2  (x_tau = S_tau - M/2).  Tridiagonal solve."""
    h = Fr(M, 2)
    if e == 0:
        return lam * (2 * h) ** 2
    a = 1 + 2 * lam  # diagonal
    b = -lam  # off-diagonal
    rhs = [Fr(0)] * e
    rhs[0] += lam * (-h)
    rhs[-1] += lam * h
    # Thomas algorithm
    cp = [Fr(0)] * e
    dp = [Fr(0)] * e
    cp[0] = b / a
    dp[0] = rhs[0] / a
    for i in range(1, e):
        den = a - b * cp[i - 1]
        cp[i] = b / den
        dp[i] = (rhs[i] - b * dp[i - 1]) / den
    x = [Fr(0)] * e
    x[-1] = dp[-1]
    for i in range(e - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    xs = [-h] + x + [h]
    val = sum(v * v for v in x) + lam * sum((xs[t + 1] - xs[t]) ** 2 for t in range(e + 1))
    return val


def psi_restricted(M, p, depths=6):
    """Psi_p(M) := min over e >= 1 and level distributions with n_0 >= 1, n_e >= 1 and
    (WLOG, see proof.md Sec. 5) no empty middle level, of
        sum_{a>=1} sum_t E(n_t,(p-1)p^(a-1)) + (1/2) sum_tau (S_tau - M/2)^2 .
    Shortest path over strictly increasing 0 < S_1 < ... < S_e < M."""
    qs = [(p - 1) * p ** (a - 1) for a in range(1, depths + 1)]
    w = [sum(E(n, q) for q in qs) for n in range(M + 1)]
    half = Fr(M, 2)
    # f[S] = min cost of a path 0 -> ... -> S (S = S_tau for some tau >= 1), including
    # level weights of the levels before S and the half-squares of the S's visited.
    f = [None] * (M + 1)
    for S in range(1, M):
        f[S] = Fr(w[S]) + (S - half) ** 2 / 2
    for S in range(1, M):
        if f[S] is None:
            continue
        for S2 in range(S + 1, M):
            v = f[S] + w[S2 - S] + (S2 - half) ** 2 / 2
            if f[S2] is None or v < f[S2]:
                f[S2] = v
    best = None
    for S in range(1, M):
        v = f[S] + w[M - S]
        if best is None or v < best:
            best = v
    return best


def main():
    ok = True
    ps = [5, 13, 17, 29, 37, 101]
    extra_lams = [Fr(1), Fr(1, 2), Fr(2), Fr(1, 1000)]
    print("=== (A) brute force over all compositions, empty levels allowed ===")
    print(" p/lam      M  e   minPhi(int)    cM^2/2    ratio   Phi_cont(e)  ok")
    worst_ratio = 1e9
    count = 0
    for lam_label, lam, q in [(f"p={p}", Fr(1, p - 1), p - 1) for p in ps] + \
                              [(f"lam={l}", l, None) for l in extra_lams]:
        c = c_of(lam)
        for M in range(1, 13):
            for e in range(1, 8):
                if comb(M + e, e) > 200000:
                    continue
                mphi, arg, mbr = brute_min(M, e, lam, q)
                count += comb(M + e, e)
                good = ge_cM2over2(mphi, M, lam)
                pc = cont_opt(M, e, lam)
                good = good and mphi >= pc and ge_cM2over2(pc, M, lam)
                if q is not None:
                    # (E) integer bracket bound: mbr >= c M^2/4 - M/2  <=>  4 mbr + 2M >= c M^2
                    t = Fr(4) * mbr + 2 * M
                    y = t / (M * M)
                    good = good and y * y + y >= lam
                ratio = float(mphi) / (c * M * M / 2)
                worst_ratio = min(worst_ratio, ratio)
                ok &= good
                if M in (4, 7, 12) and e in (1, 3, 6):
                    print(f" {lam_label:9s} {M:3d} {e:2d} {float(mphi):12.4f} {c*M*M/2:10.4f} "
                          f"{ratio:7.4f} {float(pc):11.4f}  {good}")
    print(f"compositions enumerated: {count};  min ratio minPhi/(cM^2/2) = {worst_ratio:.6f}")
    print("(A)+(C)+(E) all exact checks passed:", ok)

    print("\n=== (B) exact DP vs brute force, then larger M,e by DP ===")
    agree = True
    for p in [5, 13]:
        lam = Fr(1, p - 1)
        for M in range(1, 11):
            for e in range(1, 6):
                agree &= dp_min(M, e, lam) == brute_min(M, e, lam)[0]
    print("DP == brute force on M<=10, e<=5, p in {5,13}:", agree)
    ok &= agree
    print(" p     M   e_range   min_e minPhi_int   cM^2/2   ratio   inf_e Phi_cont  ratio")
    for p in [5, 13, 29, 101]:
        lam = Fr(1, p - 1)
        c = c_of(lam)
        for M in [20, 31, 40]:
            emax = 24 if M <= 31 else 16
            mins = []
            conts = []
            for e in range(1, emax + 1):
                v = dp_min(M, e, lam)
                pc = cont_opt(M, e, lam)
                ok &= ge_cM2over2(v, M, lam) and v >= pc and ge_cM2over2(pc, M, lam)
                mins.append(v)
                conts.append(pc)
            mi = min(mins)
            mc = min(conts)
            print(f" {p:3d} {M:4d}   1..{emax:2d}   {float(mi):14.4f} {c*M*M/2:9.4f} "
                  f"{float(mi)/(c*M*M/2):7.4f} {float(mc):14.4f} {float(mc)/(c*M*M/2):7.4f}")
    print("(B) all exact checks passed:", ok)

    print("\n=== (C) sharpness: inf over e of the continuous optimum -> c M^2/2 ===")
    for lam in [Fr(1, 4), Fr(1, 12), Fr(1, 100), Fr(1)]:
        c = c_of(lam)
        M = 2
        vals = [(e, float(cont_opt(M, e, lam)) / (c * M * M / 2)) for e in (1, 2, 3, 5, 9, 17, 33, 65)]
        print(f" lam={str(lam):6s} c={c:.6f}  Phi_cont/(cM^2/2) for e=1,2,3,5,9,17,33,65:",
              " ".join(f"{r:.6f}" for _, r in vals))

    print("\n=== (D) recursion gamma_n (proof of Lemma 4.1) ===")
    for p in [5, 13, 101]:
        lam = Fr(1, p - 1)
        c = c_of(lam)
        g = lam
        seq = [g]
        for n in range(1, 12):
            g = lam * (1 + g) / (1 + g + lam)
            seq.append(g)
        mono = all(seq[i + 1] <= seq[i] for i in range(len(seq) - 1))
        above = all(s * s + s >= lam for s in seq)  # s >= c exactly
        ok &= mono and above
        print(f" p={p:3d} c={c:.8f} gamma_0..gamma_11 = " + " ".join(f"{float(s):.8f}" for s in seq[:6]),
              "...", f"{float(seq[-1]):.8f}", " decreasing:", mono, " >=c:", above)

    print("\n=== (F) restricted minimum Psi_p(M) (prime-power depths included) ===")
    print(" p     M    Psi_p(M)    cM^2/4-M/2   E-only(p not|N)")
    for p in [5, 13, 17, 29]:
        lam = Fr(1, p - 1)
        c = c_of(lam)
        for M in [2, 3, 5, 10, 25, 50, 80]:
            ps_ = psi_restricted(M, p)
            bound = c * M * M / 4 - M / 2
            t = Fr(4) * ps_ + 2 * M
            y = t / (M * M)
            good = y * y + y >= lam
            ok &= good
            qs = [(p - 1) * p ** (a - 1) for a in range(1, 7)]
            eo = sum(E(M, q) for q in qs)
            if M in (5, 25, 80):
                print(f" {p:3d} {M:4d} {float(ps_):11.3f} {bound:12.3f} {eo:12d}   {good}")
    print("\nALL CHECKS PASSED:", ok)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
