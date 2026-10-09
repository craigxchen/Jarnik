"""Height bookkeeping for the two-logarithm forms on literal fully flipped Paley / Walsh profiles
(Section 4-5 of flipped.md).  Not an endpoint construction: these profiles are NOT on short arcs;
the script evaluates, for actual distinct split primes, the heights that enter

  * the single-label form  M Log u_a - 2 Log v_a      : L1(a)   = W_a + W_F
  * the label-pair form    M Log u_ab - 4 Log w_ab    : L2(a,b) = W_a + W_b + F_ab
  * the column character (unconditional necessary)    : (M/2) W_a + W_F - W/2   (must be >= -O(log MC))
  * Baker-type product                                : W_a * W_F

and prints kappa1* = (W/4)/min L1, kappa2* = (W/4)/min L2 (a hypothesis H_kappa with kappa < kappa* and
O(log M) error terms would exclude the profile), and min_a W_a W_F / W (Baker needs this < 1/(4 c_0)).
Logs are floating point (only ratios are reported); primes are exact (deterministic Miller-Rabin).
"""
import math, sys
import numpy as np
from fcommon import split_primes_from, paley_I, sylvester, is_prime


def canonical_assignment(M):
    return [1] + list(range(1, M))


def profile_weights(H, b, flipped_window_start, unflipped_window_start):
    """b copies per label; copy 0 of label assignment[x] is flipped at row x (copy 1 for the doubled label).
    Flipped columns get primes from one window, unflipped from another (both consecutive split primes)."""
    M = len(H)
    asg = canonical_assignment(M)
    nflip = M
    nunf = b * (M - 1) - M
    fp = split_primes_from(flipped_window_start, nflip)
    up = split_primes_from(unflipped_window_start, nunf + nflip + 5)
    up = [p for p in up if p not in set(fp)][:nunf]
    assert len(set(fp) | set(up)) == nflip + nunf
    Wa = np.zeros(M)
    f = np.zeros(M)
    used = {}
    ui = 0
    flipped_slots = {}
    for x in range(M):
        a = asg[x]
        c = used.get(a, 0)
        used[a] = c + 1
        flipped_slots[(a, c)] = x
    for a in range(1, M):
        for c in range(b):
            if (a, c) in flipped_slots:
                x = flipped_slots[(a, c)]
                w = math.log(fp[x])
                f[x] = w
            else:
                w = math.log(up[ui])
                ui += 1
            Wa[a] += w
    W = Wa.sum()
    return np.array(H, dtype=float), Wa, f, W


def analyse(name, H, b, fstart, ustart):
    Hm, Wa, f, W = profile_weights(H, b, fstart, ustart)
    M = len(H)
    WF = f.sum()
    labels = np.arange(1, M)
    L1 = Wa[1:] + WF
    # F_ab = sum_x f_x [H(x,a) != H(x,b)] = (WF - (H^T diag f H)_{ab}) / 2
    G = Hm.T @ (f[:, None] * Hm)
    Fab = (WF - G) / 2.0
    L2 = Wa[:, None] + Wa[None, :] + Fab
    iu = np.triu_indices(M, k=1)
    mask = (iu[0] >= 1)
    L2v = L2[iu][mask]
    col_char = (M / 2) * Wa[1:] + WF - W / 2
    baker = (Wa[1:] * WF).min() / W
    print(f"{name:10s} M={M:4d} W={W:9.1f} WF/W={WF/W:6.3f} "
          f"kappa1*={W/4/L1.min():6.3f} kappa2*={W/4/L2v.min():6.3f} "
          f"min col-char margin/W={col_char.min()/W:7.4f} "
          f"min_a W_a*W_F/W={baker:7.2f} (log M={math.log(M):5.2f})")
    return W / 4 / L2v.min()


def main():
    b = 5
    print("== nearby primes: all primes consecutive split primes >= M^4 (Prop P model) ==")
    for q in [11, 19, 23, 43, 67, 83, 107, 131, 163, 199]:
        H = paley_I(q)
        M = q + 1
        analyse(f"Paley-{M}", H, b, M ** 4, M ** 4 + 1)
    for t in [3, 4, 5, 6, 7]:
        H = sylvester(t)
        M = 2 ** t
        analyse(f"Walsh-{M}", H, b, M ** 4, M ** 4 + 1)
    print("== heavy flips: flipped primes >= M^(4 lam), unflipped >= M^4 ==")
    for lam in [2, 3, 4, 6]:
        for q in [43, 107, 199]:
            H = paley_I(q)
            M = q + 1
            analyse(f"P{M},l={lam}", H, b, M ** (4 * lam), M ** 4)
    print("DONE")


if __name__ == "__main__":
    main()
