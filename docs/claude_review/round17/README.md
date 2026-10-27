# Round 17 (partial snapshot): the Fourier-circle test CIRCLE(q)

**Status: in progress.** This is a safety snapshot taken after two container restarts. The judge's README will replace this
file. The three attacks and their referees are complete; the calibration (`calib.md`), its referee and the judge are still
running.

Results so far (refereed):
* `hlight.md`: **H_light for `sdc_q` is proved** for every prime `q = 3 mod 4`, `q >= 7`, with `c_0 = 1/6`, and `c_0 >= 0.2349` at
  every regime `q >= 3251` (Weil for 4-term characters). It is sharp up to `O(q^(-1/2))`. Two independent referee passes
  (`referee_hlight_pass1.md`, `referee_hlight_pass2.md`) uphold it with minor corrections. Consequence: barrier Thm 3.5 is an
  unconditional cap test for route 5.
* `cap.md`: no cap approximant found; exact Gauss-sum identity, `L^2` splitting of `Psi_A` on circles, and proofs that each
  relation family alone fails. Referee (`referee_cap.md`): survives, with the major correction that a cap may also come from
  `TC_{H'} \ TC_V`, which was not examined.
* `survive.md`: no survival certificate found; proved Liouville-free equivariant lower bound (Thm 3.1) and gap principle
  (Thm 3.4: on the circles `kappa_H` is 0 or `>= 2^(-m)/(|A| H q^m)`), and no product-type or sub-cube witness exists at degree
  `m'`. Referee (`referee_survive.md`): survives, one false number in Remark 4.4 repaired.

No growth improvement and no uniform bound: `M <= (2/3 + eps) log R/loglog R` is still the best bound. CIRCLE(q) is undecided
so far.

Repository note: `*.npz` arrays and files over 1MB from the scratch directory are excluded.
