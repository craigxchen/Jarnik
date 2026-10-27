# Referee report on round17/cap.md (CIRCLE(q) cap attempt)

Referee: adversarial, independent code only (nothing imported from `cap_checks/`). Scripts and outputs:
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round17/referee_cap_checks/`
(`regime_ref.py`, `gauss_ref.py`, `toy_ref.py`, `exact_Ni.py`, `margin_milp.py`; each `x.py` writes `x_out.txt`).

## 0. Verdict

**The headline survives:** CIRCLE(q) is undecided at every regime `q`; no cap was built and no lower bound was
proved; under round-16 README §7.4 route 5 stops; growth rate and uniform bound unchanged. Every proof in cap.md
re-derives, and every number I recomputed independently reproduces (regime arithmetic, thresholds, Prop 3.2 bounds,
Lemma 3.1 prefactor, exact trinomial ratios, slab tails, zonotope bound). Prop 4.1/4.2 (exact `L^2` splitting,
rationality of `N_i`) are correct and are additionally confirmed by an exact rational computation on a toy.

**One major error, which must be corrected before the write-up enters the record:**

* **M1 (kappa_V vs kappa_H).** Claim (1) ("a cap needs `kappa_V(hat u_theta) <= e^{-101 333}`"), Cor 3.3, §1 and §8.4
  treat `kappa_V <= T*` as *necessary* for a cap. It is only *sufficient* (`kappa_H <= kappa_V`). CIRCLE(q) is a
  `kappa_H` test, and a cap could come from `TC_{H'} \ TC_V`, which cap.md never examines. Correct form: a cap needs
  `kappa_H(hat u_theta) <= e^{-101 333}`; a cap by approximation on `A` needs `kappa_V(hat u_theta) <= e^{-101 333}`;
  all "the cap needs" statements in §3-§4 concern the second kind.

**Minor:** m1 (the "any `c_0`" floor is relative to the record's Lemma 3.4 and `C = 1`); m2 (Lemma 3.1's "no traction"
is an interpretation of an `l^1` coefficient bound); m3 (the `e^{-367}` dimension "reading" says nothing about the
special point `[psi]` and should not be used as an argument against "yes"); m4 (Remark 6.2 "proved negative
statements" rest on an `[M]` bound and on a sufficient-condition failure); m5 (the toy claim "circles have the
largest restriction ratio for `m < q`" is false at `sdc7`, `m = 5, 6`, including in the author's own table).
"Height is not the obstacle" (§6, §8.4) is a heuristic.

## 1. Regime parameters (all reproduced, independent code)

`referee_cap_checks/regime_ref.py` (own code, 50-digit `mpmath`, exact trinomials with python-flint; nothing imported
from `cap_checks/`), output `regime_ref_out.txt`.

* `q = 3251` is prime, `= 3 mod 4`, and the **first** such prime with `q >= 316 ln(8q)` (`316 ln 4r = 3212.51 <= 3251`;
  the previous candidate `q = 3203` has `316 ln 4r = 3207.8 > 3203`). `d = q`, `r = 2q = 6502`, `n = q + 1`.
* `w = 514.1998`, `j = 257`, `m' = 2057 < q`, `4j = 1028 < m'/2 = 1028.5`, `m'/(2j) = 4.002`. Since `m' = ceil(4w) >= 8j`
  always, `4j <= m'/2` holds at every regime `q` (not only at 3251).
* `q = 5479`: prime, `3 mod 4`, regime, `j = 342`, `m' = 2738 < q`. (`m'` is even here; the author's Prop 3.2 cube
  certificate then uses `alpha_x = 2` on `1369` coordinates, which my code also does; the bound `e^{-6602}` reproduces.)
* `sd k(P) = sqrt(2|P|/3) = 32.924`, `j/sd = 7.806` (`q = 3251`). Exact slab tail
  `ln Pr_T(|k(P)| > j) = -33.025`, union bound `ln(2q * tail) = -24.245`: reproduces `[e^{-33.0}, e^{-24.2}]`.
* Regime parameters are therefore really regime: `q >= 3251`, `m' < q`, `d >= 316 ln 4r`. No issue.

## 2. Claim (1): the scale, and the kappa_V / kappa_H error

**Arithmetic: correct.** With `tau_D = max((2/c_0)(2 ln(4r+5) + ln 1000), 4 ln(4r+5) + 8, 2 ln(1000 q n^2))` (barrier
Lemma 3.4, `C = 1`), the `c_0`-free terms are `48.67` and `62.337`, so `tau >= 62.337` for every `c_0`, and

| `tau` | `ln T*` (mine, `W = r tau`) | author (`W = r tau + 1/4`) |
| --- | --- | --- |
| floor `62.337` | `-101 332.92` | `-101 332.99` |
| `tau_D(0.184) = 296.09` | `-481 301.50` | `-481 301.57` |
| `tau_D(0.5) = 108.96` | `-177 121.51` | `-177 121.58` |

(The `0.06` difference is the `1/16` from `W in [r tau, r tau + 1/4]`; the author takes the worst end, which is the
right choice for a necessary condition.) `q = 5479`: `-179 354.39`, `-842 228.16`. Survival threshold
`ln N_0 > m' ln 4 + 2m' ln d + 2 ln T* = -166 545.5` reproduces.

**Major error (M1): the scale is stated for the wrong kappa.** The summary (item 1), §1 ("a cap at `q = 3251` needs
`kappa_V(hat u_theta) <= e^{-101 333}`"), Cor 3.3 ("the cap needs `kappa_V(hat u_theta) <= e^{-101 333}`") and §8.4
("Everything here is at the `kappa_V` level, which is exactly what the cap needs (`kappa_H <= kappa_V`)") all invert
the logic of `kappa_H <= kappa_V`:

* CIRCLE(q) and barrier Thm 3.5 are statements about `kappa_H` (`H = H' = e^{m'W/40}`). A cap needs
  `inf_theta kappa_H(hat u_theta) <= T*`.
* `kappa_H <= kappa_V` makes `kappa_V <= T*` **sufficient** for a cap (README §7.1: "for the cap direction `H` is
  irrelevant ... it suffices to exhibit one polynomial `P`"). It does **not** make it necessary. A cap can come from
  `TC_H \ TC_V`, i.e. from `V_{H'}` being a proper subspace of `V(A,m')` along the circle, with `kappa_V` large.
  The record says this explicitly (README §2.1 item 5: "Fakes may use `TC_H \ TC_V`, so the passage belongs to the
  barrier"; the only proved equality `kappa_V = kappa_{H*}` (referee_antitheory) has `ln H* = (E/2) ln E` with
  `E = binom(m'+d-1, d) ~ e^{3.5 * 10^3}` here, against `ln H' = 2.08 * 10^7` at the floor and `9.9 * 10^7` at
  `c_0 = 0.184`; so `V_{H'} = V(A,m')` is not known, and the codimension bound of reduction Thm 4.1 is not pointwise).
* Correct statement: **a cap needs `kappa_H(hat u_theta) <= e^{-101 333}` for some `theta` (every `c_0`); a cap
  certified by an approximant on `A` (a "`kappa_V`-cap") needs `kappa_V(hat u_theta) <= e^{-101 333}`.**
  Everything in §3-§4 of cap.md that says "the cap needs" is about `kappa_V`-caps only. The `kappa_H`-cap route
  (an integral-height obstruction: every integral `c in V(A,m')` with `|c|_inf <= H'` nearly vanishing on the
  circle) is not examined anywhere in cap.md. This does not change the verdict "undecided", but it removes the
  claim that the analysis covers the whole cap direction.
* The survival side is stated correctly (cap.md §4 "survival at the `kappa_V` level"; §6/§8.4 note that a survival
  proof needs its certificate in `V_{H'}`). The phrase "with `ln H' ~ 10^8` this is not the obstacle" (§6, §8.4) is
  a heuristic, not a proof: the existence of an integral certificate of height `<= H'` with the required value is
  exactly the open (P1)/(P2) passage (README §2.1 item 5), and should be labelled as such.

**Minor (m1): "for any value of `c_0`" is relative to the record's Lemma 3.4.** The floor `62.337` is the term
`2 ln(1000 q n^2/C)` of barrier Lemma 3.4, a *sufficient* threshold (pair-measure union bound). CIRCLE(q) as posed
fixes `tau >= tau_D`, so the statement is correct for the test as posed; it is not an intrinsic lower bound on any
possible cap (a sharper Diophantine lemma would lower `tau` and raise `T*`). It is also at `C = 1`; `C` enters only
through `O(ln C)` terms.

## 3. Claim (2): Gauss-sum identity (correct)

* Patterns: the record's `sdc_sym.py` defines `P_a = {a} u {x : chi(a - x) = 1}` and `P'_a = {a} u {x : chi(a-x) = -1}`.
  With `-1 in NQR` (`q = 3 mod 4`) these are the author's `{a} u (a + NQR)` and `{a} u (a + QR)`. Checked exactly.
* Proof re-derived: `sum_a zeta^a S_a = hat k(1)(1 + sum_{y in QR} zeta^{-y})`, `sum_{QR} zeta^{-y} = sum_{NQR} zeta^y =
  (-1 - g)/2`, so the factor is `(1 - g)/2`; `|1 - g| = sqrt(q+1)` since `g = +- i sqrt q`. The Fourier factor
  `1 + G(-b) = (1 - chi(b) g)/2` (`b != 0`), `(q+1)/2` (`b = 0`) also re-derived, so `W` is invertible over `Q`.
* `referee_cap_checks/gauss_ref.py` (own code; integer polynomial arithmetic mod `Phi_q` with python-flint): `g^2 = -q`,
  both identities (2.1), and `S_a + T_a = K + k_a`, for 25 random `k in [-5,5]^q`, all `q = 3 mod 4` primes
  `7 <= q <= 139`: all true (`gauss_ref_out.txt`). This is a finite check; the proof is what carries the claim.
* The reading "the slab bound does not constrain `hat k(1)` on `A`" is correct as stated (`2qj/sqrt(q+1) = 2.93 * 10^4 > q`).

## 4. Claim (3): slab reduction and single-family lower bounds (proofs correct; one reading overstated)

**Lemma 3.0: correct.** Modulo `pi` the `q` points `theta + 2 pi x/q` are equally spaced with gap `pi/q` (`q` odd), and an
arc of length `2 beta` contains `>= floor(2 beta q/pi)` of them. The same applies to `|c_a|/|gamma|`.

**Lemma 3.1: arithmetic correct, reading overstated (minor m2).** The bound (choose `a`, choose `2j+1` of the `m'`
factors, the rest free) is right, and `sum_{i<=q} cos(pi i/(2q)) = 2069.1508`, `ln(binom(m',2j+1) q (2069.15)^{-(2j+1)})
= -2770.23` reproduce (`q = 5479`: `-4042.72`). But the conclusion "the reduction re-expresses nothing of substance"
does not follow from an `l^1` bound on coefficients in the `S`-monomial basis: those monomials are far from
independent on `A`, the target is an `e^{-10^5}`-relative cancellation, and the ideal `(Q_a)` acts through all
multiples `h Q_a`, not only through monomial rewriting of `S_a^{2j+1}`. Lemma 3.1 proves only that *monomial
rewriting* of the slab powers changes an `e^{-2770}` share of the coefficient mass. The statement "no traction" should
be labelled an interpretation. (The author does say the joint use is the full problem, which is the correct reading.)

**Prop 3.2: correct, re-derived.** `Lambda_alpha = tensor_x Delta_x^{alpha_x}` on a product grid with `alpha_x + 1`
consecutive points kills every monomial `k^beta` with `|beta| < |alpha|` (some `beta_x < alpha_x`) and every degree-`m'`
monomial other than `k^alpha`, maps `k^alpha` to `alpha!`, so `Lambda_alpha(l_v^{m'} - P) = m'! v^alpha` with
`|Lambda_alpha|_1 = 2^{m'}`. Both supersets are legitimate: `T` is the zero set of `k_x^3 - k_x` alone; `Y =
W^{-1}({-j..j}^q)` (rational, mostly non-integral points) is the zero set of the `S`-family `Q_a` alone, and `W` is
invertible over `Q`, so degree is preserved. My recomputation: `ln kappa^{(T)} >= -4506.576` (sum of `ln cos` terms
`-86.988`), `ln kappa^{(Y)} >= -11 311.29` (`alpha = 4 x 514 + 1`), gaps to the floor threshold `96 826` and `90 022`;
`q = 5479`: `-6601.79`, `-16 417.99`. All reproduce.

**Cor 3.3: correct as a statement about approximants** (an approximant meeting `T* d^{m'}` on `A` must have error
`>= e^{-4507} d^{m'}` somewhere in `T \ A` and `>= e^{-11 312} d^{m'}` somewhere in `Y \ A`). Its first sentence
("the cap needs `kappa_V(hat u_theta) <= e^{-101 333}`") is M1: replace "the cap" by "a `kappa_V`-cap". The informal
phrase "certified using only relations true on `Y`" should be replaced by the precise form just given: an identity
`l_u^{m'} - P = sum h_a Q_a + E` with `E` small on `A` does not by itself make `E` small on `Y`.

## 5. Claim (4): apolarity and the fat-point system (correct)

* **Prop 3.4: correct.** `<f h, G> = <h, f(partial) G>` and perfectness in each degree give `(J_gen)_{m'}^perp = Inv_{m'}`;
  `<l_u^{m'}, G> = m'! G(u)`; `partial_{y_x}^3 G = 0` iff `deg_{y_x} G <= 2` iff multiplicity `>= m' - 2` at `[e_x]`;
  `(sum_{x in P} partial_x)^{2j+1} G = 0` iff `G(t 1_P + z)` has `t`-degree `<= 2j` iff multiplicity `>= m' - 2j` at `[1_P]`.
  The "base point" sentence is a necessary condition for an exact syntactic cap, and it is valid for a single `theta`
  as well: `J_gen` is translation-stable and `m' < q`, so `l_{u_theta}^{m'} in J_gen` for one `theta` forces every `M_i`,
  in particular `l_psi^{m'}`, into `J_gen tensor C` (same isotypic argument as Prop 4.1). Only the leading forms of
  `k^3 - k` and `Q_a`, `Q^-_a` are in `J_gen`; the author correctly notes that `J(A)` may be larger.
* **Prop 3.5: correct.** Modulo `(k^3)`, each `l_{1_P}^{2j+1} C[k]_{m'-2j-1}` contributes at most `H_T(m'-2j-1)`
  dimensions. Exact trinomials (mine): `ln(2q H_T(m'-2j-1)/H_T(m')) = -366.969` (`q = 3251`), `-657.674` (`q = 5479`).
  `H_A(e) > 0` for `e <= q`: `J(A)_e = C[k]_e` propagates upward, `sum_e H_A(e) = |A|`, and
  `sum_{e<q} H_T(e) = (3^q - H_T(q))/2 < 3^q/2 <= |A|` (checked exactly; `|A| >= 3^q/2` also follows directly from the
  exact union bound `mu(T \ A) <= e^{-24.2}`).
* **Minor (m3).** The "reading" (fraction `e^{-367}` of the space) is correctly labelled evidence-level in §3.3. It has no bearing on a single,
  maximally symmetric point such as `[psi]`; it should not be cited (as in §0 "Net" and §8.1) as a reason against "yes".

## 6. Claim (5): exact L^2 splitting on circles (correct; the best result of the write-up)

* **Prop 4.1, re-derived.** `l_{u_theta}^{m'} = 2^{-m'} sum_i binom(m',i) e^{i theta(2i-m')} M_i`; translation is a
  permutation of `A`, unitary on `L^2(mu_A)`, preserves `Poly_{<m'}|_A`, hence commutes with `top`; `M_i o tau =
  zeta^{2i-m'} M_i`, and the `m'+1` eigenvalues are distinct exactly when `q` does not divide `2(i - i')` for
  `0 < |i-i'| <= m'`, i.e. when `m' < q` (`q` odd prime). This holds at every regime `q` (R16-J2). Orthogonality gives (4.1);
  `N_{m'-i} = N_i` by conjugation. Correct.
* **Prop 4.2, re-derived.** `Pi` is rational; `sigma_s(m_i) = R_s m_i` with `(pi_s k)_x = k_{s^{-1} x}` a permutation of `A`
  (multiplication by `s in NQR` swaps `P_a` and `P'_{sa}`, by `s in QR` it maps `P_a -> P_{sa}`); `sigma_s` commutes with
  complex conjugation; so `N_i` is fixed by `Gal(Q(zeta)/Q)` and real, hence rational. Correct.
* **Exact independent check** (`referee_cap_checks/exact_Ni.py`, python-flint rationals; toy, not regime): `sdc_7`,
  `L = 2`, `m = 6 < 7`, `|A| = 757`, `dim V = 232`. All seven `N_i` are rational (the `zeta^1..zeta^5` coordinates
  vanish exactly), `N_i = N_{6-i}` exactly (e.g. `N_0 = 4936306540184820091642690987200/5253660409576722042453611`), and
  `4^{-6} sum binom(6,i)^2 N_i = 3467.99414530054766...` agrees with a 30-digit evaluation of `Psi_A(u_theta)` at
  `theta = 0, 0.3, 1.1` to `<= 3 * 10^{-28}` relative error.
* **Consequences: correct at the `kappa_V` level.** `kappa_V(hat u_theta) <= T*` forces every
  `N_i <= 4^{m'} binom(m',i)^{-2} d^{2m'} T*^2`, and one `N_i` above that bound forces `kappa_V(hat u_theta) > T*` on the
  whole circle; `ln N_0 > -166 545.5` reproduces. Per M1, "a cap needs all `N_i` tiny" must read "a `kappa_V`-cap
  needs ...". The "survival" bullet is correctly qualified as `kappa_V`-level survival, which is **not** survival of
  route 5 (that needs a lower bound on `kappa_H`).
* The cube comparison is legitimate: for each `i`, `mu(A) N_i^{(A)} <= N_i^{(T)}` (both are minimal `L^2` errors, the
  left one with the integrand cut to `A`), and `ln((m'!)^2 (2/3)^{m'}) = 26 447.1` reproduces.

## 7. Claim (6): zonotope degree obstruction (correct; scope of the reading)

* **Prop 6.1, re-derived.** The range of `k_0(P) + sum_{i in I} v_i(P)` over all `I` is `sum_i |v_i(P)| <= 2j`; summed over
  the `2q` patterns, `<= 4jq`; each nonzero integral `v` has `sum_P |v(P)| >= q + 1` (quoted sdc margin
  `sigma_lambda >= ||lambda||_1/2` with `lambda = (v, -sum v)`, `||lambda||_1 >= 2`; every integral `lambda` is relevant
  for `k = e_1`). So `m <= floor(4jq/(q+1)) = 1027` (`q = 3251`), `1367` (`q = 5479`). Reproduced.
* **Independent check of the margin input** (`margin_milp.py`, HiGHS, proved optimal, entries `|v_x| <= 3`, WLOG
  `v_0 >= 1` by translation and sign): `min sum_P |v(P)| = q + 1` exactly for `q = 7, 11, 19, 23, 31, 43`, attained by
  `e_x` and by pairs `e_x - e_y`. So the constant `q + 1` is sharp. This is evidence for bounded entries; the proof
  relies on the quoted, refereed margin.
* **Scope.** Prop 6.1 is about single zonotope certificates (products of `m` difference operators). Since
  `H_A(m') > 0` (Prop 3.5 gives `J(A)_{m'} != C[k]_{m'}`, so `V(A,m') != 0`), certificates of order `m'` exist; they are
  just not zonotopal. "The survival direction is blocked" (summary item 6) should read "zonotope certificates cannot
  reach order `m'`; no other certificate family with explicit values is known". The author's §6 text says this; the
  summary compresses it.
* **Remark 6.2 (minor m4).** Item 1 rests on an Alon-Furedi support bound quoted from memory `[M]`, so "proved negative
  statement" overstates; and the sub-box variant shows only that the Hoeffding-based *sufficient* condition fails
  (`64.25` against `1140`), not that the mechanism fails. Item 2: the threshold `3^{-3m'}/4 = e^{-6780.9}` against
  `mu(T \ A) >= e^{-33.03}` reproduces; the `p`-biased part is labelled heuristic, correctly.

## 8. Toy evidence (reproduces; one summary claim is false)

`referee_cap_checks/toy_ref.py` (own code, float64, pivoted QR; **evidence**; the cube `Psi_T` formula
`(m!)^2 [t^m] prod_x (1 + (2/3) u_x^2 t + (1/18) u_x^4 t^2)` was validated against a direct computation on `{-1,0,1}^5`).
My controls are 4 coordinate-permuted circles and 4 uniform `u in [-1,1]^q` with my own seed.

| box | `m` | `dim V` | `rho` on circles (8 points) | permuted circles | uniform `u` | author's table |
| --- | --- | --- | --- | --- | --- | --- |
| sdc7, `L=2`, `|A| = 757` | 3 | 721 | 0.4718 (const) | 0.250-0.470 | 0.114-0.365 | 0.472 |
| | 5 | 483 | 0.4023 (const) | 0.187-**0.405** | 0.039-0.369 | 0.402 |
| | 6 | 232 | 0.1955 (const) | 0.136-**0.247** | 0.015-**0.214** | 0.196 |
| | 7 `= q` | 44 | 0.0259-0.0339 (varies) | 0.033-0.074 | 0.0025-0.043 | 0.026-0.034 |
| | 8 | 1 | 7.17e-5 - 7.33e-5 | 9.6e-5 - 3.4e-3 | 5.4e-5 - 1.3e-3 | 7.0e-5 |
| sdc11, `L=2`, `|A| = 9373` | 3 | 9295 | 0.0328 (const) | 0.013-0.026 | 0.015-0.027 | 0.0328 |
| | 5 | 8140 | 0.0355 (const) | 0.0034-0.015 | 0.0053-0.021 | 0.0355 |
| | 6 | 5897 | 0.0289 (const) | 0.0015-0.011 | 0.0027-0.017 | 0.0289 |
| | 7 | 2162 | 0.0134 (const) | 0.0004-0.0050 | 0.0009-0.0093 | 0.0134 |

* `|A|`, every `dim V` and every circle value reproduce. The identity (4.1) holds to `<= 2.4 * 10^{-15}` relative error at
  every `m < q`, and `Psi_A` is constant on the circles to `<= 10^{-14}` for `m < q` (spread `0.23` at `m = q = 7`).
* **Minor error (m5).** The summary says "for `m < q`, the circle directions have the largest restriction ratio of all
  tested directions". This is false for `sdc7`: at `m = 6` the author's own table shows permuted circles up to `0.264` and
  uniform `u` up to `0.312` against the circle's `0.196`; my controls also beat the circle at `m = 5` and `m = 6`. It is
  true in the `sdc11` rows (author's and mine). Correct statement: "at `sdc11` (`m = 3..7`) the circles have the largest
  ratio among the tested directions; at `sdc7` they are typical". Also, at `sdc7`, `m = 8`, I find the uniform controls
  as low as `5.4 * 10^{-5}`, below the circle; "the circle is 10-100 times below the uniform controls" depends on the
  sample. None of this is regime evidence (the author says so).

## 9. kappa_V vs kappa_H, and evidence vs proof (summary)

* **M1** (§2 above) is the one substantive error: every "the cap needs" statement in cap.md is about approximants on
  `A`, i.e. `kappa_V`-caps. The test CIRCLE(q) is about `kappa_H` with `H' = e^{m'W/40}`, and a `kappa_H`-cap
  (through `TC_{H'} \ TC_V`) is neither excluded nor examined. §8.4 states the inequality backwards in effect
  ("`kappa_V` level, which is exactly what the cap needs").
* Survival is handled correctly at the `kappa_V` level and correctly flagged as insufficient for `kappa_H`; "height is
  not the obstacle" is a heuristic.
* Character sums: only the classical quadratic Gauss sum `g^2 = chi(-1) q` is used, and it is proved inline. No Weil-type
  bound is used. H_light is correctly labelled unproved and enters only through `tau_D`, which is handled by taking
  the `c_0`-free floor.
* Labels: proved items are proved (with M1 corrected); the toy table is labelled evidence; the "readings" after Prop 3.5
  and Lemma 3.1 are interpretations (m2, m3) and should be labelled as such in the summary too.

## 10. Corrected statement

At `C = 1`, `T = 1`, `sdc_q` (`q = 3 mod 4` prime, regime `q >= 3251`, `d = q`, `r = 2q`), barrier scale
`W = r tau`, `tau >= tau_D` of barrier Lemma 3.4 (so `tau >= 62.337` for every `c_0`), `T* = eps_0/(2 sqrt(m'+1)
e^{pi^2/(4q)})`:

1. **CIRCLE(q) is undecided** at every regime `q`. No cap and no lower bound. Route 5 stops (README §7.4). Growth rate
   (`M <= (2/3 + eps) log R/loglog R`) and the uniform bound are unchanged.
2. **Scale (proved arithmetic).** `ln T* <= -101 333` at `q = 3251` for every `c_0` (`-481 302` at `c_0 = 0.184`). A cap
   needs `kappa_H(hat u_theta) <= e^{-101 333}` for some `theta`; a cap by an approximant on `A` (via
   `kappa_H <= kappa_V`) needs `kappa_V(hat u_theta) <= e^{-101 333}`. cap.md analyses only the second kind; the
   `kappa_H`-only route (`TC_{H'} \ TC_V`) is unexamined.
3. **Gauss-sum identity (proved).** `(1 - g) hat k(1) = 2 sum_a zeta^a S_a`, `g^2 = -q`; `l_{u_theta} = sum_a c_a S_a` with
   `|c_a| <= 2/sqrt(q+1)`. (Exact checks `q <= 107` author, `q <= 139` referee.)
4. **Single relation families fail for approximants (proved).** Any degree-`< m'` approximant within `T* d^{m'}` on `A`
   has error `>= e^{-4507} d^{m'}` somewhere on `{-1,0,1}^q \ A` and `>= e^{-11 312} d^{m'}` somewhere on
   `W^{-1}({-j..j}^q) \ A` (`q = 3251`). Monomial rewriting of slab powers touches `<= e^{-2770}` of the `l^1`
   coefficient mass of `(sum c_a S_a)^{m'}` (proved; "no traction" is an interpretation).
5. **Fat-point form (proved).** `l_u^{m'} in (J_gen)_{m'}` iff every `G` of degree `m'` with multiplicity `>= m'-2` at the `q`
   coordinate points and `>= m'-2j` at the `2q` pattern points vanishes at `u`; `(J_gen)_{m'}` has relative dimension
   `<= e^{-367}` modulo `(k^3)`; `H_A(e) > 0` for `e <= q`.
6. **`L^2` splitting (proved; exact toy check).** For `m' < q`, `Psi_A(u_theta) = 4^{-m'} sum_i binom(m',i)^2 N_i`
   (constant on each circle), `N_i = N_{m'-i}` rational. `kappa_V`-level survival follows from one bound
   `ln N_0 > -166 545`; a `kappa_V`-cap needs every `N_i <= 4^{m'} binom(m',i)^{-2} d^{2m'} T*^2`. Neither bears directly
   on `kappa_H` survival.
7. **Zonotope bound (proved, from the refereed sdc margin).** Every zonotope certificate in `A` has order
   `<= floor(4jq/(q+1)) < 4j <= m'/2`; certificates of order `m'` exist but none with explicit value is known.
8. **Toys (evidence).** On `sdc11` (`L = 2`, `m = 3..7`) the circles have the largest restriction ratio of the tested
   directions; on `sdc7` they are typical (beaten at `m = 5, 6`). No near-zero appears; no regime information.
