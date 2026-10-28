# Round 17 — calibration of the circle test CIRCLE(q) on small sdc_q (EVIDENCE ONLY)

Status: final for this task. Unrefereed. Everything in §4-§7 is evidence about toy boxes unless marked proved/exact.

## 0. Summary

**Status: calibration evidence, plus one small proved reduction. No growth improvement, no uniform bound; CIRCLE(q)
is not decided.** `M <= (2/3 + eps) log R/loglog R` stands.

**Proved (Prop 2.1).** For any `AGL(1,q)`-invariant box and every `m < q`, `Psi_A` (the squared `L^2(mu_A)` distance of
`l_u^m` from `Poly_{<m}`) is **exactly constant** on the Fourier circle `C_1`:
`Psi_A(u_theta) = 4^{-m} sum_i binom(m,i)^2 ||R_i||^2`, where `R_i` is the residual of `khat_1^i khat_{-1}^{m-i}` inside one
translation eigenspace, and every `||R_i||^2` is **rational**. Hence `inf_theta kappa_V(hat u_theta) >= Psi_A(C_1)^{1/2}/d^m`
rigorously. This turns the circle test into one number per `(q, L, m)` and makes exact calibration cheap. (Lemma 7.1, also
proved: `rho(A) >= d + 1` for every half-cube, so regime orders satisfy `m'/rho(A) < 0.633`.)

**Computed** (exact multimodular for `L = 1`, all primes `7 <= q <= 103`, and `L = 2`, `q = 7, 11` (all `m`) and `q = 19, 23` (`m <= 4`); complex128, validated
against the exact values to `>= 10` digits, for `L = 2`, `q = 19` (`m = 5, 6, 7`) and `q = 23` (`m = 5, 6`);
LP for `kappa_V`; LLL for `kappa_H`). `q = 31, 43` at `L = 2` were not reachable (enumeration of `A`); `q = 31, 43`
are covered at `L = 1`.
* The circle is **never** in `TC_V` below saturation, and is not even close: the relative residual
  `rel = Psi_A(C_1)/E l^{2m}` decays like the **uncut cube's** on the same circle (`rel_A/rel_cube in [0.29, 1.31]` for
  `m <= rho - 2`), i.e. at `0.2 - 0.86` per degree, with a single faster drop at the last degree before saturation
  (`sdc11`, `m = 8`, `dim V = 5`: `rel = 6.0e-6`).
* The circle is not more tangent than random directions (above the random range for `sdc11`, `sdc19`, and `sdc23` except `m = 5`; inside it for
  `sdc7`; at most `12%` below it for the pair boxes `L = 1`).
* `kappa_V(hat u_theta)` is constant on the circle to within `0.6%` (`sdc7`, `L = 2`) -- far better than barrier
  Thm 3.2's `(m+1)^{-1/2}`.
* `kappa_H`: `V_H = V(A,m)` already for `H <= 4` in every computed case with `m < q`; no `H`-gap is visible.
* Against the threshold: each toy circle meets `eps_0/(2 sqrt(m+1) e^{pi^2/(4q)})` only at `tau* <= 2.22`, while the
  barrier needs `tau >= tau_D = 164 - 221`.
* **Scaling shape (evidence).** In every toy below saturation,
  `inf_theta kappa_V(hat u_theta) ~ k_lo = (rel_A(m) E_A l^{2m})^{1/2}/q^m ~ (rel_A(m) (2m-1)!!)^{1/2} (sigma_A/q)^m`, with
  `sigma_A^2 = E_A l^2 = Psi_A(C_1)|_{m=1}` (`2.19 - 3.35` for `L = 2`, `~ 1` for `L = 1`) and `rel_A(m) ~ rel_cube(m) ~ e^{-(0.2 .. 0.6) m}`
  (e.g. `sdc19`, `L = 2`, `m = 6`: model `1.6e-5`, `k_lo = 1.25e-5`). So `ln inf kappa` is **linear in `m`** with slope
  `-(ln(q/sigma_A) + O(1))`, nearly independent of `L` beyond `sigma_A`, and converges in `q` at fixed `(L, m)` after the
  `q^{-m}` normalisation; the threshold shape `eps_0/(2 sqrt(m+1))` contributes only `O(ln m)`, so the comparison is
  decided entirely by `eps_0 ~ e^{-W/4} = e^{-q tau/2}`, which no toy approaches (`tau* <= 2.22`).
* All exact circle data of round 16 are reproduced (`sdc7`: `Phi_c = 7717.5`; `sdc11`: `dim V(A,8) = 5`, only the
  trivial eigenspace; `sdc19`: `m = 5, 6` deficiencies) and made quantitative (§6).

**Extrapolation (labelled evidence, §7).** A "yes" through `kappa_V` at `q = 3251` needs `ln Psi_A(C_1) <= -9.3e5`,
i.e. a relative residual decaying by `>= 466` per degree (`>= 9321` at `q ~ 10^6`), against `<= 1.5` in every toy and
`~ ln 2` for the cube. If regime circles behave like the toys at comparable `m/rho`, CIRCLE(q) is a robust **no** by a
factor `~ e^{4.8e5}` in `kappa` at the first regime prime, and the gap grows with `q`. This does **not** decide CIRCLE:
the cap could still come from an `H`-gap (`kappa_H << kappa_V`) or from an exact regime identity with no toy shadow;
and a proved "no" needs a lower bound for one `||R_i||` at regime size plus (P1)/(P2). Read with the judge's rule:
**route 5 is undecided after this calibration**; the evidence points to "no" (route 5 survives the circle test) but
supplies no certificate.


## 1. Objects and conventions (recalled from round16)

Labels: *proved* (all quantifiers, here); *exact* (integer / rational / mod-`p` computation, stated how); *evidence*
(floating point, or extrapolation). **Everything in §4-§7 is calibration evidence about toy boxes; no toy is a regime
box** (the first regime `sdc_q` has `q = 3251`, R16-J2).

* `sdc_q` (`q` prime, `q = 3 mod 4`): coordinates `x in F_q` (base row `inf` omitted, `d = q`), the `2q` patterns
  `P_a = {a} u {x : chi(a-x) = 1}`, `P'_a = {a} u {x : chi(a-x) = -1}` (`chi` = Legendre symbol), as in
  `round16/barrier_checks/flib.py` and `referee_barrier_checks/sdc_sym.py`.
* Toy slab box at half-width `L` (the regime uses `L = j = floor(w/2)`):
  `A = A_L(q) = {k in {-1,0,1}^q : |k(P)| <= L for all 2q patterns}`, `U = |A|`, `mu_A` uniform.
  `A` is invariant under `AGL(1,q)` acting on coordinates (referee_barrier §4, proved there for all `q`).
* `khat_b(k) = sum_x zeta^{bx} k_x` (`zeta = e^{2 pi i/q}`); translation `(gk)_t = k_{t-1}` gives `khat_b(gk) = zeta^b khat_b(k)`.
* Circle `C_1`: `u_theta(x) = cos(theta + 2 pi x/q)`, so `<k, u_theta> = Re(e^{i theta} khat_1(k))`. `hat u = u/|u|_inf`.
* For `m >= 1`: `Psi_A(u) = min_{P in Poly_{<m}} E_{mu_A}(l_u^m - P)^2` (probability normalisation),
  `Cheb_A(u) = dist_inf(l_u^m|_A, Poly_{<m}|_A) = d^m kappa_V(u)` (nearstrata Prop 2.1), and (reduction Thm 2.1(2))

  ```text
  Psi_A(u)^{1/2}  <=  Cheb_A(u)  <=  (U Psi_A(u))^{1/2}.                                            (1.1)
  ```
* `V(A,m)` = measures on `A` annihilating `Poly_{<m}|_A`; `dim V(A,m) = U - dim Poly_{<m}|_A`;
  `h_A(e) = dim Poly_{<=e}|_A - dim Poly_{<e}|_A` (Hilbert function of `R[k]/J(A)`); `rho(A)` = least `m` with
  `dim V(A,m) = 0`.
* `kappa_H` uses `V_H` = real span of integer `c in V(A,m)` with `|c|_inf <= H`; `kappa_H <= kappa_V`, and
  `d^m kappa_H(u) >= Psi_{A,H}(u)^{1/2}` with `Psi_{A,H} = Psi_A` whenever `V_H = V(A,m)` (reduction §4).
* Cap threshold of CIRCLE(q): `inf_theta kappa_H(hat u_theta) <= eps_0/(2 sqrt(m'+1) e^{pi^2/(4q)})`,
  `eps_0 = C d e^{-W/4}/(m'+1)`, `W = 2q tau`, `tau = K ln n >= tau_D` (barrier Lemma 3.4).

**Toy boxes used.** `L = 1`: `A_1(q) = {0} u {+-e_x} u {e_x - e_y : x != y}` exactly (`U = 1 + q(q+1)`; enumerated for all
`q <= 107` below; it is the "pair box"). `L = 2`: `U = 757, 9373, 264291, 818755` for `q = 7, 11, 19, 23` (the first
three agree with the record's counts). Every `k in A_L` has `|supp k| <= 4L^2 + (2L-1)^2/q` (proof: with `S = sum k`,
`k(P_a) - k(P'_a) = (chi * k)(a)` and `k(P_a) + k(P'_a) = S + k_a`, so `|(chi*k)(a)| <= 2L`; Parseval and `|chi-hat(b)|^2 = q`
for `b != 0` give `sum_a (chi*k)(a)^2 = q|supp k| - S^2`; and `|S| (q+1)/2 = |sum_a k(P_a)| <= qL`). This is used only to
prune the enumeration (`enumA.c`). The `L = 2` enumeration for `q = 31, 43` was not completed (see §3).


## 2. Exact reduction to the `Z_q`-isotypic components (proved)

**Proposition 2.1.** Let `q` be an odd prime, `A` any translation-invariant subset of `{-1,0,1}^{F_q}` that is also
invariant under all multiplications `x -> tx`, `t in F_q^*` (e.g. `A_L(q)` of `sdc_q`), and `1 <= m < q`. Put
`F_i = khat_1^i khat_{-1}^{m-i}` (`0 <= i <= m`) and let `R_i` be the `L^2(mu_A)`-residual of `F_i|_A` modulo `Poly_{<m}|_A (x) C`.
1. `Psi_A(u_theta) = 4^{-m} sum_{i=0}^m binom(m,i)^2 ||R_i||^2` and `E_{mu_A} l_{u_theta}^{2m} = 4^{-m} sum_i binom(m,i)^2 ||F_i||^2`
   for **every** `theta`: both are constant on `C_1`.
2. `||R_i|| = ||R_{m-i}||`, and every `||R_i||^2` and `||F_i||^2` is a **rational number**; so `Psi_A(C_1) in Q`.
3. `R_i` lies in the translation eigenspace of character `zeta^{2i-m}`, and `R_i` is computable inside that eigenspace
   (functions on the free translation orbits, plus the fixed point `0` for the trivial character).
4. For the normalised direction: `Psi_A(hat u_theta) = Psi_A(C_1)/|u_theta|_inf^{2m}` with `cos(pi/(2q)) <= |u_theta|_inf <= 1`,
   so `inf_theta Psi_A(hat u_theta) = Psi_A(C_1)` (attained at `theta = 0`), and by (1.1)

   ```text
   inf_theta kappa_V(hat u_theta)  >=  Psi_A(C_1)^{1/2} / d^m .                                      (2.1)
   ```

*Proof.* (1) `<k,u_theta> = (e^{i theta} khat_1 + e^{-i theta} khat_{-1})/2`, so `l_{u_theta}^m = 2^{-m} sum_i binom(m,i) e^{i theta(2i-m)} F_i`.
`F_i(gk) = (zeta khat_1)^i (zeta^{-1} khat_{-1})^{m-i} = zeta^{2i-m} F_i(k)`. The orthogonal projection `Pi` onto `Poly_{<m}|_A`
commutes with the unitary operator `f -> f o g` (which preserves `A`, `mu_A` and `Poly_{<m}`), so `R_i = F_i - Pi F_i` is an
eigenvector with eigenvalue `zeta^{2i-m}`. For `0 <= i < i' <= m` the exponents differ by `2(i'-i) in (0, 2q)`, not divisible by
the odd prime `q`; eigenvectors of a unitary operator for distinct eigenvalues are orthogonal. Hence the residual of
`l^m_{u_theta}` is `2^{-m} sum_i binom(m,i) e^{i theta(2i-m)} R_i` with pairwise orthogonal summands of `theta`-free norms; the same
for `F_i`. (2) `F_{m-i} = conj(F_i)` and `Pi` is real, so `R_{m-i} = conj(R_i)`. `Pi` is a rational matrix and `F_i(k) in Z[zeta]`,
so `||R_i||^2 = U^{-1} sum_k R_i(k) conj(R_i(k)) in Q(zeta) cap R`. For `t in F_q^*` let `sigma_t(zeta) = zeta^t` and
`(pi_t k)_y = k_{t^{-1} y}`; then `sigma_t(khat_{+-1}(k)) = khat_{+-1}(pi_t k)`, so `sigma_t(F_i) = F_i o pi_t`. `pi_t` preserves `A`,
hence commutes with `Pi`, so `sigma_t(R_i) = R_i o pi_t`; conjugation is `sigma_{-1}` and commutes with `sigma_t`; summing over
`k in A = pi_t A` gives `sigma_t(||R_i||^2) = ||R_i||^2`. A Galois-fixed element of `Q(zeta)` is rational. (3) is the
eigenvector statement of (1). (4) `l_{hat u}^m = l_u^m/|u|_inf^m`. The `2q` points `theta + 2 pi x/q` and `theta + pi + 2 pi x/q`
are equally spaced at distance `pi/q`, so one of them is within `pi/(2q)` of `0`, giving `|u_theta|_inf >= cos(pi/(2q))`;
`|u_0|_inf = 1`. Then (1.1) and `Cheb_A = d^m kappa_V`. QED

*Remarks.* (a) This sharpens barrier Thm 3.2 on `C_1`: there `inf >= (m+1)^{-1/2} sup` for `kappa`; here `Psi_A` (the `L^2`
form) is **exactly constant** on `C_1`, for every `m < q` (the referee's `m < N` relaxation is exactly what is needed). It
does not make `kappa_V` constant (Chebyshev is not an `L^2` quantity; the LP values of §4.1 and §5.1 show it varies, though by `<= 15%`). (b) By the `AGL(1,q)`
invariance the same values hold on every circle `C_a`. (c) Prop 2.1 is what makes the calibration cheap: one number per
`(q, L, m)`, computed in `(m+1)/2` eigenspaces of dimension `~U/q`. (d) (2.1) is the only rigorous `kappa` statement the
computation gives directly; a lower bound for `kappa_H` additionally needs `V_H = V(A,m)` (§5.6).


## 3. Computation method and certification

All scripts are in `calib_checks/` (outputs in `calib_checks/out/`). Nothing below is used in a proof except Prop 2.1.

1. **Enumeration** (`enumA.c`, exact integer DFS with slab and support pruning). Writes one representative per
   translation orbit. Counts reproduce the record (`757, 9373, 264291` for `L = 2`, `q = 7, 11, 19`).
2. **Exact isotypic engine** (`exact_iso.py`). For primes `p = 1 mod q` (`zeta -> z in F_p`, complex conjugation
   `-> z^{-1}`), the eigenspace of character `zeta^c` is spanned by `P_c(k^alpha) = sum_s zeta^{-cs} (g^s k)^alpha` over
   monomial orbit representatives `alpha in {0,1,2}^q`, `|alpha| < m` (support `<= max supp A`). An independent subset `S`
   (pivots of the rref mod `p`) is fixed at the first prime; `||R_i||^2 = U^{-1}(F*WF - F*WM_S G^{-1} M_S*WF)`,
   `G = M_S*WM_S` (Hermitian Gram, weights `W = q` per free orbit, `1` for `k = 0`). `det G != 0 mod p` at every prime
   certifies that `S` is independent over `Q(zeta)`. The rationals `||R_i||^2, ||F_i||^2` (Prop 2.1(2)) are recovered by CRT
   and rational reconstruction, accepted when stable over 3 further primes. **Status: exact modulo the standard
   multimodular caveat** (the char-0 rank of the full column set equals the mod-`p` rank; checked equal at the first two
   primes, and rank mismatches at any later prime are reported -- none occurred). Totals `dim Poly_{<m}|_A` (all `q`
   eigenspaces, one prime) give `dim V(A,m)` and `h_A`. For the large boxes (`sdc19`, `sdc23`, `L = 2`) the same engine
   is run with smaller primes (`exact_iso_big.py`: `2^23 < p < 2^24`; `exact_iso_big23.py`: `2^22 < p < 2^23`), so that
   every `int64` dot product (length `<= R + 1` orbits, entries `< p^2`) is exact; the guard `(R+1) p^2 < 2^63` is asserted.
3. **Independent check without symmetry** (`full_check.py`): full monomial matrix on all of `A`, pivoted QR, floating
   point, at five `theta`. Agrees with 2 to 10-12 digits on `sdc7` (`L = 1`, `m = 2, 3`), and confirms the
   `theta`-independence of `Psi_A` on `C_1` numerically.
4. **Floating isotypic engine** (`float_iso.py`, complex128): same eigenspace columns, degree-block CGS2 + SVD rank
   truncation at absolute tolerance `10^{-9}` on unit columns. Also computes `Psi_A` at `6` random `u` (`|u|_inf = 1`,
   uniform in `[-1,1]^q`), by splitting `l_u^m` into all eigenspaces. Validated against 2: identical dimensions and
   `Psi_A(C_1)` to `>= 10` digits on `sdc7` (`L = 2`, `m <= 6`) and `sdc11` (`L = 2`, see §4). **Evidence.**
5. **Cube comparison** (`cube_psi.py`, 50 digits): `Psi_cube(u) = (m!)^2 [t^m] prod_x (1 + (2/3)u_x^2 t + u_x^4 t^2/18)`
   (reduction Thm 5.1(a)) and `E_cube l_u^{2m}`; constant on `C_1` (checked at two `theta`).
6. **Chebyshev / `kappa_V`** (`cheb_lp.py`): LP (HiGHS) on all of `A` over the `theta`-grid of the fundamental domain
   `[0, pi/(2q)]` (the symmetries `theta -> -theta`, `theta -> theta + 2 pi/q` and `u_{theta + pi/q} = -u_theta o (translation)`
   are box automorphisms). **Evidence** (floating LP). The rigorous statement available is (2.1).
7. **`kappa_H`** (`kappa_h.py`): integer basis of `V(A,m)` (exact `fmpz` nullspace) and LLL (flint); if every reduced basis
   vector has `|c|_inf <= H_red` then `V_H = V(A,m)` for `H >= H_red`, hence `kappa_H = kappa_V` there (exact for the
   stated `H`). Feasible only for small `dim V`.

**Not done / infeasible.** `L = 2` enumeration for `q = 31, 43` was not completed (the pruned DFS is too slow; `|A|`
is expected in the `10^7 - 10^8` range); `q = 23, L = 2` is covered only for small `m` (eigenspace dimension `35598`).
`L = 1` was run for every prime `q = 3 mod 4` with `7 <= q <= 103`.


## 4. Data

### 4.1 Main table (every `(q, L, m)` with `m < q` and `dim V(A,m) >= 2`)

Columns: `U = |A|`; `dim V = dim V(A,m)` (exact mod `p` for `E` rows; floating rank for `F` rows, which agree with the
record's mod-`p` deficiencies where those exist); `Psi_A(C_1)` (probability normalisation, constant on `C_1`; `E` = exact
rational, decimal shown, full rationals in `out/exact_*.txt`; `F` = complex128, agreeing with `E` to `>= 10` digits on
every overlap); `rel = Psi_A(C_1)/E_{mu_A} l^{2m}` (also constant on `C_1`); `Psi_A/Psi_cube` on `C_1`; `rel` of the uncut
cube on `C_1`; the range of `rel` over 6 random `u` (`|u|_inf = 1`); `k_lo = Psi_A(C_1)^{1/2}/d^m`, a **rigorous** lower
bound for `inf_theta kappa_V(hat u_theta)` whenever `Psi_A(C_1)` is exact (2.1); the LP value of `inf_theta kappa_V(hat u_theta)`
over the grid with `inf/sup` in brackets; `tau*` (§5.7) from `k_lo` and from the LP value; `tau_D(q)`.

<!--TABLE-->
| q | L | U | m | dim V | Psi_A(C_1) | rel | -ln(rel)/m | Psi_A/Psi_cube | cube rel | random-u rel range | k_lo = sqrt(Psi)/d^m | LP inf kappa_V (inf/sup) | tau* (k_lo / LP) | tau_D |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 7 | 1 | 57 | 1 | 56 | 0.982456 (E) | 1 | 0.000 | 0.421 | 1 | [1, 1] | 1.416e-01 | - | 0.52 | 164 |
| 7 | 1 | 57 | 2 | 49 | 1.0611 (E) | 0.5237 | 0.323 | 0.116 | 0.627 | [0.508, 0.64] | 2.102e-02 | 3.687e-02 (0.903) | 0.89 / 0.73 | 164 |
| 7 | 1 | 57 | 3 | 21 | 1.20888 (E) | 0.2218 | 0.502 | 0.0273 | 0.327 | [0.194, 0.238] | 3.206e-03 | 4.993e-03 (0.856) | 1.30 / 1.18 | 164 |
| 11 | 1 | 133 | 1 | 132 | 0.992481 (E) | 1 | 0.000 | 0.271 | 1 | [1, 1] | 9.057e-02 | - | 0.52 | 174 |
| 11 | 1 | 133 | 2 | 121 | 1.124 (E) | 0.533 | 0.315 | 0.0466 | 0.642 | [0.509, 0.593] | 8.762e-03 | 1.587e-02 (0.960) | 0.83 / 0.72 | 174 |
| 11 | 1 | 133 | 3 | 55 | 1.19408 (E) | 0.2081 | 0.523 | 0.00561 | 0.356 | [0.153, 0.285] | 8.210e-04 | 1.413e-03 (0.940) | 1.18 / 1.08 | 174 |
| 19 | 1 | 381 | 1 | 380 | 0.997375 (E) | 1 | 0.000 | 0.157 | 1 | [1, 1] | 5.256e-02 | - | 0.42 | 185 |
| 19 | 1 | 381 | 2 | 361 | 1.17453 (E) | 0.5414 | 0.307 | 0.0156 | 0.653 | [0.536, 0.655] | 3.002e-03 | 5.465e-03 (0.986) | 0.66 / 0.60 | 185 |
| 19 | 1 | 381 | 3 | 171 | 1.17254 (E) | 0.197 | 0.542 | 0.000926 | 0.375 | [0.182, 0.28] | 1.579e-04 | 2.856e-04 (0.980) | 0.93 / 0.86 | 185 |
| 23 | 1 | 553 | 1 | 552 | 0.998192 (E) | 1 | 0.000 | 0.13 | 1 | [1, 1] | 4.344e-02 | - | 0.39 | 189 |
| 23 | 1 | 553 | 2 | 529 | 1.18716 (E) | 0.5437 | 0.305 | 0.0106 | 0.655 | [0.554, 0.664] | 2.060e-03 | 3.746e-03 (0.991) | 0.60 / 0.55 | 189 |
| 23 | 1 | 553 | 3 | 253 | 1.16586 (E) | 0.1942 | 0.546 | 0.000502 | 0.38 | [0.198, 0.274] | 8.874e-05 | 1.621e-04 (0.986) | 0.83 / 0.78 | 189 |
| 31 | 1 | 993 | 1 | 992 | 0.998993 (E) | 1 | 0.000 | 0.0967 | 1 | [1, 1] | 3.224e-02 | - | 0.33 | 195 |
| 31 | 1 | 993 | 2 | 961 | 1.20292 (E) | 0.5466 | 0.302 | 0.00584 | 0.658 | [0.552, 0.61] | 1.141e-03 | - | 0.50 | 195 |
| 31 | 1 | 993 | 3 | 465 | 1.15679 (E) | 0.1906 | 0.552 | 0.000195 | 0.385 | [0.201, 0.255] | 3.610e-05 | - | 0.70 | 195 |
| 43 | 1 | 1893 | 1 | 1892 | 0.999472 (E) | 1 | 0.000 | 0.0697 | 1 | [1, 1] | 2.325e-02 | - | 0.27 | 202 |
| 43 | 1 | 1893 | 2 | 1849 | 1.21579 (E) | 0.549 | 0.300 | 0.00304 | 0.661 | [0.568, 0.607] | 5.963e-04 | - | 0.41 | 202 |
| 43 | 1 | 1893 | 3 | 903 | 1.1488 (E) | 0.1877 | 0.558 | 7.04e-05 | 0.389 | [0.213, 0.244] | 1.348e-05 | - | 0.56 | 202 |
| 47 | 1 | 2257 | 1 | 2256 | 0.999557 (E) | 1 | 0.000 | 0.0638 | 1 | - | 2.127e-02 | - | 0.25 | 204 |
| 47 | 1 | 2257 | 2 | 2209 | 1.21865 (E) | 0.5495 | 0.299 | 0.00254 | 0.661 | - | 4.997e-04 | - | 0.39 | 204 |
| 47 | 1 | 2257 | 3 | 1081 | 1.14695 (E) | 0.1871 | 0.559 | 5.35e-05 | 0.39 | - | 1.032e-05 | - | 0.53 | 204 |
| 59 | 1 | 3541 | 1 | 3540 | 0.999718 (E) | 1 | 0.000 | 0.0508 | 1 | - | 1.695e-02 | - | 0.22 | 209 |
| 59 | 1 | 3541 | 2 | 3481 | 1.22494 (E) | 0.5507 | 0.298 | 0.00161 | 0.662 | - | 3.179e-04 | - | 0.33 | 209 |
| 59 | 1 | 3541 | 3 | 1711 | 1.1428 (E) | 0.1857 | 0.561 | 2.65e-05 | 0.392 | - | 5.205e-06 | - | 0.46 | 209 |
| 67 | 1 | 4557 | 1 | 4556 | 0.999781 (E) | 1 | 0.000 | 0.0448 | 1 | - | 1.492e-02 | - | 0.20 | 212 |
| 67 | 1 | 4557 | 2 | 4489 | 1.22789 (E) | 0.5513 | 0.298 | 0.00125 | 0.663 | - | 2.468e-04 | - | 0.30 | 212 |
| 67 | 1 | 4557 | 3 | 2211 | 1.14081 (E) | 0.185 | 0.562 | 1.8e-05 | 0.393 | - | 3.551e-06 | - | 0.42 | 212 |
| 71 | 1 | 5113 | 1 | 5112 | 0.999804 (E) | 1 | 0.000 | 0.0422 | 1 | - | 1.408e-02 | - | 0.19 | 213 |
| 71 | 1 | 5113 | 2 | 5041 | 1.22912 (E) | 0.5515 | 0.298 | 0.00111 | 0.663 | - | 2.199e-04 | - | 0.29 | 213 |
| 71 | 1 | 5113 | 3 | 2485 | 1.13997 (E) | 0.1847 | 0.563 | 1.5e-05 | 0.394 | - | 2.983e-06 | - | 0.40 | 213 |
| 79 | 1 | 6321 | 1 | 6320 | 0.999842 (E) | 1 | 0.000 | 0.038 | 1 | - | 1.266e-02 | - | 0.18 | 215 |
| 79 | 1 | 6321 | 2 | 6241 | 1.23121 (E) | 0.5519 | 0.297 | 0.000901 | 0.663 | - | 1.778e-04 | - | 0.27 | 215 |
| 79 | 1 | 6321 | 3 | 3081 | 1.13853 (E) | 0.1843 | 0.564 | 1.08e-05 | 0.394 | - | 2.164e-06 | - | 0.37 | 215 |
| 83 | 1 | 6973 | 1 | 6972 | 0.999857 (E) | 1 | 0.000 | 0.0361 | 1 | - | 1.205e-02 | - | 0.17 | 217 |
| 83 | 1 | 6973 | 2 | 6889 | 1.23211 (E) | 0.5521 | 0.297 | 0.000816 | 0.664 | - | 1.611e-04 | - | 0.26 | 217 |
| 83 | 1 | 6973 | 3 | 3403 | 1.13791 (E) | 0.1841 | 0.564 | 9.33e-06 | 0.395 | - | 1.866e-06 | - | 0.36 | 217 |
| 103 | 1 | 10713 | 1 | 10712 | 0.999907 (E) | 1 | 0.000 | 0.0291 | 1 | - | 9.708e-03 | - | 0.15 | 221 |
| 103 | 1 | 10713 | 2 | 10609 | 1.23555 (E) | 0.5527 | 0.296 | 0.00053 | 0.664 | - | 1.048e-04 | - | 0.22 | 221 |
| 103 | 1 | 10713 | 3 | 5253 | 1.1355 (E) | 0.1833 | 0.566 | 4.83e-06 | 0.396 | - | 9.752e-07 | - | 0.30 | 221 |
| 7 | 2 | 757 | 1 | 756 | 2.19155 (E) | 1 | 0.000 | 0.939 | 1 | [1, 1] | 2.115e-01 | - | 0.40 | 164 |
| 7 | 2 | 757 | 2 | 749 | 10.0733 (E) | 0.6771 | 0.195 | 1.1 | 0.627 | [0.545, 0.652] | 6.477e-02 | 2.061e-01 (1.000) | 0.57 / 0.24 | 164 |
| 7 | 2 | 757 | 3 | 721 | 60.4324 (E) | 0.3731 | 0.329 | 1.36 | 0.327 | [0.239, 0.352] | 2.266e-02 | 6.595e-02 (0.997) | 0.75 / 0.44 | 164 |
| 7 | 2 | 757 | 4 | 644 | 348.048 (E) | 0.1545 | 0.467 | 1.49 | 0.148 | [0.0855, 0.171] | 7.770e-03 | 2.091e-02 (0.998) | 0.96 / 0.67 | 164 |
| 7 | 2 | 757 | 5 | 483 | 1427.15 (E) | 0.04015 | 0.643 | 1.16 | 0.0579 | [0.0238, 0.0653] | 2.248e-03 | 5.359e-03 (0.996) | 1.23 / 0.98 | 164 |
| 7 | 2 | 757 | 6 | 232 | 3467.99 (E) | 0.005779 | 0.859 | 0.565 | 0.0196 | [0.00482, 0.0179] | 5.006e-04 | 1.073e-03 (0.995) | 1.59 / 1.38 | 164 |
| 11 | 2 | 9373 | 1 | 9372 | 2.76848 (E) | 1 | 0.000 | 0.755 | 1 | [1, 1] | 1.513e-01 | - | 0.42 | 174 |
| 11 | 2 | 9373 | 2 | 9361 | 15.9589 (E) | 0.6756 | 0.196 | 0.661 | 0.642 | [0.601, 0.653] | 3.302e-02 | 1.654e-01 (0.973) | 0.59 / 0.30 | 174 |
| 11 | 2 | 9373 | 3 | 9295 | 131.81 (E) | 0.3954 | 0.309 | 0.62 | 0.356 | [0.308, 0.355] | 8.626e-03 | 4.756e-02 (0.960) | 0.76 / 0.44 | 174 |
| 11 | 2 | 9373 | 4 | 9020 | 1416.47 (E) | 0.2182 | 0.381 | 0.64 | 0.178 | [0.144, 0.172] | 2.571e-03 | 1.365e-02 (0.947) | 0.91 / 0.61 | 174 |
| 11 | 2 | 9373 | 5 | 8140 | 17047.8 (E) | 0.1072 | 0.447 | 0.671 | 0.0822 | [0.058, 0.075] | 8.107e-04 | <= 3.516e-03 (theta=0 only) | 1.07 / 0.81 | 174 |
| 11 | 2 | 9373 | 6 | 5897 | 167404 (E) | 0.03658 | 0.551 | 0.546 | 0.0349 | [0.0186, 0.0257] | 2.310e-04 | - | 1.26 | 174 |
| 11 | 2 | 9373 | 7 | 2162 | 949954 (E) | 0.006483 | 0.720 | 0.253 | 0.0137 | [0.00352, 0.00487] | 5.002e-05 | - | 1.50 | 174 |
| 11 | 2 | 9373 | 8 | 5 | 30136.9 (E) | 6e-06 | 1.503 | 0.000662 | 0.0049 | [2.73e-07, 2.48e-06] | 8.099e-07 | - | 2.22 | 174 |
| 19 | 2 | 264291 | 1 | 264290 | 3.25808 (E) | 1 | 0.000 | 0.514 | 1 | [1, 1] | 9.500e-02 | - | 0.36 | 185 |
| 19 | 2 | 264291 | 2 | 264271 | 20.5063 (E) | 0.6589 | 0.209 | 0.272 | 0.653 | [0.636, 0.65] | 1.254e-02 | - | 0.51 | 185 |
| 19 | 2 | 264291 | 3 | 264081 | 183.089 (E) | 0.3811 | 0.322 | 0.145 | 0.375 | [0.351, 0.369] | 1.973e-03 | - | 0.66 | 185 |
| 19 | 2 | 264291 | 4 | 262770 | 2106.12 (E) | 0.2089 | 0.392 | 0.0793 | 0.2 | [0.175, 0.193] | 3.521e-04 | - | 0.81 | 185 |
| 19 | 2 | 264291 | 5 | 255816 | 27393.3 (F) | 0.1039 | 0.453 | 0.0421 | 0.101 | [0.0805, 0.0938] | 6.684e-05 | - | 0.95 | 185 |
| 19 | 2 | 264291 | 6 | 225835 | 347191 (F) | 0.04294 | 0.525 | 0.0194 | 0.0488 | [0.0315, 0.04] | 1.252e-05 | - | 1.10 | 185 |
| 19 | 2 | 264291 | 7 | 133989 | 3.45409e+06 (F) | 0.0124 | 0.627 | 0.00648 | 0.0225 | [0.00996, 0.0107] | 2.079e-06 | - | 1.27 | 185 |
| 23 | 2 | 818755 | 1 | 818754 | 3.35384 (E) | 1 | 0.000 | 0.437 | 1 | [1, 1] | 7.962e-02 | - | 0.33 | 189 |
| 23 | 2 | 818755 | 2 | 818731 | 21.6665 (E) | 0.6583 | 0.209 | 0.194 | 0.655 | [0.642, 0.656] | 8.799e-03 | - | 0.47 | 189 |
| 23 | 2 | 818755 | 3 | 818455 | 199.177 (E) | 0.3814 | 0.321 | 0.0857 | 0.38 | [0.359, 0.378] | 1.160e-03 | - | 0.61 | 189 |
| 23 | 2 | 818755 | 4 | 816178 | 2336.06 (E) | 0.2075 | 0.393 | 0.0383 | 0.205 | [0.185, 0.204] | 1.727e-04 | - | 0.75 | 189 |
| 23 | 2 | 818755 | 5 | 801757 | 31253.9 (F) | 0.1035 | 0.454 | 0.0165 | 0.106 | [0.0882, 0.106] | 2.747e-05 | - | 0.88 | 189 |
| 23 | 2 | 818755 | 6 | 727445 | 414429 (F) | 0.0436 | 0.522 | 0.00618 | 0.0523 | [0.0373, 0.0423] | 4.349e-06 | - | 1.02 | 189 |
<!--/TABLE-->

`L = 1`: `rho(A_1) = 4` for every `q` (`dim V(A_1, 4) = 0`), so only `m <= 3` qualify. `L = 2`: `rho = 9` (`q = 7, 11`);
for `q = 19` the runs reached `m = 7` (exact for `m <= 4`, `F` for `m = 5, 6, 7`; `dim V(A,7) = 133989`, i.e. `51%` of `U`, still far from saturation; per-eigenspace `dim Poly_{<7,c} = 6858` for `c = 0` and `c = 1` (out of `13911`, `13910`); `max |supp k| = 8` on `A_2(19)`).
`q = 23`, `L = 2`: exact for `m <= 4` (`out/exact_q23_L2.txt`, `22` primes at `m = 4`), `F` at `m = 5, 6` (agreeing with the exact values to `10` digits for `m <= 4`); `dim V(A,6) = 727445` (`89%` of `U`), per-eigenspace `dim Poly_{<6,c} = 3970` (`c = 0, 1`). Two concurrent runs (`q = 19`, `m = 7`; `q = 23`, `m = 6`) were killed by the container's memory limit (kernel memcg OOM, `~4.5 GB` each); `q = 19`, `m = 7` was then rerun alone (`out/float_q19_L2_m7run4.txt`, `682 s`, with 2 random `u` instead of 6). `q = 23`, `m = 6` was likewise rerun alone (`out/float_q23_L2_m6run4.txt`, `889 s`, 2 random `u`).

### 4.2 Per-eigenspace residuals (sample)

`||R_i||^2 = ||R_{m-i}||^2`, and for `m < q` they are all strictly positive in every toy below saturation; e.g.
`sdc19`, `L = 2`, `m = 5`: `||R_i||^2 = 1.0426e6, 1.8228e5, 8.4259e4` (`i = 0, 1, 2`); the extreme monomials
`khat_1^m`, `khat_{-1}^m` always carry the largest residual and the balanced ones `|khat_1|^{m}` the smallest. At `sdc11`,
`m = 8` all `R_i` vanish except `R_4` (the trivial-character one): the near-saturation toys concentrate the residual in
the trivial eigenspace, as barrier §6 found for the tangent forms.


## 5. Scaling in `m`, `L`, `q`; comparison with the cap threshold (evidence)

Read off from the tables of §4 (exact entries marked `E`, floating `F`).

1. **Constancy on the circle.** `Psi_A` is exactly `theta`-free (Prop 2.1). The LP values show that `kappa_V(hat u_theta)` is
   *almost* constant too: `inf/sup >= 0.994` on `sdc7, L = 2` (`m = 2..6`), `>= 0.947` on `sdc11, L = 2` (`m = 2..4`),
   `>= 0.856` on every `L = 1` box. Barrier Thm 3.2 only guarantees `(m+1)^{-1/2}` (`0.38..0.58` here). So "somewhere =
   everywhere" holds on these circles with essentially no loss. The LP optimum `Cheb_A` exceeds `Psi_A^{1/2}` by a factor
   `1.56 - 5.67` over all LP cases, against the a-priori room `U^{1/2}` (`7.5 - 97`) of (1.1).
2. **The circle is not a special (more tangent) direction.** Relative to 6 random `u` (`|u|_inf = 1`, uniform in the cube):
   on `sdc11, L = 2` the circle's `rel` exceeds the largest random value at every `m = 2..8` (factor `1.03 - 1.43`, and
   `2.4` at `m = 8`); on `sdc7, L = 2` it is above the random range for `m = 2, 3` and inside it for `m = 4, 5, 6` (between
   `1.2 x` the smallest and `0.32 x` the largest random value at `m = 6`); on the pair boxes `L = 1` it sits at the bottom
   of the random range, at most a factor `0.88` below it (`q = 43`, `m = 3`). On `sdc19, L = 2` the circle is above the
   random range at every `m = 2..7` (factor `1.01 - 1.15`; only 2 random `u` at `m = 7`); on `sdc23, L = 2` it is above it for `m = 2, 3, 4, 6` (factor
   `1.00 - 1.03`; 2 random `u` at `m = 6`) and inside it at `m = 5` (`0.1035` against `[0.0882, 0.1060]`). No toy shows a tendency of the Fourier
   circle to approach `TC_V` faster than a random direction; for large `q` at fixed `L` the circle and the random directions
   become indistinguishable at this resolution (both sit at the cube value, item 5).
3. **Dependence on `m`.** `-ln(rel)/m` grows slowly with `m` and is of the same size as for the **uncut cube** on the
   same circle: `rel_A/rel_cube in [0.29, 1.31]` for every `m <= rho(A) - 2` with `m < q` (`sdc7`, `sdc11`, `L = 2`), and
   `rel_A >= rel_cube` for `m <= 6` at `q = 11`. The per-degree decay `-ln(rel)/m` is `0.2 - 0.86` up to `m = rho - 2`; the
   only fast drop is at the last degree before saturation (`sdc11`, `m = 8 = rho - 1`: `1.50` per degree, `rel = 6.0e-6`),
   where `dim V(A,m) = 5` and the toy is underdetermined in the record's sense.
4. **Dependence on `L`.** `L = 1` (the pair box, `rho = 4`) decays faster at equal `m` than `L = 2` (`rho = 9` for
   `q = 7, 11`): `rel = 0.553` vs `0.676` at `m = 2`, `0.183 - 0.222` vs `0.373 - 0.395` at `m = 3`. In both cases the decay
   is governed by `m/rho(A)`, i.e. by the distance to saturation, not by `m/(2L)`.
5. **Dependence on `q`.** At fixed `(L, m)` the values converge as `q` grows: `L = 1` gives `rel -> 0.5527, 0.1833`
   (`m = 2, 3`, `q = 103`), monotone in `q` from `q = 7`. For `L = 2` the `q = 7, 11` values of `rel` agree to within a factor
   `1.41` for `m <= 4` and differ by `2.7` at `m = 5` (both boxes have `rho = 9`, but `sdc7` has `U = 757` against `9373`).
   Going on to `sdc19` (`F`) the `L = 2` values stabilise and lock onto the **uncut cube**: `rel(C_1) = 0.659, 0.381, 0.209,
   0.104, 0.0429, 0.0124` for `m = 2..7`, i.e. `rel_A/rel_cube = 1.01, 1.02, 1.04, 1.03, 0.88, 0.55`, and within a factor `1.17` of `sdc11` for
   every `m <= 6` (at `m = 7` the slab starts to bite at `sdc19` too: `rel_A/rel_cube = 0.55`, `dim V(A,7)/U = 0.51`; the
   `sdc11` value `0.0065` is lower still, `sdc11` being two degrees from saturation). `sdc23` (`E` for `m <= 4`, `F` at `m = 5, 6`) continues the pattern: `rel = 0.658, 0.381, 0.2075, 0.1035, 0.0436`
   (`m = 2..6`), `rel_A/rel_cube = 1.00, 1.00, 1.01, 0.98, 0.83`; at `m = 6`, `-ln(rel)/m = 0.551, 0.525, 0.522` for `q = 11, 19, 23`; at `m = 5` the per-degree decay `-ln(rel)/m` is `0.447, 0.453, 0.454`
   for `q = 11, 19, 23` (converging in `q`). In words: at fixed `(L, m)` and growing `q` the slab constraint becomes invisible to the
   circle at degree `m` (the box is far from saturation, `dim V(A,m)/U >= 0.85` at `sdc19`, `m = 6`), and the circle's
   distance from `Poly_{<m}` is the cube's, `rel ~ e^{-(0.2 .. 0.63) m}` (`0.627` per degree at `sdc19`, `m = 7`), slowly approaching the cube's `ln 2`-per-degree
   asymptotics.
6. **`kappa_H`.** Where computed (`L = 1`, `q <= 23`, `m = 2, 3`; `sdc7`, `L = 2`, `m = 5, 6, 8`) an LLL-reduced integer basis
   of `V(A,m)` has height `H_red = 1` (`L = 1`), `3, 4` (`sdc7`, `m = 5, 6`) and `42` (`sdc7`, `m = 8 > q`). Hence
   `kappa_H = kappa_V` for every `H >= H_red` (exact), and already `kappa_H(hat u, H = 1) >= 0.28 sqrt(Psi_A)/d^m` on
   `sdc7, m = 5, 6`. No toy shows an `H`-gap on the circle at heights `>= 4` (`m < q`).
7. **Against the threshold shape.** Writing `kappa = eps_0(W)/(2 (m+1)^{1/2} e^{pi^2/(4q)})` and solving for `W` gives the
   barrier scale `W* = 2q tau*` at which each toy circle would sit exactly on the CIRCLE threshold (with `m` in place of
   `m'`). Using the rigorous lower bound `k_lo` of (2.1), `tau* <= 2.22` in every toy (and `<= 1.38` with the LP value),
   while the barrier needs `tau >= tau_D = 164 - 221` (`q = 7..103`). In logarithms: `ln inf kappa_V(C_1)` exceeds the
   threshold at the barrier scale by `q(tau_D - tau*)/2 >= 568` (`q = 7`) and `>= 945` (`q = 11`). The toy circles
   are nowhere near a cap; `ln inf kappa_V` decreases only linearly in `m` (slope `-1.1` to `-1.6` per degree for `L = 2`,
   except `-4.1` at the last degree before saturation of `sdc11`): the shape `-m ln d + (1/2) ln E_A l^{2m} + O(m)` forced by
   the normalisation `d^m`, with no sign of the super-linear decay a cap needs (§7).


## 6. Consistency with barrier.md's exact `sdc7 / sdc11 / sdc19` data

| record datum (barrier §6, referee_barrier §6) | this round | agreement |
| --- | --- | --- |
| `sdc7`, `L = 2`, `m = 8`: `|A| = 757`, `dim V = 1`, the octic `c` translation-invariant, `Phi_c = -7717.5` constant on `C_1, C_2, C_3`; circles not in `TC` | `full_check.py 7 2 8`: rank 756, `dim V = 1`; `Psi_A = 23.1272016733` at all 5 `theta` (so constant, although `m = 8 > q` is outside Prop 2.1); LP `Cheb = 7.65625` at every `theta`; `kappa_h.py`: primitive integer `c` with `|c|_inf = 42`, `|c|_1 = 1008`, and `7.65625 x 1008 = 7717.5 = |Phi_c|` | **yes** (exact `Phi_c` reproduced) |
| `sdc7`, `L = 2`: `psi_1 notin TC(A,i)` for `i <= 6` (author, mod `p`) | `Psi_A(C_1) > 0` exactly for `m = 1..6` (rationals in `out/exact_q7_L2.txt`) | **yes** |
| `sdc11`, `L = 2`, `m = 8`: `dim V = 5`, deficiencies `delta_0 = 5`, all others `0`; `F_4` not in span, all other `F_i` in span; circles not in `TC` | exact: `dim Poly_{<8}|_A = 9368`, `dim V = 5`, `delta_0 = 5`, `delta_c = 0` (`c != 0`); `||R_i||^2 = 0` for `i != 4` and `||R_4||^2 = 32923557125687554191360/81681557465859169 ~ 4.0307e5`; `Psi_A(C_1) = 1406646117430498531125/46675175694776668 ~ 3.0137e4` | **yes**, and now quantitative |
| `sdc19`, `L = 2`, `m = 5, 6`: no `F_i` in `Poly_{<m}`; deficiencies `13464` (`m = 5`), `11886/11887` (`m = 6`) | `float_iso_sym.py 19 2`: per-eigenspace ranks `dim Poly_{<5,c} = 446` (`c != 0`), so deficiency `13910 - 446 = 13464`; `dim Poly_{<6,c} = 2024` for `c = 0` and `c != 0`, deficiencies `13911 - 2024 = 11887` and `13910 - 2024 = 11886`; every `||R_i||^2 > 0` (`m = 5`: `1.0426e6, 1.8228e5, 8.4259e4`; `m = 6`: `3.1272e7, 4.2831e6, 1.4277e6, 1.0218e6`), `Psi_A(C_1) = 2.7393e4` (`m = 5`), `3.4719e5` (`m = 6`) | **yes** (deficiencies identical to `barrier_checks/sdc_circles_q19_L2_out.txt`; floating, so evidence; quantitative now) |
| `sdc11`: `rho(A) = 9` implicit (`dim V(A,8) = 5`) | exact `dim V(A,9) = 0` | **yes** |

So all exact circle data of the record are reproduced, and each "circle not in `TC`" becomes a number: the
circles sit at relative residual `rel = 6.0 x 10^{-6}` (`sdc11`, `m = 8`, the last degree before saturation) and higher
everywhere else. The record's only circle zeros-by-saturation are the trivial ones (`m >= rho(A)`).


## 7. Extrapolation to the regime (LABELLED EVIDENCE)

### 7.1 What a "yes" would require, through `kappa_V` (exact bookkeeping)

`regime_need.py` (output `out/regime_need.txt`; `C = 1`, `c_0 = 0.184`, `tau = tau_D` of barrier Lemma 3.4, `d = q`). By
(2.1), `inf_theta kappa_V(hat u_theta) <= T := eps_0/(2 sqrt(m'+1) e^{pi^2/(4q)})` **forces**
`ln Psi_A(C_1) <= 2(ln T + m' ln d)` (and is implied by that bound minus `d ln 3`, by (1.1)):

| `q` | `j` | `m'` | `m'/q` | `tau_D` | `W = 2q tau_D` | `ln T` | needed `ln Psi_A(C_1)` | heuristic `ln E_A l^{2m'}` | needed `-ln(rel)/m'` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 3251 | 257 | 2057 | 0.633 | 296.1 | 1.93e6 | -4.813e5 | <= -9.293e5 | 2.94e4 | **>= 466** |
| 5479 | 342 | 2738 | 0.500 | 307.4 | 3.37e6 | -8.422e5 | <= -1.637e6 | 4.14e4 | >= 613 |
| 10007 | 475 | 3803 | 0.380 | 320.5 | 6.42e6 | -1.604e6 | <= -3.138e6 | 6.10e4 | >= 841 |
| 100003 | 1648 | 13191 | 0.132 | 370.6 | 7.41e7 | -1.853e7 | <= -3.676e7 | 2.58e5 | >= 2806 |
| 1000003 | 5638 | 45107 | 0.045 | 420.6 | 8.41e8 | -2.103e8 | <= -4.194e8 | 1.04e6 | >= 9321 |

The column "heuristic `ln E_A l^{2m'}`" uses a Gaussian model of variance `q/3` (the cube value of `Var <k,u_theta>`; regime
`A` has density `>= 1/2`); it only converts the needed `Psi` into a needed relative residual `rel = Psi/E l^{2m'}`, and the
conclusion does not depend on it within any reasonable factor (`ln E` is `< 3.2%` of `|needed ln Psi|` in every row).
Even with the Diophantine bookkeeping improved to its floor `tau ~ 2 ln(4r)` (the prime-window requirement), the
needed decay at `q = 3251` is still about `30` per degree (`tau = 20.3`).

### 7.2 What the toys show

Across all toys with `dim V(A,m) >= 2` and `m < q` (§4):
* the relative residual on the circle satisfies `-ln(rel)/m <= 0.86` whenever `m <= rho(A) - 2`, and `<= 1.51` overall
  (the single exception being the last degree before saturation of `sdc11`);
* this is the uncut cube's behaviour on the same circle (`rel_A/rel_cube in [0.29, 1.31]`), whose Gaussian-model decay
  tends to `ln 2 = 0.69` per degree;
* the circle is not more tangent than a random direction (§5.2);
* `V_H = V(A,m)` already at heights `H <= 4` (`m < q`), so in the toys `kappa_H` on the circle is `kappa_V` up to a factor
  `<= 4` even at `H = 1`.
The toy ratios `m/(2L)` reach `2` (`sdc11`, `m = 8`; `1.75` at `sdc19`, `m = 7`), against the regime's `m'/(2j) in [4, 4.04]`; the toy
ratios `m/rho(A)` cover `[0.1, 0.9]`, while in the regime `m'/rho(A) < m'/(q+1) <= 0.633` (proved below) and `-> 0`.

**Lemma 7.1 (proved).** If `A subset {-1,0,1}^d` has `|A| >= 3^d/2`, then `rho(A) >= d + 1`.
*Proof.* `Poly_{<rho}|_A = R^A` needs `#{alpha in {0,1,2}^d : |alpha| <= rho - 1} >= |A| >= 3^d/2` (the reduced monomials of degree
`< rho` span `Poly_{<rho}` on the ternary cube). The distribution of `|alpha|` under the uniform measure on `{0,1,2}^d` is
symmetric about `d`, so `#{|alpha| <= d - 1} < 3^d/2`. Hence `rho - 1 >= d`. QED
(Regime boxes have `|A| >= 3^d/2`, t1rigour §2; so `m'/rho(A) < m'/(d+1)`, which is `0.633` at `q = 3251` and decreases.)


### 7.3 Extrapolated reading

**Labelled extrapolation (evidence, not a proof).** If the regime circle behaves like every toy at comparable
`m/rho(A)` -- relative residual decaying at a rate of order `ln 2` to `1` per degree -- then at `q = 3251`

```text
ln inf_theta kappa_V(hat u_theta)  ~  (1/2)(ln E_A l^{2m'} - c m') - m' ln d  ~  (1/2)(2.94e4 - 2.06e3) - 2057 x 8.087  ~  -2.9e3   (c ~ 1)
ln T                               =  -4.81e5 ,
```

so the circle would miss the cap threshold by a factor `e^{-4.8e5}`, i.e. by `~ 470` per degree of `m'`, and the gap
widens with `q` (§7.1 table: `>= 9321` per degree at `q ~ 10^6`). Through `kappa_V`, CIRCLE(q) would then be a robust
**no**. Three caveats keep this from being more than evidence:
1. **The `H`-gap.** A "yes" needs `kappa_H` small, not `kappa_V`. The toys have `V_H = V` at heights `<= 4`, which
   suggests (P1)-type behaviour on the circle; nothing is known at the regime's `ln H' = m'W/40`.
2. **Ratio and density.** No toy reaches `m/(2L) = 4` with `dim V >= 2`, and no toy has density `>= 1/2`
   (`|A|/3^q <= 0.35`). The extrapolation uses `m/rho(A)` as the governing parameter, which the data support (§5.3-5.4),
   but this is not proved.
3. **Exact vanishing is invisible to calibration.** A degree-`m'` identity in `J(A)` along `psi_1` would be a
   regime phenomenon with no toy shadow; Prop 2.1 reduces it to `||R_i|| = 0` for all `i`, which the toys never show
   below saturation.

**What a proof of "no" would need**, as the data indicate: a lower bound
`Psi_A(C_1) >= e^{-o(qtau)}` (the needed `ln Psi` is `-9.3e5` at `q = 3251`; the toys suggest the truth is `~ +2.7e4`), or,
through Prop 2.1, `||R_i||^2 >= e^{-o(q tau)}` for **one** `i` (all `||R_i||` enter with positive weights) -- a lower bound
for the distance of the single eigenfunction `khat_1^i khat_{-1}^{m'-i}` from `Poly_{<m'}` in one translation eigenspace,
together with (P1)/(P2) to pass from `V` to `V_{H'}`. The toys say the margin for such a bound is enormous (about `9.6e5` in `ln Psi`, i.e. a
factor `e^{4.8e5}` in `kappa`, at `q = 3251`), so a crude certificate may suffice; none was attempted here.



## 8. Scripts and outputs (`calib_checks/`)

| file | content |
| --- | --- |
| `enumA.c` | exact enumeration of `A_L(q)` (translation-orbit representatives to `data/A_q{q}_L{L}.bin`) |
| `exact_iso.py` | exact multimodular isotypic computation of `||R_i||^2`, `||F_i||^2`, `Psi_A(C_1)`, dims (`out/exact_q*_L*.txt`) |
| `full_check.py` | independent full-space floating check (no symmetry), five `theta`, random `u`, small LP |
| `float_iso.py` | complex128 isotypic engine with random-`u` comparison (`out/float_q*_L*.txt`) |
| `cube_psi.py` | `Psi_cube` and `E_cube l^{2m}` on `C_1`, 50 digits (`out/cube_q*.txt`) |
| `cheb_lp.py` | `Cheb_A = d^m kappa_V` on a `theta`-grid by LP (`out/cheb_q*_L*_m*.txt`) |
| `kappa_h.py` | integer basis of `V(A,m)`, LLL heights, `kappa_H` lower bounds (`out/kappaH_*.txt`) |
| `regime_need.py` | regime bookkeeping of §7.1 (`out/regime_need.txt`) |
| `float_iso_sym.py` | complex128 engine using the `AGL(1,q)` symmetry (only eigenspaces `c = 0, 1` built); `sdc19`, `sdc23` runs (`out/float_q19_L2*.txt`, `out/float_q23_L2*.txt`) |
| `exact_iso_big.py`, `exact_iso_big23.py` | `exact_iso.py` with primes `< 2^24` / `< 2^23` (exact `int64` dot products) for `sdc19`, `sdc23` (`out/exact_q19_L2*.txt`, `out/exact_q23_L2.txt`) |
| `tables.py`, `insert_table.py` | assemble the §4 table from the outputs and insert it here |

Data files under `data/` are regenerated by `enumA` (largest: `q = 23`, `L = 2`, 0.8 MB).
