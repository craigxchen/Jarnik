# Referee report on round11/subcluster.md (adversarial)

Referee: adversarial review of `round11/subcluster.md` (cited **[S]**), against `round10/determinant.md` (**[D]**) and
`round10/referee_determinant.md` (**[RD]**). I re-derived every proof and re-did the computations with my own
code. I did not import the author's helpers. The only shared code is the round-10 referee's fake-profile
constructors, used read-only. Scripts and outputs are in `round11/referee_subcluster_checks/` (`rc_*.py` writes
`rc_*_out.txt`).

## 0. Verdict summary

**The result survives as stated at the top level: no growth improvement and no uniform bound.** Every theorem
I checked holds, except Lemma 5.3(i)-(ii): its proof fails under the stated ZE-tor definition. Separately, the
hypothesis of the conditional Theorem 5.4 is over-quantified. Both have easy fixes, given in §7 and §10. The interpretation and the computational claims are weaker than the write-up says.

Correct as proved (I re-derived each):
* Lemma 1.0, Lemma 1.1, Lemma 1.2;
* Proposition 2.1, Proposition 2.2;
* Theorem 3.1, Proposition 3.2 (i)-(iii);
* Theorem 4.1 (a)-(b), Corollary 4.2, Proposition 4.3;
* Theorem 4.4 (the new signed mean-width bound), which I also tested exactly; Corollary 4.5;
* Theorem 5.4 as a logical implication.

Problems, most severe first (full list in §9):

1. **(Major) ZE-col(univ), the hypothesis of the headline conditional Theorem 5.4, is quantified over every profile
   and every nondegenerate point of `L`.** In that form it is not merely unproved: ZE-col is false for thin
   relation families. I give a rigorous rank-3 example (§7.2): one integer measure whose zero set has
   nondegenerate points arbitrarily close to `e`. Nothing in [S] rules out profiles whose universal family is
   thin on `L`. **Fix.** Theorem 5.4's proof needs ZE-col only for *robustly visible* profiles: those with
   `rho_hat >= 300 ln 4r`, or with `M >= (rho_hat+1)(rho_hat+2)` and the whole degree-2 kernel admissible. It also
   needs it only at the actual position.
2. **(Major) The ZE-tor definition and Lemma 5.3(i)-(ii).** At the fixed radius `1/(4D)`, "coincides with a finite
   union of rational subtori through `e`" allows non-identity sheets, and torsion cosets of intersections. The
   proof step "near `e` a closed subgroup coincides with its identity component" is false at a fixed radius.
   Proposition 5.1 and the ZE-tor half of Theorem 5.4's discussion inherit the gap. **Fix.** State ZE-tor
   locally, in the lifted ball, as a finite union of `V_i ∩ ball` with `V_i` rational subspaces containing `1`.
3. **(Major, interpretive) "Nondegenerate sub-clusters give the same system in other coordinates, so they add
   nothing" is unsupported.** It also contradicts [S]'s own Remark (i) of §1.2. Theorem 3.1(c) concerns *one*
   sub-cluster. In the full-rank case (`rho_hat = M - 1`), which is the case [D] Theorem 5.7 reduces to through its
   nondegenerate sub-cluster, the push-forward is the identity. The intersection over many nondegenerate `J` is then exactly the question of the round, and
   Theorem 3.1 says nothing about it. The headline "reduces exactly to ZE-col" is, in that case, a restatement.
4. **(Major, computational claim) "Below `n ~ 2137` the honest combined system is empty" is false as a general
   statement.** The author tested only `n = 4, 5`. By exact computation (§8.4), the [D] slab-box with `T = n` is
   LP-supercritical for every `n` from 440 up to the end of my range (2400) at `r = n`, so `V(A,m) != 0`. Its integer
   generators are admissible in the sense of [S]'s own Lemma 1.1 for every `n` from 468 up to 2400 at
   `W = 3n ln n - 14n`. The number 2137 is only [D]'s factor-2 margin `m >= 4w`.
5. **(Moderate) The §6.5 evidence is structurally uninformative.**
   * "Fakes and random positions certified off every locus" is automatic: the loci are Lebesgue-null.
   * "No tangent-cone zeros" is *exactly* true in all 560 systems, but for a trivial reason. There `m = 6 < 2T + 1`,
     so no polynomial of degree `<= m` vanishes on `A` (§8.3).
   * In the regime `m >= 0.7 n(n-1) >> 2n + 1`, and the tangent cone always contains the sparse and column
     directions.
6. **(Minor)**
   * The "Reading" of Proposition 3.2 (products are no stronger than characters) holds only at leading order.
   * The summary's "fibres are finite" (Proposition 2.2 allows whole-circle fibres).
   * Proposition 5.1's equivalence is close to tautological.
   * `R >= R_0(C)` must refer to the normalised radius, with a constant bound below it.
   * The word "certified" in §6.5, and the "verified with fractions" header of `lowrank_profiles.py`.
   * "`n = 4..6`" was tested only at `n = 4, 5`.

## 1. Setup and definitions checked

* **Windings and (1.1).** From `arg z_x = <2a_x - e, phi> + (pi/2) j_x` and `arg z_x = theta_0 + delta_x mod 2pi`,
  one gets `delta_x = <hat a_x, hat phi> - c_0` with `hat phi = (2 phi, -pi/2)`, `c_0 = theta_0 + <e, phi>`.
  * The bound `|kappa_x - kappa_b| <= 2 sum e_p + 1` follows from `0 < phi_p < pi/2` and `Theta < 1`. Correct.
* **`L` and `rho_hat`.** Write `hat A` for the matrix with rows `hat a_x`. Then
  `dim L = 1 + rank(hat A - 1 hat a_b^T) = 1 + dim span{hat a_x - hat a_b}`, so `rho_hat` is the affine rank of
  `{hat a_x}` and `rho_hat <= rank + 1`. Correct. The coordinates used in Theorem 4.1(b) are consistent with this.
* **`ord_L`, `LP`.** `ord_L(c) >= ord(c)` because `P o hat pi` has degree `deg P`. On the slice,
  `<k, beta_p> = (2/e_p) <k, a_p>`, so `LP(c) = (2/W) sum_p log p h_{supp c}(a_p)`. Correct.
* **Normalisation.** [S] works with normalised clusters (`N` odd, all primes split). Under `z = g z'` the arc
  constant becomes `C' = C/sqrt|g| <= C`, but the normalised radius may be *small*. See §7.3.

## 2. Lemma 1.0 / 1.1 (universal admissibility, nondegeneracy relative to L)

**Lemma 1.0: correct.** If `u = sum v_x a_x != 0`, then `alpha = prod pi_p^{u_p}` has
`alpha/conj(alpha) = prod (pi_p/conj pi_p)^{u_p}`. This has nonzero `pi_p`-valuation for some `p` (`p` odd split,
so `pi_p` and `conj pi_p` are non-associate). So it is not a unit, `2 arg alpha notin (pi/2)Z`, and
`<v,delta> != 0`. This is the `(pi/4)Z` gap, correctly cited, as [RD] required for [D] Lemma 3.2. Minimality of `L`:
* any rational `V` containing `1` and `delta` has `V^perp` spanned by integer vectors in `Z^M_0` that annihilate
  `delta`;
* these lie in `L^perp`, so `L subset V`.

**Lemma 1.1: correct.** I redid the chain.
* Taylor: `sum c_k <hat pi(k), hat phi>^j = 0` for `j < ord_L(c)`, which gives
  `|g_c| <= |c|_1 (D Theta/2)^m/m!`.
* Forced divisibility: `|g_c|^2 = Norm F_c/N^D >= e^{-(W/2) LP(c)}`.
* `ln Theta <= ln C' - W/4`.
* The contradiction with (1.2) follows.

Two remarks, neither an error:
* the valuation of `G_k(z)` depends on `k` only through `hat pi(k)`, so `LP(hat pi_* c)` could replace `LP(c)` in (1.2);
* the lemma holds with any `m <= ord_L(c)`, which Theorem 4.1(b) uses implicitly.

**Lemma 1.2: correct.** Each sub-cluster `J` is a cluster with the same `N`, `W`, `r` and arc, so [D] Theorem 4.1
applies. Sub-cluster generators also satisfy (1.2): `f(m)` is increasing once `W/2 > 2 ln(DC'/2)`.

**Remark (i) of §1.2 is correct and important.** `V(A_J, m_J)` is not contained in `V(A_{J'}, m_{J'})`, so
sub-cluster systems are genuinely new relations. This is in tension with the later slogan "nondegenerate
sub-clusters add nothing" (§4 below, issue 3).

## 3. Props 2.1-2.2 (dimension and fibres of the combined locus)

**Proposition 2.1: correct.**
* Each `Z_{P_J}` is real algebraic in `(cos, sin)` of the differences, including the other base points. It is
  null, so its dimension is `<= n - 2`.
* At a smooth top-dimensional point of `Z_comb^{(n)}`, `n - 1` of the functionals `delta_x - delta_b` are
  independent on the tangent space. So `pi_J` is a submersion for `J = {b, x_1, ..., x_{n-1}}` with base `b`, which
  is a member of the family, and its image has interior.
* As the author says, the statement is irrelevant for one point.

**Proposition 2.2: correct.** For `x` not the base point, `g_c` is a Laurent polynomial in `zeta_x` with exponents in
`[-T, T]`, so the fibre is the whole circle or has `<= 2T` points.
* *Wording.* The task summary's "its fibres are finite" drops the whole-circle alternative. The write-up itself
  states it.

## 4. Thm 3.1 / Prop 3.2 (push-forward, no product certificates)

**Theorem 3.1: correct.**
* (a) `<k, delta'(s)> = <hat pi(k), s>` on `Z^M_0`, and distinct characters are independent.
* (b) `ord_L(c) = ord(hat pi_* c)` holds by definition, and the widths shrink.
* (c) If `{hat a_x}_J` is affinely independent, then `hat pi` is injective on `Z^J_0`, with a linear left inverse.

**The "Consequence" paragraph and the slogans are overstated (issue 3, major).**
* Theorem 3.1(c) says that *one* nondegenerate sub-cluster's system, pushed to `L`, is isomorphic to itself. It
  says nothing about the *intersection* of several such systems.
* In the full-rank case (`L = R^M`) the push-forward is the identity. [D] Theorem 5.7 reduces squarefree
  clusters to this case: Lemma 3.1 gives rank `>= M/4 - 2`, and Lemma 3.2 a nondegenerate full-rank sub-cluster.
* There `Z_comb ∩ L = Z_comb`, and its intersection over `J` is precisely the object of the round. [S]'s own
  Remark (i) of §1.2 shows it is genuinely smaller than any single `Z_{P_J}`.
* Therefore the slogans "nondegenerate sub-clusters add nothing" (§0 table; task summary item 3) and "the
  question reduces exactly to ZE-col" carry no content in the full-rank case. In that case ZE-col(comb) *is* the
  question "does `Z_comb` contain nondegenerate points near `e`", word for word.
* What is true: the push-forward isolates the extra structure of *degenerate* sub-clusters (fibre sums of `c`).
  This is all that §3 proves.

**Proposition 3.2: (i)-(iii) correct.**
* (i) is the iterated-difference identity.
* (ii) `mu_alpha` annihilates degree `< |alpha|` and pairs with `prod ell_i^{alpha_i}` to give `prod alpha_i! != 0`,
  so `ord = |alpha|`. Also `h_{supp c}(a_p) >= h_Pi(y_p) = sum alpha_i |<v_i,a_p>|`.
* (iii) is (ii) inserted into (1.2).
* The deduction that an admissible product push-forward cannot exist for an actual cluster (some `zeta_i = 1`,
  hence `lambda_i = 0`) is correct.

**The "Reading" is a leading-order statement only (minor).** The right side of (iii) is at most
`2|alpha| ln(2|alpha|/(eDC')) + O(ln|alpha|)`.
* For a push-forward supported in `hat pi({|k|_1 <= D})`, `|alpha|` can be as large as about `D Y/min|lambda_i|`,
  with `Y <= 2W/ln 5 + 1`. So the per-unit slack can reach `+2 ln(W/(eC'))`, and products are admissible whenever
  the `alpha`-average of `X_{v_i}` lies below `W/2 + O(log W)`.
* Characters exclude only `X_v < W/2 - 2 log(C|v|_1/(2 sqrt 2))`.
* So inside the window `[W/2 - O(log), W/2 + O(log W)]`, product certificates could in principle exclude profiles
  that characters do not.
* That window is the scale at which the *constant* in `log R/loglog R` is decided: a Plotkin-type argument with
  `d_xy >= W/2 + K log W` gives `M <~ W/(2K log W)`.
* Since `ln(|alpha|/D) <= ln W + O(1)`, the extra slack is `O(log W)` per unit. So an `o(log R/loglog R)`
  improvement from products is not expected. But "single characters already decide" in that window, and "the R1
  full-box mechanism never appears", are not proved.
* Correct reading: product push-forward certificates can only move the character threshold by `O(log W)` per unit
  of `|alpha|`. This is a heuristic statement about the order of growth, not a theorem.

## 5. Thm 4.1 / Cor 4.2 (visibility, unconditional bound for blind clusters)

**Theorem 4.1(a): correct.** A maximal affinely independent `J_0` has `rho_hat + 1 >= 300 ln 4r + 1 >= 417 > 64`
points. [D] Theorem 4.1 applies to it, `V_H(J_0) != 0`, and Theorem 3.1(c) gives `hat pi_* c != 0`. So `g_c|_L` is a
nonzero trigonometric polynomial.

**Theorem 4.1(b): correct; every constant re-derived.**
* *Form and support.* With `D = 2` and support `{e_x - e_b}` (including `k = 0` for `x = b`),
  `G_{e_x - e_b} = z_x conj(z_b)` and `G_0 = N = z_b conj(z_b)`. So `F_c = conj(z_b) sum_x c_x z_x`.
* *Equations.* `ord_L >= 3` amounts to `E = binom(rho_hat + 2, 2)` integer equations in `M` unknowns, in
  `rho_hat` coordinates restricted from `Z^{r+1}`. The entries are bounded by `B = Y^2`, with
  `Y <= 2 sum e_p + 1 <= 2W/ln 5 + 1 <= 2W`.
* *Siegel.* Since `(2MHB+1) <= (2MB+1)(H+1)`, the condition `(H+1)^{M-E} > (2MB+1)^E` suffices, and it holds for
  `H = 2MB + 1` when `M >= 2E`. Then `|c|_1 <= M(2MB+1) <= 3M^2 B`.
* *Admissibility.* `LP(c) <= sum omega_p (max_x beta_{p,x} - min_x beta_{p,x}) <= 2`. The right side
  `f(m) = m(W/2 - 2 ln C') + 2 ln m! - 2 ln|c|_1` is increasing. `f(3) > W` is equivalent to
  `W > 12 ln C' + 8 ln M + 8 ln W + 4 ln 2`, which is implied by the stated `W > 8 ln M + 8 ln(2W) + 12 ln C'`.
* *Distinct points.* `hat pi_* c` has distinct atoms: the `a_x` are distinct, since `z_x/z_y` would otherwise be a
  unit on a short arc.

So actual clusters with `M >= (rho_hat+1)(rho_hat+2)` satisfy a nontrivial relation `sum c_x z_x = 0` with
`ord_L >= 3` and small height. (This is the classical "short-arc points satisfy small linear relations"
mechanism.)

**Corollary 4.2: correct**, including `(rho_hat+1)(rho_hat+2) < (300 ln 4r + 2)^2` and `r <= W/ln 5`.
* The phrase "holds for `R >= R_0(C)`, because `M <= W`" must be read for the *normalised* cluster (§7.3).
* "Blind" is a property of the profile together with the windings. The corollary is unconditional, but it
  applies only to blind clusters. Nothing proves that actual clusters are blind. That is exactly what ZE-col
  would supply.

**Heuristic remark on sumsets: correctly labelled.** The Blichfeldt-type bound
`|P ∩ Z^d| <= d! vol(P) + d` is used in its correct direction.

**Side observation (not an error).** The threshold `300 ln 4r` in Theorem 4.1(a) is not intrinsic to [S]'s
framework.
* Lemma 1.1 admits a slab-box generator as soon as `m > 2w` and `W(m/2 - w)` beats the error terms.
* By §8.4, this happens from `n ≈ 440-470` at `r = n`, against `300 ln 4r + 1 ≈ 2700`.
* So Corollary 4.2 and Theorem 5.4 could be stated with a smaller (exactly computable) threshold. This does not
  change any order of magnitude.

## 6. Thm 4.4 / Cor 4.5 (signed mean-width inequality)

**Theorem 4.4: correct. I re-derived the induction and tested it exactly.**
* *Base cases.* `rho = 0` is trivial. For `rho = 1`, `reg = |S| - 1 <= max - min = Omega'`.
* *Key inequality.* For slices `c <= c'`:
  * `h(<eps',y'> + y_1) >= (c' + b_2) - (c + a_1)` and `h(<eps',y'> - y_1) >= (a_2 - c) - (b_1 - c')`;
  * these add to `2(c' - c) + l_c + l_{c'}`, and averaging gives `K >= (c' - c) + (g_c + g_{c'})/2`;
  * `c = c'` gives `g_c <= K`.
* *Counting and peeling.* `|C_v| <= K - v + 1` (slices sit at integer heights). The peeling lemma
  `reg(S) <= max_i (reg S_(i) + i - 1)` is the refereed Lemma 2.1 of `threshold/construct.md`. I re-derived it:
  interpolate on `S_(1)`, then factor out `(y_1 - c_(1))`.
* *Comparison with Theorem B.* `MW_rho(S) = E_eps (1/2)[h(<1+eps,.>) + h(<1-eps,.>)] >= E_eps h(<eps,.>)` by
  subadditivity. This matches `threshold/construct.md` §2.4 (`MW_n = 2 E_P w_P`, with `1_P = (1+eps)/2`). So
  Theorem 4.4 does strengthen Theorem B.
  * Note that `Omega'` is *not* `GL_rho(Z)`-invariant. The bound depends on coordinates, and in the slice
    formulation it depends on the dropped base coordinate.
* *Novelty.* I cannot check the literature offline. The proof is the Theorem B induction with signed weights,
  so "new" should be read as "new in this project".

*Exact test* (`rc_thm44.py`, my own exact-rational rank code; output `rc_thm44_out.txt`):
* 950 random sets in dimensions 1-4, and 50 hill-climbing restarts maximising `reg - Omega'` in dimensions 2-3:
  **no violation**. The best hill-climb value is `reg - Omega' = 0`.
* `l^1`-balls (`rho = 2, 3`; `R = 2, 3`): `reg = Omega' = 2R`, with `MW = 6, 9, 7, 10.5`. Boxes: `reg = Omega' = sum a_i`.
  This reproduces [S]'s equality claims and its `MW = 10.5` example.

**Corollary 4.5: correct, and vacuous for actual clusters, as the Caveat says.**
* *The chain.* `4w_J <= m_J <= ord(c) <= ord_L(c) = ord(mu) <= reg(supp mu) <= Omega'(supp mu) <= LP(mu)/gamma <= LP(c)/gamma <= 2w_J/gamma`
  forces `gamma <= 1/2`. `ord(mu) <= reg(supp mu)` holds because a nonzero `mu` cannot annihilate all functions on
  its support.
* *Weights in the Example.* With primes in `[P, P^{1.9}]`, `gamma >= 1/1.9 > 1/2`. But pair separation then
  fails, since the pair sums are `2 sum` of the weights of the separating classes and can drop well below `1`.
  The write-up concedes a loss "bounded by the weight imbalance". Only nearly equal weights give both properties.
* *Windings.* Lambda's coordinates must be given by columns, which requires `kappa` to be an integral affine
  function of the columns. [RD] §6 gives `kappa ∈ span{a_p} + R1` with *rational* coefficients, when the rank is
  small. Integrality is an extra assumption. Harmless, because the example has no admissible sub-cluster anyway.
* *Why the Corollary is restricted to slab-box systems.* With `gamma ∈ (1/2, 1)` the geometry is blind to
  slab-box systems (margin `m >= 4w`) but not necessarily to the universal family (margin `m > LP(1+o(1))`).
  Only `gamma = 1` makes it blind to everything. The Corollary is correctly restricted, but §6.2's "universal
  slab space" computation uses `gamma = 1` only.
* *§6.2 itself is fine.* The push-forward-zero claims are certified over `Q` by `V(S_slab, m) = 0`. This holds
  because `reg(l^1`-ball of radius `w) = 2w < m` in every run. I checked that every run has `m >= 2w + 1`.

## 7. Prop 5.1 / Thm 5.4 (rank target, conditional result)

**Proposition 4.3 (checked here because §5 uses it): correct.**
* `sum_{p in S} h_{A_J}(a_p) <= 3 w_J < 4 w_J <= m_J`, so [D] Proposition 5.2 with `kappa = 0` applies.
* The Caveat's deduction "cluster in `K_S` implies `rank <= 3`, hence `M <= 83`" silently uses Lemma 1.0:
  `delta ∈ K_S ∩ L` is nondegenerate in `L`, so `L subset span{a_p : p in S} + R1`. The deduction is correct.

### 7.1 ZE-tor and Lemma 5.3(i)-(ii): the definition must be local (major)

ZE-tor says: in the ball `max|delta'_x - delta'_b| < 1/(4D)`, `Z_F` coincides with a finite union of rational
subtori through `e`.

* A rational subtorus `K = V_K/(V_K ∩ 2pi Z^M)` of large height meets a fixed ball in **many parallel sheets**
  `(V_K + 2pi v) ∩ ball`, not only in `V_K ∩ ball`.
* The intersection of two subtori through `e` is a closed subgroup whose non-identity components (torsion
  cosets) can lie in the ball. Example: in `T^2`, `K_1 = {(t,0)}` and `K_2 = {(t, Nt)}` meet in the points
  `(2 pi j/N, 0)`. These lie within `1/(4D)` of `e` once `N > 8 pi D`.
* So Lemma 5.3(i)'s step "near `e`, a closed subgroup coincides with its identity component, a rational subtorus"
  is false at the fixed radius `1/(4D)`. It holds only in a neighbourhood that depends on the subgroup.
* Lemma 5.3(ii)'s step "`delta'` lies in `K`, the smallest rational subspace containing `delta'` is `L`, so
  `L subset K`" needs `delta'` on the *identity sheet* `V_K`. If `delta'` lies on another sheet, nothing follows.
* For actual `delta` (`|delta| <= Theta ≈ e^{-W/4}`), the extra sheets matter only when the height of `K` exceeds
  about `1/Theta`. The known height bounds for torsion cosets in a variety of degree `D = 2n(n-1)` are of order
  `D^{O(n)} = e^{O(n log n)}`, and that is comparable to `e^{W/4} >= e^{0.75 n log n}`. So the issue is not
  automatically harmless.

**Fix (sufficient for every use in [S]).** Define ZE-tor in the lifted ball `B ⊂ R^M/R1`, radius `1/(4D) < pi`:
`Z_F ∩ B = ∪_i (V_i ∩ B)` with finitely many rational subspaces `V_i` containing `1`.
* Then (i) is trivial: intersections of such unions are such unions, and `pi_J^{-1}(V ∩ B) = (V × R^{M-|J|}) ∩ B`.
* (ii) holds by the author's argument.
* Proposition 5.1, and the remark "for clusters visible to `Z_comb`, ZE-col(univ) follows from ZE-tor", then hold
  as stated. This local form is also the natural reading of `round10/README.md` §8.1 ZE-real, and of [RD]'s warning
  about translates `gH`.

**Proposition 5.1: correct under the local ZE-tor.** The equivalence (i) ⇔ (ii) is close to tautological. Under
ZE-tor, every actual cluster in range is blind, so both (i) and (ii) say "every actual cluster in range has
`rho_hat <= theta ln W`". The content lies in the *sufficient* version B'(theta), which quantifies over all
pair-separated exponent profiles and windings. The write-up says this, and B is correctly left open.
* The squarefree remark (`M < 1200 ln 4r + 8` under ZE-tor) is correct: [D] Lemma 3.1 needs
  `W' >= 4M log^+ C'`, which holds for large `R`.
* Remark 5.2 (ZE alone gives `rho_hat < 300 ln 4r` on all norms) is correct. `J_0` is nondegenerate in [D]'s sense
  by Lemma 1.0 and Theorem 3.1(c).

### 7.2 ZE-col is false for thin families; Theorem 5.4's hypothesis must be restricted (major)

ZE-col(F) (for a profile): if `L ⊄ Z_F`, then no `delta' ∈ L` with `|delta'| <= Theta` that is nondegenerate
relative to `L` lies in `Z_F`. Theorem 5.4 assumes ZE-col(univ) **for every profile in the regime of Theorem 4.1**.

**Proposition R (ZE-col fails for a one-measure family; exact, `rc_onerel.py`).** Take on `Lambda = Z^3`

```text
mu = 8[0,0,0] - 8[1,0,0] - 3[0,1,0] - 7[0,0,1] + [1,1,0] + 2[2,0,1] + 3[1,2,-1] + 4[0,-1,2]      (ord mu = 2),
```

and let `F` consist of `mu` and its translates. Then `Z(ghat_mu)` contains nondegenerate points arbitrarily close to `0`.

*Proof.*
* *Expansion.* `ghat_mu(s) = -Q(s)/2 - i C(s)/6 + O(|s|^4)`, where
  * `Q = 4s1^2 + 14s1s2 + 2s1s3 + 14s2^2 - 28s2s3 + 14s3^2`;
  * `C = sum mu_y <y,s>^3`.
* *A common zero of the leading forms.* In the chart `s1 = 1`, the resultant of `Q` and `C` in `s3` is
  `288 (5488X^4 + 5504X^3 - 4032X^2 + 256X + 497)`.
  * The quartic is irreducible over `Q` and has two real roots.
  * At a real root `x`, the common zero is `u = (1, x, y(x))` with
    `y = -(76832x^3 + 92834x^2 - 40481x)/24777 + 6361/49554`.
  * Both `Q` and `C` vanish there identically modulo the quartic.
* *Transversality.* The Jacobian determinant reduces to a nonzero cubic in `x`, so it is nonzero at every root.
* *A curve of zeros.* `Re ghat` is even and `Im ghat` is odd in `s`. So
  `G(t,v) = (Re ghat(tv)/t^2, Im ghat(tv)/t^3) = (-Q(v)/2 + O(t^2), -C(v)/6 + O(t^2))` is real-analytic near
  `(0,u)`. The implicit function theorem on `T_u S^2` gives zeros `s(t) = t v(t)` with `v(0) = u/|u|`.
* *Nondegeneracy.* `(1, x, y)` is `Q`-linearly independent: `y ∉ Q + Qx`, because its `x^2`, `x^3` coefficients are
  nonzero and `[Q(x):Q] = 4`. Hence for each nonzero `y ∈ Z^3`, the function `t -> <y, s(t)>` is analytic and not
  identically zero. So `s(t)` is nondegenerate for all but countably many `t`. QED

The 60-digit Newton continuation in `rc_onerel_out.txt` follows both curves to `|ghat| ~ 1e-60` at
`t = 1e-2, ..., 1e-8`. This is an illustration only; the proof is the exact algebra above.

**Consequence.**
* ZE-col(F) is a statement about the *thickness* of `F|_L`. A family of codimension `< rho_hat` on `L` generically
  violates it.
* [S] never argues that `univ|_L` is thick whenever `L ⊄ Z_univ`. A profile that is visible through a sparse set
  of admissible measures (for example near the visibility threshold) would make the hypothesis of Theorem 5.4
  false, and the theorem vacuous.
* The same applies to Remark 5.5 when only the Siegel vector, rather than all of `C_3`, is admissible.

**Fix (all that Theorem 5.4's proof uses).** Let **ZE-col\*(univ)** assert the ZE-col(univ) conclusion only for
the actual position `delta`, and only for profiles that are *robustly visible*, in one of two ways:
* (V1) `rho_hat >= 300 ln 4r`. Visibility then comes from `V_H(J_0)`, of dimension `>= |A_{J_0}|/2`, injectively
  pushed forward.
* (V2) `M >= (rho_hat+1)(rho_hat+2)` and `4 ln|c|_1 < W - 12 ln C' + 4 ln 6` for a basis of `C_3 ∩ Z^M`. Then all
  `M - E'` (at least `E`) linear-form relations are admissible.

*Proof that ZE-col\* suffices.*
* If (V1) or (V2) held, the cluster would be visible (Theorem 4.1), and ZE-col\* would be contradicted by
  `delta ∈ Z_univ` (Lemma 1.2) together with Lemma 1.0.
* So `rho_hat < 300 ln 4r`. Then either `M < (rho_hat+1)(rho_hat+2)`, or (V2) fails only through the height
  condition.
* The height condition holds for `W >= K_0 E ln(MW)`, which is `O((log W)^3) << W` here (Remark 5.5's
  Bombieri-Vaaler bound). This holds for `R >= R_1(C)`.

### 7.3 Theorem 5.4: normalisation and constants (minor)

* *Normalisation.* Theorem 5.4 is stated for normalised clusters with "`R >= R_0(C)`". After normalisation the
  radius `R' = R/|g|` can be small even when the original `R` is huge. For `R' < R_0(C)`, the honest statement is
  `M <= r_2(R'^2) <= K(C)`.
  * The conclusion for actual circles is `M <= max(K(C), (300 ln(4W'/ln 5) + 2)^2)`, with `W' <= 2 log R`. That is
    still `O_C((loglog R)^2)`.
* *Hypothesis.* With ZE-col replaced by ZE-col\* (§7.2), Theorem 5.4 is a correct implication.
* *Comparison with [D] Theorem 5.7.* The comparison "covers all norms, needs no rank lemma" is fair. "Its
  hypothesis is different, not weaker" is also fair.

**Lemma 5.3(iii): correct.** `Z_univ subset Z_P` (whole cluster); nondegenerate relative to `R^M` means
nondegenerate; and `Theta <= 1/(4D)` in the regime.

## 8. Computations (re-done independently)

All scripts are my own: exact rational ranks, or ranks mod `p = 1000003`. Full rank mod `p` certifies full rank
over `Q`; nothing here relies on a rank *deficiency* mod `p`.

### 8.1 Theorem 4.4 (exact): see §6. No violation in 950 random sets and 50 adversarial hill-climbs; the equality cases are reproduced.

### 8.2 §6.4 "honest small models" (exact, `rc_honest.py`)

*What the author tested.* `honest_small.py` takes, for each code, only the **first 6 subsets in lexicographic
order**, with **base point = last point**, at `n = 4, 5`. That makes 24 systems per profile and 120 in total.
* The text says "`n = 4..6`" and concludes that "the combined system of an honest profile is empty". That is a
  generalisation from 24 of up to `C(M,n)·n·3` systems per profile, and `n = 6` was not run.

*My run.* I used my own pair-separated codes (minimum distance `>= r/2`, equal weights), *random* `(J, base)`
including all base points, and `n = 4, 5, 6`.
* `M = 6..12`: 60 + 40 + 30 systems per code. `n = 6`: 12 systems at `M = 8` and at `M = 12`.
* `V(A, floor(LP)+1) = 0` in every case (full rank mod `p`; `rc_honest_out.txt`).

So the author's finding is reproduced and extended, as **evidence for these slab-boxes at these sizes**.

*What the finding does not show.*
* "No relation is forced, `Z_comb` is the whole torus" is true **by definition** at these sizes, since
  `|J| - 1 < 300 ln 4r`. The computation does not bear on it.
* The universal family of Lemma 1.1 (other supports, other `T, w`) was not tested.

### 8.3 §6.5 tangent cones decided exactly; the evidence is structurally uninformative (`rc_tc_exact.py`)

For each of the 280 `(J, base)` systems of each §6.5 model, I rebuilt `A`: `n = 4`, `T = 3`, median slab cut, and
`m` = the largest integer with `4 binom(m+2, 3) <= |A|`, which gives `m = 6`. I then computed the rank of the
evaluation matrix of *all* monomials of degree `<= m` on `A`.

| model | systems | `m` | `|A|` | no degree-`<= 6` polynomial vanishes on `A` |
| --- | --- | --- | --- | --- |
| skew double core, `q = 7` | 280 | 6 | 225 | 280/280 |
| round-5 Paley-I, `b = 8`, `q = 7` | 280 | 6 | 225-239 | 280/280 |

*Consequence (exact).* Suppose a form of degree `m` agreed on `A` with a polynomial of degree `< m`. Their
difference would be a nonzero polynomial of degree `m` vanishing on `A`. So the single tangent cone is **exactly
`{0}`** in all 560 systems, and so is the combined one.
* The BFGS statements ("no zeros; `min tau = 9.6e-4`", from 12 sampled systems) are therefore true. They hold
  for a trivial reason, and they did not need floating point.

*Why this carries no information about the regime.*
* Here `m = 6 < 2T + 1 = 7`.
* In the regime, `T = n` and `m >= 0.7 n(n-1) >> 2n + 1`. Then `prod_{j=-T}^{T} (k_x - j)` (degree `2T + 1 <= m`)
  vanishes on `A`.
* So the tangent cone always contains the sparse directions `e_x` ([D] Proposition 5.1). It also contains the
  column directions `beta_p` whenever `e_p w < m`: `<k, beta_p>` takes at most `e_p w + 1` values on `A`.
* A model with `m <= 2T` is in a different structural class from the regime. Its trivial tangent cone is no
  evidence for ZE.
* Likewise, "100/100 fake and random positions certified off every single base locus" is automatic: every `Z_P`
  is Lebesgue-null ([D] Proposition 4.2), and [RD] §4 made this point.

On the word "certified":
* `certify_off` uses a floating-point `c = (I - P) f`. Mathematically the exact `c` has leading term
  `||c||^2/m! > 0`, so the test is a well-conditioned float evaluation of an exact criterion (`tau ~ 1e-3`).
* It is not an interval certificate. The real span of `V_H` was not checked to equal `V(A,m)` at these sizes.
  [RD] checked it only at other sizes.
* The label "evidence" in the section header is right. "Certified" in the bullet is too strong.

### 8.4 Where the universal family starts (exact binomials, `rc_threshold.py`)

Take the [D] slab-box (`T = n`, `w = n sqrt(8(n-1) ln 4r)`), with exact `m` (largest `m` with
`8 binom(m+n-2,n-1) <= (2n+1)^{n-1}`):

| `r` | first `n` with `m > 2w` (so `V(A,m) != 0`, `LP <= 2w < m`) | `m >= 4w` ([D]'s margin) |
| --- | --- | --- |
| `n` | **440**, and every `n` in `[440, 2400)` | 2137 |
| `n ln n` | 563 | 2679 ([D]) |
| `n^2` | 881 | 4278 ([D]) |

At `r = n`, with `W = W_min(n) = 3n ln n - 14n` (the `two_thirds` lower bound, `C = 1`), the [D] integer generators
satisfy [S]'s admissibility (1.2) for **every `n` in `[468, 2400)`**. Here `LP <= 2w` and
`|c|_1 <= U^3 n^{2(m-1)}`.
* Hence the task summary's "Below `n ~ 2137`, all 120 tested pair-separated sub-cluster systems are empty" is
  wrong in its quantifier. "Below the regime the honest combined system is empty" (§0 item 8) holds only in the
  vacuous sense that `Z_comb` is defined from `|J| - 1 >= 300 ln 4r` (about 2717 at `r = n = 2137`).
* For the universal family that [S] itself introduces (Lemma 1.1), honest relations exist from `n ≈ 470`.

### 8.5 Other numerical claims

* `lowrank_profiles.py`:
  * The header says "exact LP ... verified with fractions". The code contains no fraction verification; it is a
    float HiGHS LP over a finite test set `[-6,6]^rho`.
  * Restricting to a finite test set *over*-estimates the min-norm. So the non-existence conclusions
    (`s = 0.8, 0.667 < 1` at `M = 6, 7`, and by monotonicity in the point set at all larger `M`) remain valid
    modulo float error, with large margins.
  * It is evidence, not exact. Minor.
* §6.1 (cube profile): the positive dimensions and memberships are mod `p`, and labelled so. Not re-run; not
  used in any proof.
* §6.6: correctly labelled inconclusive.

## 9. Issues list

| # | severity | where | issue | fix |
| --- | --- | --- | --- | --- |
| 1 | major | §5.2 ZE-col, Thm 5.4 | ZE-col(univ) is assumed for *every* profile and *every* nondegenerate point of `L`. ZE-col is false for thin families (Prop. R, §7.2, exact), and thickness of `univ|_L` is never shown. | Assume ZE-col\* only for robustly visible profiles (V1)/(V2), and only at the actual `delta` (§7.2). Theorem 5.4 then holds as stated. |
| 2 | major | §5.2 ZE-tor, Lemma 5.3(i)-(ii), Prop 5.1 | At the fixed radius `1/(4D)`, subtori have extra sheets and intersections have torsion cosets. "Near `e` a closed subgroup = its identity component" is false at that radius, and "`L subset K`" needs the identity sheet. | Define ZE-tor locally (lifted ball, rational subspaces through `0` containing `1`); then (i), (ii), Prop 5.1 hold. |
| 3 | major (interpretive) | §0 item 3, §3 Consequence, table; task summary | "Nondegenerate sub-clusters ... add nothing" is unsupported and contradicts Remark (i) of §1.2. In the full-rank case (`rho_hat = M - 1`, to which [D] Theorem 5.7 reduces) the push-forward is the identity, and the intersection over nondegenerate `J` is the round's question itself. "Reduces exactly to ZE-col" is then a restatement. | Say: Thm 3.1 isolates the extra content of *degenerate* sub-clusters. For nondegenerate ones it gives nothing either way, and the full-rank intersection question is untouched. |
| 4 | major (claim) | §0 item 8, §6.4, task summary | "Below `n ~ 2137` the honest combined system is empty" was tested only at `n = 4, 5` (first 6 subsets, one base point). Exactly, slab-boxes are LP-supercritical for every `n` in `[440, 2400)` and Lemma-1.1-admissible for every `n` in `[468, 2400)` (`r = n`). | Restrict the claim to "the tested slab-box systems at `n <= 6`". Note that `Z_comb` is trivial below `300 ln 4r` by definition, and that `univ` is nonempty from `n ≈ 470`. |
| 5 | moderate | §6.5, summary | "Fakes certified off every locus" is automatic (null loci). "No tangent-cone zeros" is exactly true in all 560 systems (§8.3), but only because `m <= 2T`; in the regime `m >> 2T + 1` and the tangent cone always contains the sparse and column directions. | Label §6.5 as structurally uninformative about ZE. Replace BFGS with the exact criterion. |
| 6 | minor | Prop 3.2 "Reading", §0 item 3 | "No stronger than characters" holds only at leading order. Inside the `O(log W)`-per-unit window products can exclude profiles that characters do not. That window is where the *constant* of `log R/loglog R` lives. | State as a leading-order and heuristic order-of-growth comparison. |
| 7 | minor | Thm 5.4, Cor 4.2 | `R >= R_0(C)` must be the normalised radius; small normalised radii need `M <= K(C)`. | Add `max(K(C), ·)`. |
| 8 | minor | Prop 5.1 | The equivalence (i) ⇔ (ii) is near-tautological under ZE-tor; the content is in B'. | Say so. |
| 9 | minor | task summary item 2 | "Its fibres are finite": Prop 2.2 allows whole-circle fibres. | Add "or full". |
| 10 | minor | Cor 4.5 Example | Primes in `[P, P^{1.9}]` give `gamma > 1/2` but break pair separation (pair sums can drop to `≈ 0.69`). Integral dependence of `kappa` on the columns is assumed. | Note that only near-equal weights give both. Harmless (no admissible sub-cluster). |
| 11 | minor | §6.2 `lowrank_profiles.py`, §6.5 | The "verified with fractions" header is false (float LP). "Certified" in §6.5 is a float criterion, not an interval certificate. | Wording. |

No theorem is wrong as proved, apart from the quantifier and definition defects in items 1-2, which have local
fixes. No numerical evidence is used inside a proof. The places where evidence or a restatement is *presented as*
a finding (items 3-5) are the main corrections.

## 10. Corrected statement

**No growth improvement and no uniform bound.** Intersecting the forced-divisibility rigidity relations over
sub-clusters does not decide whether actual positions avoid the base locus. What is proved, for a normalised
cluster (`N` odd, all primes split, `W = log N`, arc constant `C`, `C' = max(C,1)`), is the following.

1. **(Lemmas 1.0-1.2.)**
   * `delta` lies in the column space `L = span(1, columns of (a_x, kappa_x))` and is nondegenerate relative to `L`.
   * Every integer `c` on `Z^M_0` satisfying (1.2) with `m = ord_L(c)` has `sum c_k e^{i<k,delta>} = 0`. All [D]
     sub-cluster slab-box systems (`|J| - 1 >= 300 ln 4r`, `|J| >= 64`, every base point) are instances.
   * Sub-cluster systems are genuinely new relations (§1.2 Remark (i)).
2. **(Props 2.1-2.2.)**
   * For each admissible size `n`, the combined locus has dimension `<= n - 2` for every `M >= n`.
   * Each fibre over `n - 1` coordinates is either the full circle or has `<= 2T` points.
   * Neither statement bears on one actual point.
3. **(Thm 3.1, Prop 3.2.)**
   * On `L`, each relation is the trigonometric polynomial of the push-forward `hat pi_* c`, with `ord_L(c) = ord(hat pi_* c)`.
   * For an affinely independent `J`, `hat pi_*` is injective on `Z^J_0`.
   * This isolates the extra content of *degenerate* sub-clusters. It says nothing about intersections over
     nondegenerate ones; in the full-rank case (`rho_hat = M - 1`, the case [D] Theorem 5.7 reduces to) that intersection is the open
     question itself.
   * An admissible product push-forward `mu_alpha` can exist only if
     `sum_i alpha_i (X_{v_i} - W/2) < -2 ln|c|_1 - 2|alpha| ln(DC'/2) + 2 ln|alpha|!`, and it cannot exist for an
     actual cluster. So product certificates act only within `O(log W)` per unit of `|alpha|` of the character
     threshold. That they are no stronger than characters at the level of growth order is heuristic.
4. **(Thm 4.1, Cor 4.2; unconditional.)** Assume `W >= 16 ln M + 8 log^+ C + 14` and
   `W > 8 ln M + 8 ln(2W) + 12 ln C'`.
   * If `rho_hat >= 300 ln 4r`, then `L ⊄ Z_comb`.
   * If `M >= (rho_hat+1)(rho_hat+2)`, then some `c` supported on `{e_x - e_b}` with `|c|_1 <= 12 M^2 W^2` and
     `ord_L(c) >= 3` gives `sum_x c_x z_x = 0` with `hat pi_* c != 0`; hence `L ⊄ Z_univ`.
   * A cluster with `L subset Z_univ` therefore has `M < (rho_hat+1)(rho_hat+2) < (300 ln(4W/ln 5) + 2)^2`.
   * Nothing shows that actual clusters are blind.
5. **(Thm 4.4; exact tests agree.)**
   * `reg(S) <= E_eps h_S(<eps,.>)` for every finite `S subset Z^rho`, with equality on boxes and `l^1`-balls.
   * It strengthens Theorem B of `threshold/` and is coordinate-dependent.
   * Corollary 4.5: cross-polytope column geometries with `gamma > 1/2` are blind to all slab-box sub-cluster
     systems. The known examples have no admissible sub-cluster, so this is a mechanism, not a statement about
     actual clusters.
6. **(Prop 5.1, under local ZE-tor.)** Assume, for every single admissible slab-box system,
   `Z_{P_J} ∩ B = ∪ (V_i ∩ B)` in the lifted ball `B` of radius `1/(4D)`, with `V_i` rational subspaces
   containing `1`.
   * Then every actual cluster with `M >= 300 ln 4r + 1` is blind to `Z_comb`, so `rho_hat < 300 ln 4r`.
   * Whether the intersection forces `rho_hat <= theta ln W` is the open combinatorial question B(theta); its
     sufficient form B'(theta) is pure linear algebra.
7. **(Thm 5.4, corrected hypothesis.)** Assume **ZE-col\***: for every actual cluster whose profile is robustly
   visible, its position `delta` does not lie in `Z_univ`. Robustly visible means (V1) `rho_hat >= 300 ln 4r`, or
   (V2) `M >= (rho_hat+1)(rho_hat+2)` with a basis of the degree-2 kernel `C_3` admissible. Then for `R >= R_1(C)`:

   ```text
   M  <=  max( K(C), (300 ln(4W/ln 5) + 2)^2 )  =  O_C((loglog R)^2)
   ```

   This holds on all norms. ZE-col in the original form (all profiles, all nondegenerate points of `L`) is false
   for thin relation families (Proposition R, §7.2), so the restriction is necessary for the hypothesis to be
   credible.

**Computations.**
* *Exact.*
  * Theorem 4.4: 950 random sets and adversarial searches, no violation.
  * In the §6.5 models, every single tangent cone is exactly `{0}` (560/560 systems). This is because
    `m = 6 <= 2T`, so it carries no information about ZE in the regime.
  * The tested honest slab-box systems at `n <= 6` have `V(A, floor(LP)+1) = 0` (the author's 120, plus my 544
    random `(J, base)` samples).
  * The [D] slab-box is LP-supercritical for every `n` in `[440, 2400)` and Lemma-1.1-admissible for every `n` in
    `[468, 2400)` (`r = n`). So "below `n ~ 2137` the honest system is empty" is false for the universal family.
* *Evidence.* Fake and random positions lie off the loci, which is automatic since the loci are null. The
  non-honest supercritical model (§6.6) is inconclusive.

## 11. Scripts (`round11/referee_subcluster_checks/`)

| script | what it does | status |
| --- | --- | --- |
| `rc_lib.py`, `rc_thm44.py` | exact `reg(S)` (fraction-free rank), `Omega'`, `MW`; 950 random sets, hill-climbing, equality cases | exact |
| `rc_honest.py` | own pair-separated codes; random `(J, base)`; `n = 4, 5, 6`; `V(A, floor(LP)+1) = 0` | exact (full rank mod `p`) |
| `rc_tc_exact.py` | §6.5 models: no polynomial of degree `<= m` vanishes on `A`, hence the tangent cone is `{0}` in 560/560 systems | exact (full rank mod `p`) |
| `rc_threshold.py` | exact-binomial thresholds for `m > 2w`, Lemma 1.1 admissibility, `m >= 4w` | exact `m`; real-number comparison of explicit bounds |
| `rc_onerel.py` | Proposition R: a one-measure family on `Z^3` with nondegenerate zeros arbitrarily near `e` | exact algebra (sympy) plus a 60-digit illustration |
| `notes.md` | referee scratch notes | - |
