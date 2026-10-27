# Referee report on round16/barrier.md (adversarial)

Setting checked throughout: T = 1, w = sqrt(8 d ln 4r), j = floor(w/2), m' = ceil(4w), regime d = n-1 >= 316 ln 4r,
barrier scale W = r tau, tau = K ln n; fake profile sdc_q (n = q+1, r = 2q, mu = 1, k = e_1). Scripts (my own code,
nothing imported from barrier_checks/) are in `referee_barrier_checks/`; `x.py` writes `x_out.txt`.

## 0. Verdict

**The headline survives: no growth improvement, no uniform bound, and route 5 is neither capped nor unblocked.** The
proofs of Thm 2.3, Prop 2.5, Prop 3.1, Thm 3.2, Cor 3.3 and Lemma 3.4 are correct up to minor, repairable slips. Three
claims in §0 and the Status block are **logically overstated**, and they must be corrected before the write-up enters
the record (§8 gives the corrected statement):

1. "AC is strictly stronger than MIN^pair" and "a proof of AC would decide the barrier question for all fake boxes":
   both are false as stated. AC controls kappa_V, while MIN^pair and the fake barrier are about kappa_H <= kappa_V. So
   AC implies neither of them unless V_H = V(A,m'), which is open. The author's own Cor 2.2(3) says this, so the
   summary contradicts the body.
2. "Collision strata never give barrier points for fakes" (§0 item 3, §5 item 2, "never the issue"): Prop 2.4 proves
   only that the **Lemma 2.0 transversal ceiling cannot certify** kappa <= eps_0 near a stratum with a light relation.
   It does not exclude near-zeros or exact zeros at those distances that come from other mechanisms.
3. "The circle case reduces to one explicit algebraic question": this holds only in the "yes" direction. An exact
   "no" (circle not in TC) does not settle the barrier test, which is quantitative at scale eps_0 = e^{-Theta(q log q)}.

No fake-pattern-box counterexample (barrier point) was found. My adversarial tests (§6) support the author's data and
hypotheses: exact AGL(1,q) symmetry, sigma_xy = 1, and a global MILP minimum for non-pair slacks.

## 1. Prop 2.1 / Cor 2.2 (AC is box-free)

The proof is correct and trivial: S = {|Q| <= delta} in one direction, a best approximant in the other, then
finite-set LP duality. Two problems remain.

* **Major (logic).** Cor 2.2(2) says AC on the gap regions "would decide the barrier question of fake Thm 2.1 for all
  fake pattern boxes at once". The barrier question is whether fakes have kappa_H(hat u) <= eps_0, where kappa_H uses
  Step 0's V_H. AC gives kappa_V^{(A)} > eps_0, and since kappa_H <= kappa_V this leaves kappa_H <= eps_0 possible.
  AC therefore decides only the H = infinity version. For the same reason "AC is strictly stronger than MIN^pair" is
  false: AC does not imply MIN^pair at all (the author's Cor 2.2(3)), and no example is given of MIN^pair holding
  while AC fails. Correct phrasing: *AC implies the H = infinity form of MIN^pair; it implies MIN^pair only together
  with V_H = V(A,m') on the relevant set.*
* **Scale remark (added; this is not an error).** For every u there is the Chebyshev witness
  Q = 2^{1-m'} s^{m'} T_{m'}(<k',u>/s), with s the median of |<k',u>| over T (about 0.55 |u|_2). It has top form
  l_u^{m'} and satisfies |Q| <= 2 (0.275 sqrt d)^{m'} on half the cube. So AC fails for **every** u unless
  W/4 >= m' ln(sqrt d/0.275) + O(log d), i.e. unless roughly W >= 2 m' ln d + 5.2 m'. At the relevant scale
  W >~ n log n (the steelman's t = e^{-Theta(n log n)}) this is satisfied, so AC is not trivially false. But this
  scale condition must be stated as a hypothesis of AC; barrier.md never states the scale of W in AC.

## 2. Thm 2.3 (only light characters need Diophantine control)

I re-derived the proof against fake.md §2.1.
* Class (1'): the constant 12 * 1.01 * e^{1/16} = 12.9 < 13 is correct. The second term of (2.2') gives
  eta ||lambda||_1 <= |lambda.tilde u|/2.
* Class (3'): the bound 24 n (1 + 1/(Delta ||lambda||_1)) e^{-w/2}/eta is fake.md's. Substituting gives
  <= 384 e^{1/8} n^{3/2}/(theta C ell_0) e^{-tau(sigma - r/2)/2} <= 450 (...), which is correct. The group ratio
  e^{2n/r - tau/2} <= 1/2 holds for r >= 2n-2, n >= 11, tau >= 4.
* **Minor slip.** The first heavy group is bounded by exactly 450/1800 = 1/4. The extra factor e^{-tau/2} claimed in
  the text does not follow from the definition of sigma*. The geometric sum is still <= 1/2, so the conclusion stands.
* The thresholds are O(log n): class 2 needs only w(v(lambda)) >= tau (v(lambda) is a nonzero integer vector by (AI)),
  and the prime window needs tau >~ 2 ln r.

**Verdict: correct.** "Heavy relations are free" is genuinely proved.

## 3. Prop 2.4 and Prop 2.5 (rational strata)

**Prop 2.4.** The algebra is correct (log of c_1^{m'} eta^o > eps_0, W >= r tau). The **interpretation is overstated**:
* The ceiling of Lemma 2.0 is an *upper* bound. Showing that it exceeds eps_0 shows only that this certificate fails;
  it does not show kappa_H > eps_0. The order-o transversal coefficient a_o(y,v) can itself be small or vanish for
  special (y,v), which is a higher-order stratum, and other real components of TC_H can pass through the collision
  neighbourhood.
* Prop 2.4 assumes u* satisfies (2.2'), which is the *sufficient* condition of Thm 2.3. To speak about *all*
  gap-admissible fakes, one should use the necessary pair inequality: gap-admissibility is equivalent to the
  character inequalities (outside.md Lemma 0.1, quoted in fake.md §1). That gives
  |hat u_x - hat u_y| >~ (2/C) e^{-tau sigma_xy/2}, the same order with a different constant. The same remark applies
  to row 10 of §4 (sum-form directions): the exclusion is valid only through the necessary form.
* Corrected statement: *at the barrier scale, for sdc_q (sigma_xy = 1), no transversal-vanishing certificate along a
  stratum carrying a light relation lambda_0 with 2 o sigma_{lambda_0} < r - 4 ln d/tau can certify a barrier point.*

**Prop 2.5.** The construction is correct: y_1 avoids pair zeros (Schwartz–Zippel count); y_2 uses base-B digits;
v uses base B' with the l1 bound (|lambda_x| <= 2r+2 > B'/2 is harmless, because the total l1 bound already prevents
carries). Light lambda in L_T are integral for sdc, since all fractional characters are irrelevant (fake.md §1).
* **Minor gap in tau_L.** eta = 26 e^{-tau sigma_min/2} (2r+3)^{d-1}/C contains the factor (2r+3)^{d-1}. Both the pair
  check (1/(3bTh_L) - 2 eta) and the light-non-pair check (epsilon/(2bTh_L) - eta(2r+2)) need
  tau c_0 r/2 >= (d-1) ln(2r+3) + ..., i.e. an extra term (2/(c_0 r))(d-1) ln(2r+3) + O(1), which is about
  ln(2r+3)/c_0. tau_L omits it, and it matters when b = dim L is small. It is still O(log n), so the conclusion is
  unaffected.
* Prop 2.5 needs rho_1 >= 1/2 (to have ||lambda||_1 <= 2r+2 for light lambda). This holds for sdc and pal; state it.

## 4. §3 (products, cyclic symmetry, circles)

* **Prop 3.1** (additivity of nu on products; nu >= 2 min(j, s_b)): proof correct. J(B) = (k_x^3, k_y) is a correct
  Gröbner computation. **Minor error:** "m' <= 8j+4 (regime)" is false. In the sdc regime m' > 8j+4 for 4446 primes
  q < 2*10^5 (first q = 3271: j = 257, m' = 2064 > 2060). The correct bound is m' <= 8j+8, with max m' - 8j = 8
  attained (`regime_arith_out.txt`). "At most 4 blocks of size >= j" still holds for j >= 4.
* **Thm 3.2**: proof correct. **Strengthening:** |f_c|^2 has only *even* frequencies 2l', |l'| <= m. For odd N the
  orbit average equals the mean whenever **m < N**, not only when 2m < N. The same relaxation applies to the
  eigenvalue separation in Cor 3.3. Numerically (`thm32_check_out.txt`, sdc7 with L = 2, N = 7) the orbit mean of
  |Phi|^2 over ||f||_2^2 is 1.000000 for m = 4, 5, 6 (where 2m >= 7) and 0.859 at m = 7 = N. This matters because
  **the author's claim "2m' < q holds in the regime" is false at the regime start:** the first regime sdc prime is
  q = 3251 with m' = 2057, and 2m' < q first holds at q = 5479. With the relaxation, m' < q holds throughout and
  Thm 3.5 covers all regime q.
* **Cor 3.3**: correct. **Strengthening** (`sdc_sym_out.txt`, q <= 83, and a two-line proof): sdc_q is invariant under
  the full affine group AGL(1,q). x -> tx with t a residue maps P_a -> P_{ta}; t = -1 maps P_a -> P'_{-a}. Hence all
  Fourier circles C_a are images of one another under box automorphisms, and kappa_H is *exactly* equal on all of
  them (quantitatively, not only for exact zeros as with Galois). AGL(1,q) is 2-transitive on F_q, which explains
  sigma_xy = const.
* **Lemma 3.4**: correct. The relevant lambda are integral, beta_lambda = 0 iff lambda|F_q is constant, and the norm
  bound and the measure count are as stated. Note that the fractional character (2/(q+1))(1_{F_q} - q e_inf), the
  generator of L_T/Z^n_0, has sigma = r/2 exactly. It vanishes identically on every circle, and it is harmless only
  because it is irrelevant for k = e_1. Any other choice of windings that makes it relevant would make it light with
  an exact zero on the circles. Worth stating.
* **Thm 3.5**: (i)/(ii) and the "Consequence" are correct (modulo H_light). **Major (interpretation):** the paragraph
  "turns the barrier question ... into one explicit algebraic question" and §8 item 1 ("Exactly: ...") overstate. A
  "yes" (circle in TC_V, hence in TC_H) caps route 5. A "no" leaves the quantitative test
  inf_{C_1} kappa_H <= eps_0/(2 sqrt(m'+1)) completely open. Exact non-membership plus a Liouville bound at
  theta = 0 (u in Q(zeta+zeta^{-1}), algebraic integer after scaling by 2^{m'}) gives only
  kappa_H >= exp(-(q/2)(m' ln d + q ln 3 + ln H) - O(m')) = e^{-Theta(q^2)}. The barrier needs e^{-Theta(q log q)}.
  Moreover, all exact circle data (sdc7, sdc11 at L = 2, m = 8; sdc19 at m = 5, 6) are exact-zero statements at toy
  (L, m), not regime (j, m').
* **Sylvester aside** (end of §3.3): "MIN^pair would fail for Sylvester-type data at every W" needs
  eps_xy >= 2 ln(n^2/c) for all pairs. Sylvester data have zero pair margin, so this is unsupported. It should be
  removed or conditioned.

## 5. §4–§5 (table, separation, (P1)/(P2))

* Row 7(a)/(b): proofs check (the RNC base-locus count gives e <= 8; the coordinate-point RNC is
  t -> (c_x/(t - a_x)), and the Möbius map forces at most 2 values of a_x).
* Row 10 needs the necessary form of the light-character inequality (see §3).
* §5 item 3 is stated as "iff", but the third alternative is "a kind not yet identified", so the equivalence carries no
  content. It should be stated as a sufficient trichotomy.
* (P1)/(P2) are correctly labelled as unproved targets.

## 6. Independent computations (all mine)

| check | result | file |
| --- | --- | --- |
| sdc_q built from the Paley-I definition: translation and **all** F_q^* multiplications preserve the pattern set | true, q = 7..83 (q ≡ 3 mod 4) | sdc_sym_out.txt |
| sigma_xy for all pairs incl. base | exactly 1, q = 7..83 (author: q <= 43) | sdc_sym_out.txt |
| exhaustive min sigma/r, norm 4 (e+e-e-e) | 0.2143, 0.2273, 0.1842, 0.1957 (q = 7, 11, 19, 23); matches author | hlight_ref_out.txt |
| exhaustive min sigma/r, norm 6 (not computed by author) | 1^3 vs -1^3: 0.357, 0.227, 0.289, 0.283; other types >= 0.54 | hlight_ref_out.txt |
| exhaustive norm 8 (q = 7, 11) | min 0.357, 0.409 | hlight_ref_out.txt |
| **global MILP** min sigma over integral non-pair lambda, \|lambda_x\| <= 3, support >= 3 | optimum = the 4-term type, proved optimal by HiGHS: sigma/r = 0.2143, 0.2273, 0.1842, 0.1957, 0.2097, 0.1977 (q = 7..43) | hlight_milp_out.txt |
| Thm 3.2 averaging identity | holds for m < N (even when 2m >= N), fails at m = N | thm32_check_out.txt |
| sdc7, L = 2, m = 8 | \|A\| = 757, dim V = 1 (exact, fmpz nullspace); c translation-invariant and even; Phi_c = -7717.5 constant on C_1, C_2, C_3 (50 digits): circles not in TC | sdc7_octic_out.txt |
| regime arithmetic | first regime sdc q = 3251; 2m' < q only from q = 5479; m' <= 8j+8 sharp | regime_arith_out.txt |

So H_light(c_0) holds empirically with c_0 = 0.184 for q <= 43 (B = 3), which supports "o >= 3" in Prop 2.5 as
**evidence**. I did not re-run sdc11 at m = 8 (a 31824 x 9373 nullspace), sdc19 at m = 5, 6, or the M8 searches;
those remain the author's (mod-p, respectively floating-point) data.

**Counterexample search.** I tried four further sources of barrier points:
* symmetric directions of AGL(1,q): multiplicative-character circles have u_0 = u_inf, a base collision;
* fixed spaces of subgroups: constant on cosets, so collisions;
* the integer ramp u_x = x + 1: collision-free, but its 4-term relations are light and it has no visible cheapness;
* the half-cube freedom in AC: cheap integer v plus a small perturbation, which needs eps_xy >= W/(2o) for most
  pairs, and that is excluded by the record's separation radius.

None of these yields a barrier point. No counterexample found.

## 7. Literature (checked by search)

* Meka–Nguyen–Vu, Theory of Computing 12 (2016) art. 11: confirmed.
* Kane, Ann. Probab. 45(3) (2017) 1612–1679, arXiv 1204.0543: confirmed. The journal title is "...poorly
  anticoncentrated polynomials of Gaussians and applications to the study of polynomial threshold functions"; the
  number of structural polynomials is bounded in terms of the degree. This is consistent with the author's reading.
* Barber, arXiv 2104.15057: the title is "Small sums of **five** roots of unity" (BLMS 2023).
* arXiv 2607.06098 ("Small sums of roots of unity"): exists.
* Konyagin–Lev and the "k <= 4" statement: not re-checked; label [M].
* Alon–Füredi min-support 3^{s - m/2} on {-1,0,1}^s: correct as an [M] consequence.

## 8. Corrected statement

At T = 1, regime parameters, barrier scale W = r tau, tau = K ln n, and for sdc_q (r = 2q, k = e_1):
(a) AC(u) is equivalent to "no half-cube S has kappa_V^{(S)}(u) <= eps_0". It is meaningful only for
    W >~ 2m' ln d + 5m'. It implies the H = infinity form of MIN^pair; it implies neither MIN^pair nor a decision of
    the kappa_H barrier unless V_H = V(A,m').
(b) Thm 2.3 holds as stated: only light characters need (2.2').
(c) Transversal-vanishing certificates along strata with a light relation lambda_0 cannot certify barrier points when
    2 o sigma_{lambda_0} < r - 4 ln d/tau. For sdc_q this covers every collision stratum for q > m' + O(1). Near-zeros
    near collision strata from other mechanisms are **not** excluded.
(d) Prop 2.5 holds, with tau_L enlarged by about ln(2r+3)/c_0 and with rho_1 >= 1/2 assumed.
(e) Thm 3.2 and Cor 3.3 hold for m < N with N odd. sdc_q is AGL(1,q)-invariant, so kappa_H is identical on all
    circles C_a. Under H_light, if inf_{C_1} kappa_H <= eps_0/(2 sqrt(m'+1) e^{pi^2/(4q)}), route 5 is capped, for
    every regime q (not only q >= 5479). An exact "no" to C_1 ⊂ TC does **not** settle this.
(f) No barrier point is constructed and none is excluded. Growth and uniformity are unchanged.
