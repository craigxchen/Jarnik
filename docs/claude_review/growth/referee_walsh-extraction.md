# Referee report: `walsh-extraction.md` (Theorems A, F, B, B'', C; Proposition P; Sections 7-10)

**Verdict.** The rigorous results survive: Theorems A, F, B, B'', C and Lemmas 7.2-7.3. Their proofs
are correct as written. I re-derived every step by hand, and the item-333 steps that Theorem C
reuses are re-checked below. Exact tests are in `referee_walsh_checks/`.

Several surrounding claims are overstated:

* the "barrier is *exactly* non-Walsh core with almost all rows flipped" characterisation is not proved;
* most of Proposition P is a restatement of items 334/391-393;
* the comparisons with items 276/277 understate what those items already prove;
* Theorem C's growth statement follows from items 278 + 333 without Theorem B.

No general growth improvement and no uniform theorem is claimed, and none is obtained.

## 1. Theorem B (pure Hadamard cores): correct

Steps checked:

* **Pinning.** `H^t H = M I` and balanced columns give `phi_a = (pi/2M) m_a + eta_a` with
  `|eta_a| <= Delta/2`; `theta'` cancels.
* **Two-term composites.** For `mu = lambda_a - lambda_b`:
  * the column image is `e_a - e_b`;
  * `||mu||_1 = 1`, since distinct columns differ in exactly `M/2` rows;
  * `beta = G_a conj(G_b)` is a conjugate-primitive nonunit of odd norm.
* **Lemma 1.1.** `|Y| = |beta| |sin(arg - j pi/2)| <= |beta| e` and `Y != 0`, so
  `V >= 2 log(2/(||mu||_1 Delta))`. This is exactly the reverse of (F.1).
* **Sidon step.**
  * For 4-term combinations of light labels, `||mu||_1 <= 4` and `V < 4 theta_L`.
  * Residues `r_a = m_a mod M` with `r_a + r_b = r_c + r_d` and `{a,b} != {c,d}` would force
    `mu^t k` to be an integer.
  * So `{r_a}` is a Sidon set in `Z/M`, which gives `|B|(|B|-1) <= M-1`.
* **Heavy labels.** At most 8 when `W > 36 log^+(2C)`.
* **Lemma B.1.** Correct. The bound is tight at `gamma_+ = 0`: Walsh-16 with weight on the 8 labels
  `a.e = 1` has every `G_xy <= 0` and `|A| = 8 = M/2`.
* **Corollary B.2.**
  * The orders failing `M/2 > s_M` are exactly 4, 8, 12, 16, 20, 24 (checked to 5000).
  * The 14 smallest split primes give `sum log = 53.17 > 36 log 2 = 24.95`.
  * So `M <= 24` at `C <= 1` holds.
* **Corollary B.3.** Inequality (6) of `inert_prime_cofactor_bound.md`, `F(M) <= b_M log R + binom(M,2) log C`,
  applied with primitive radius `e^{W/2}`, gives `W >= 8F(M)/M - 4(M-1) log^+ C`.
  * Primitive radius is correct: every nonconstant block occurs in both orientations, so the cluster gcd is `g`.
  * With `F(M) >= (1/4+o(1)) M^2 log M`, this gives `M_1(C) = C^{2+o(1)}`.
  * `M_1(C)` should be defined as the threshold beyond which the inequality holds for *all* larger orders.

**Hypothesis missing from the summary.** (B.1) needs `M > 2 s_M`, i.e. `M >= 28`. The summary's
"any pure Hadamard-core cluster satisfies W <= max(...)" drops this. For `M <= 24`,
`M - 2 s_M <= 0` and the formula is meaningless.

**Literal tests** (`ref_pinning_literal.py`). On 32 literal Gaussian cores (Walsh 8/16, Paley 12/20)
with random units and actual split-prime blocks:

* the pinning holds and (1.2) holds for every 2-term and 4-term combination, with the actual width
  `Delta` and recovered winding integers;
* Lemma 1.1's consequence holds at all 2338 residue coincidences;
* a search over 286k order-4 pure cores finds short-arc rectangles at `C ≈ 2.01`, e.g. blocks
  `12+5i, 46+19i, 56+23i`. Their block angles lie within `Delta/2` of `(pi/8)Z`, as the pinning requires.

**Novelty relative to the notes.** No note uses the cross-label residue comparison or a Sidon bound in
`Z/M`. I grepped "Sidon", "residue", "pigeonhole", "winding" and `pi/(2M)`. The grid pinning itself is
`fixed_hadamard_roth_phase_obstruction.md` (5). The comparison of *different* labels through `m_a mod M`
is new and is the paper's best idea. It replaces a degree-dependent single-target estimate (Roth/Matveev)
by an axis estimate of degree one.

But the claim needs qualifying:

* Item 276 (`effective_hadamard_short_character_phase_gap.md`) already proves an effective uniform
  bound for pure cores of **any** normalised Hadamard matrix, Paley included. It uses the integral
  character `lambda = H_0 e_j` plus Matveev, in the unit regimes `C <= sqrt 2` (or `C <= 2` with a common
  unit).
* Item 278 covers pure Walsh for every `C`.

So what is new is:

1. every `C` for non-Walsh cores;
2. explicit small thresholds (`M <= 24` at `C <= 1`, versus `2^82`-size cutoffs);
3. an elementary method.

"Pure Paley passes every integral character inequality, so Theorem B is a genuinely non-integral
obstruction" is true only for the Liouville-level half-height inequality. An integral character combined
with Matveev already excludes pure Paley at `C <= sqrt 2`.

## 2. Theorem B'' (row-deleted cores): correct

Steps checked:

* `|A| >= (M-f)/3` uses `W >= 8M log^+ C`.
* The heavy count is at most `32L/7`.
* `||mu||_1 <= 2L`.
* The certificate vanishes on `F`, has zero sum, and has column image `v - v'`.
* `mu^t k` is an integer when records agree.
* There is no hidden content issue: no nonconstant column is constant on `U`, since `f < M/2`, so the
  sub-cluster gcd is still `g`.
* The pigeonhole inequality holds with margin at every multiple of 4 in `[4096, 70000]`
  (`ref_constants.py`). It in fact already holds at `M = 256`.

**Misstated comparison.** "Item 277 allows only a bounded number of deleted rows" is inaccurate.
For `k = 0`, item 277's record argument (`N-1 > 2^d`, cutoff `U(0,d) = max(2d+1, 2^d+2, (256K)^2+1)`)
already allows `d ~ log_2 N` deletions, in the unit regimes. Theorem B'' still improves this
substantially: `M/(40 log M)` deletions, every `C`, elementary.

## 3. Theorem C (arbitrary flipped row set): correct; novelty modest

**Step (ii) is fine.**

* On an affine flat `Y` inside the unflipped set, each physical column restricts to an affine
  character.
* Merged blocks over distinct split primes are conjugate-primitive, and constant restrictions become content.
* So `Y` carries a pure Sylvester core of order `S_0`, excluded by Corollary B.2 (`C <= 1`) or B.3.
* I re-derived the density-lemma recursion `alpha_{i+1} >= alpha_i^2/2` and the threshold
  `2(K/M)^{1/K}`.

**Step (iii) is the claim that item 333 runs "verbatim" with unflipped rows as bad rows.**
I re-derived item 322's identity (1) from Theorem A with flip masses `f_x = 0` on unflipped rows:

```text
V_c = W0/2 - kappa/2 + Psi/2 + 2h Fcell/M - 2F/M + Fp - 2Fmatch
```

It holds exactly, including on cosets that contain unflipped rows. One point is not addressed in the
write-up but is true: the coset composite is always a nonunit, because every label has `b >= 5`
positive-weight copies and `ell != 0`. So (4) is valid on every coset and the averaging (6)-(8) goes
through.

The remaining ingredients also survive:

* (2a) uses only `F_a <= W_a`;
* the bad-row average gains `N|U|/M <= 2.23 N M^{-1/32}`, which tends to 0 even at `N <= 2^{log* M}`;
* the surviving aligned flat consists of flipped rows with flipped prime log `>= (log M)/2`.

`ref_partial_flip_identity.py` checks all of this exactly (Fractions) on 246 coset certificates, at
flip densities 0 to 1 and `t = 3, 4, 5`:

* the pair average `E_x G_{x,x+d} = qhat(d)`;
* identity (1) against literal physical columns;
* the nonunit composite;
* the identity behind (2a).

**Novelty.** Item 278 (k = 0), applied to the unflipped sub-cluster with density `delta`, gives
`|U| < 2J 2^{(s-1)/2^{s-1}} M^{1-2^{1-s}}`. That already makes the unflipped density tend to 0, so
growth for arbitrary `F` follows from items 278 + 333 by the same patch, without Theorem B. What is
genuinely new:

* the bad-row absorption itself, a short routine observation not recorded in the notes;
* the explicit exponent `1 - 1/32` at `C <= 1`.

The write-up does acknowledge that (ii) is the k = 0 case of item 278.

**Small gap for Theorem 8.2.** Definition 8.1 frames have between 5 and `B` primes per label, and
prime-power blocks. Theorem C and item 333 are stated with exactly `b` copies of distinct primes.
The proof only uses capacity `<= B`, `r >= 5(M-1)` and distinct primes, so this adapts, but it should
be said explicitly.

## 4. Proposition P and the "exact barrier" claim: overstated

* **Parts 1-3 and 5 are already in the notes.**
  * Pair: item 334.
  * All integral/half-integral characters for every capacity-two assignment, `b >= 5`, `M >= 36`:
    items 391-393.
  * Item 392 states explicitly that comparing flipped columns with their unflipped copies makes `c/2`
    cover the whole saturated rational lattice. That is Lemma F.1; see also item 321.
  * Part 3 ("every rational character") is therefore a restatement.
  * Item 4 (aggregate collisions) is a correct direct count.
* **The new content** is the contrast: pure Paley fails rational characters (Theorem B), and Paley with
  `<= M/(40 log M)` flipped rows fails them for `M >= 4096` and `M >= C^{O(1)}` (Theorem B'').
* **"The barrier is exactly 'non-Walsh Hadamard core with almost all / at least M/(40 log M) rows flipped'"
  is not established:**
  1. For Paley with `M/(40 log M) < f < M` flipped rows, neither direction is proved. Lemma F.1 needs
     every row flipped. B'' needs `(2L+1)^f M < binom(n,L)`. A Siegel-lemma certificate on
     `ker(H_F)` at `f ~ M/2` has `||v||_1 ~ M^{3/2}`, far above half-height.
  2. Proposition P is proved only for prime-order Paley (`q = 3 mod 4`) with capacity two. "Non-Walsh"
     cores such as `H_{2^s} ⊗ P_q` contain exact Walsh sub-cubes of order `2^s` in each fibre. Nothing
     shows such cores pass the aligned-flat tests, so the dividing line is not "Walsh vs non-Walsh".

  Correct form: within prime-order Paley capacity-two profiles, the fully flipped ones pass and the
  `<= M/(40 log M)`-flipped ones fail; the intermediate range and other Hadamard types are open.
* **Table 7.4 is inconsistent with this.** It marks pure Paley "integral characters: pass". That holds
  only at the Liouville level, since item 276 excludes pure Paley with an integral character and
  Matveev.

## 5. Smaller points

* **Section 10, "Equivalently, by Lemma 7.2".** The pinning system is necessary, not equivalent.
  * Pinning errors `|eta_a| <= Delta/2` reconstruct row phases only to `sum_a |eta_a|`, up to `(M-1)Delta/2`.
  * The `m_a` must also lie in the image `H^t k` of an integer vector.
* **Section 7.2, "per-pair collisions with `q <= 0.9(b-4) log M` are satisfied for any residues".**
  This is true only if each prime is counted at the first level, i.e. `2 theta(Q)`. With full
  valuations `v_q`, higher powers are not controlled by residues.
* **Section 9.2.** "Odd quartic support `<= 4`" is necessary but not sufficient for a flipped 2-frame.
  * 64 of 4096 four-odd-prime pattern tuples admit no distinct flip points (`ref_frame_criterion.py`).
  * Definition 8.1 needs `>= 5` primes on each of the 3 labels, i.e. `>= 15` varying primes. That is
    impossible for most census `k` (8 to 16).
  * So the census measures a weaker notion than EH's frames. The conclusion "data give no support to
    EH" is a heuristic remark at `M <= 16`, not evidence about the asymptotic hypothesis.
* **Theorem A** is a correct two-line identity. I re-derived (A.2) and the aligned-flat specialisation
  matching item 322(2). It is a reformulation, not a new constraint.
* **Lemma 7.3** (triangle-area sector uniqueness) and **Lemma 7.2** (flip-shifted pinning) are correct.
* **Author's checker part 3** (`check_hadamard_residue_pigeonhole.py`) is vacuous. Random cores have
  `log C_act ~ W/4`, so (B.1) holds trivially. The real content is in parts 1, 2, 4 and 5.

## 6. Corrected statement

1. **Theorem B.** Let `H` be any normalised Hadamard matrix of order `M >= 28`. Any pure Hadamard-core
   cluster on an arc of length `<= C sqrt R` satisfies (B.1). Hence `M <= 24` for `C <= 1`, and
   `M < max(28, C^{2+o(1)})` in general. The proof is elementary and effective.
   * New relative to item 276 (any Hadamard type, but Matveev and `C <= sqrt 2`/`2`): every `C`,
     explicit small constants.
   * New relative to item 278: non-Walsh types.
2. **Theorem B''.** For `M >= 4096`, `f <= M/(40 log M)` deleted or modified rows, and
   `W >= max(8M log^+ C, 32 log^+(LC))` (automatic for `M >= C^{6+o(1)}`), no such configuration exists.
   This extends item 277's `~log_2 M` deletions in unit regimes.
3. **Theorem C** holds as stated. The growth part also follows from items 278 + 333 via the same
   bad-row absorption. Theorem B contributes the explicit exponent `1 - 1/32` at `C <= 1`.
4. **Proposition P.** Fully flipped capacity-two prime-order Paley profiles pass pair, all
   integral/rational character (from items 334/391-393), aggregate-collision and capacity tests at
   `W = Theta(M log M)`. Pure Paley and Paley with `<= M/(40 log M)` flipped rows fail rational
   characters. Whether intermediate flip counts or other non-Walsh Hadamard types are barriers is open.
5. **Theorem A and Lemmas 7.2-7.3** are correct. Section 10's system is a necessary condition only.

## Files

`referee_walsh_checks/` contains:

* `ref_partial_flip_identity.py` (Theorem C step iii, exact);
* `ref_pinning_literal.py` (pinning, (1.2) and Lemma 1.1 on literal Gaussian cores; order-4 short-arc search);
* `ref_frame_criterion.py` (Section 9.2 criterion);
* `ref_constants.py` (Corollary B.2 and Theorem B'' integers).

All pass. The repository `/home/user/Jarnik` was not touched.
