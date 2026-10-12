# Round 8: the growth rate (refereed)

Each write-up of this round was attacked by an adversarial referee whose report sits next to it
(`referee_*.md`). **Every statement below is the referee-corrected version.** The question of the round is the
growth rate only: does any input, possibly a standard conjecture, give `M = o(log R/loglog R)` for actual circles,
even with a loglog gain? Constants are quoted only where they measure a barrier.

Conventions:
* `N = R^2`, `W = log N = 2 log R`. An arc of length `<= C sqrt R` has angular width `Delta <= C e^(-W/4)`.
* `w(v) = sum_j |v_j| log p_j`, and `supp v` is the set of primes with `v_j != 0`.
* *Fake circles* (angle models) have genuine primes, a genuine exponent profile and real angles
  (`round4/outside.md`, `round5/construct.md`).
* Labels: **proved** means proved with all quantifiers and refereed. **Draft** means unrefereed. **Evidence** means a
  finite computation that no proof uses.
* Paths of earlier rounds are relative to `docs/claude_review/`. The unrefereed round-7 draft `round7/p3G.md` exists
  only in the session scratchpad; `round7/p3lower.md` is in the repository.

## 0. Answer to the round's question

**No growth improvement for actual circles was proved this round, either unconditionally or from a standard
conjecture.** Unconditionally, `M <= (2/3 + eps) log R/loglog R` (`two_thirds/proof.md`) is still the best bound.

* **Lang-Waldschmidt** (`lw.md`, survives as a barrier). Under the angle-model inputs, it gives no growth through
  any single linear form (representations over `Q(i)`), in either the max form or the product form, for any
  `eps >= 0`. At a fixed `eps` it gives
  none through all forms at once, but this is proved only in the max form with `Log i` among the logarithms. Open:
  the full conjecture with `eps -> 0`, and the standard product form simultaneously.
* **Simultaneous Lang-Waldschmidt** (`weakest.md`). The claimed conditional growth theorem
  `SLW^log => M = O((log R)^(1-eta))` collapses: the hypothesis is **false** for genuine algebraic numbers
  (polynomial templates, `referee_weakest.md` Props R1-R2). With coefficient-weighted heights the implication is
  repaired, but the repaired hypothesis is new, open, false for `k <= 5` forms, and supported by no heuristic.
* **Unconditional structure** (`stability.md`, `weakest.md` Thm 1.3). Two reductions hold for actual circles:
  * the uniform theorem is equivalent to bounded rank of the exponent lattice;
  * near-extremal clusters live on `M^(O(1))`-smooth norms with `Theta(M)` prime factors.

  The fake circles share both properties (full rank, smooth norms), so neither reduction separates actual clusters
  from fakes.
* **The only unconditional growth statement is restricted.** `M = o(log R/loglog R)` along circles with
  `log P_+(N)/loglog R -> infinity`. This is a two-line consequence of Cilleruelo-Cordoba (1992), and every fake of
  size `O(M log M)` violates its input.

The single most promising open route (Section 3) is the lower-bound problem for the class P3-all at polynomial
moduli: genuine residues with the prime-level parity law, imposed on every monomial. It is the only local class in
the barrier map that is valid for actual circles, is not known to be capped, and needs no conjecture.

## 1. Referee-corrected statements

### 1.1 `lw.md`: Lang-Waldschmidt and the skew double core

**Referee verdict (`referee_lw.md`): survives as a barrier, with two major corrections of scope.**
* The write-up called the max form "the Lang-Waldschmidt conjecture". The standard statement is the product form,
  and the two forms are incomparable at equal `eps` when `n` is unbounded.
* The barrier holds relative to the angle-model input class only.

The proofs re-derived without error are Lemma 2.2, Theorem 5.1, Theorem 6.1, `flipped.md` Theorem L, and the
`exact.md` Thm 4.1(iii) margin `Phi(c) >= 2(M-1) + n`. The exact minimum `s = M/2 + 1` was re-certified by modular
rank.

*Hypotheses.*
* `LW_max(eps, C)`: for all `n` and all `alpha_k in Q(i)^*`,
  `|Lambda| >= C^n B^(-(n-1+eps)) prod_k A(alpha_k)^(-(1+eps))`, where `A(alpha) = max(H(alpha), e^(|Log alpha|))` and
  `H` is the naive height.
* *The Lang-Waldschmidt conjecture* is the product form `|Lambda| >= C(eps)^n B/(|b_1 ... b_n| A_1 ... A_n)^(1+eps)`.
  Sources: Lang 1978, pp. 212-217; Waldschmidt, Mosc. Math. J. 4 (2004), Conj. 2.5, stated for rational integers;
  Waldschmidt 2000, Ch. 1. The section numbers are unverified (network blocked).

**Theorem A (single forms; proved).**
1. *Valuation lemma.* For every representation of a character form `sum_j v_j Log u_j + b Log i` by numbers of
   `Q(i)^*`, `sum_k log H(alpha_k) >= w(supp v)`. Single numbers attain this.
2. Under either form of LW, a single linear form can exclude an `M`-point cluster on an arc `C_0 sqrt R` only through
   a character with `w(supp v) < W/4 + log(1/C_0)`.
3. It does exclude the cluster when some character has all of the following, for `M >= M_1(eps, C(eps), C_0)`:
   * `K = O(W/(log M)^3)` distinct values;
   * coefficients `e^(O((log M)^2))`;
   * `w(supp v) < W/(4(1+eps)) - (1+eps) W/log M`.

   The write-up's error bound `W/(2 log M)` was false; this is the corrected form.
4. `round4/flipped.md` Theorem L is correct. Under LW, the genuine versions of the `construct`, `lower` and `sharp`
   (P3-char) fake profiles do not exist at Paley orders `M = q + 1 >= M_0` (ineffective). At other `M` these fakes
   are row-deleted sub-profiles, which Theorem L does not cover.

**Theorem B (skew double core; proved).** Take the Paley skew double core of `round5/exact.md` §4, of order
`M = q + 1` with `q = 3 mod 4` prime.
* Every nonzero zero-sum `c` is nonzero on at least `M/2` of the `2q` types. This is a two-basis uncertainty
  principle for `K, K^T`, proved mod `q`.
* At `M = 8, 12` the exact minimum is `M/2 + 1`, certified by rank. The minimisers are exactly the `2q` type
  vectors.
* With primes in `[P, P e^(1/(8q))]`, every character has `w(supp v) >= W/4 + (1/2) log P - 1/16`.
* So for this profile no single LW inequality, in either form, for any `eps >= 0` and any representation by numbers
  of `Q(i)^*`, beats the arc's a priori bound.

**Theorem C (simultaneous fakes at fixed `eps`; proved for the max form with `i` in the list only).** Let
`0 < eps <= 1`, `0 < C_0 <= 1` and `M = q + 1 >= M_2(eps, C_0)`. There are angle models on the skew double core with:
* `M` points on an arc `C_0 sqrt R`;
* every gap of `round4/outside.md` Lemma 0.1;
* Haar-generic angles;
* `LW_max(eps, C)` for all `C <= 1`, for every form in the logarithms of `Gamma_J` whose list contains `i`.

Their size is `W <= (max(8, 4/eps) + o(1)) M log M + 8 M log(1/C_0)`. Hence these inputs cannot prove
`M <= (min(1/4, eps/2) - eta) log R/loglog R`. The referee shows that the same method lowers `4/eps` to
`2/eps + o(1)`, which gives `min(1/4, eps)`. Point-level residue data can be added only through the unrefereed
`exact.md` Thm 5.1; the constant then becomes `max(12, 4/eps)`. Parts 1-2 of
Theorem 6.1 (the class `A_0`: profile, all monomial gaps, Haar-generic angles, with `eps = 1` giving `8`) are
self-contained.

**Conditional barrier (Prop 7.2, Lemma 7.1).** Call a family *`mu`-spread* if `w(supp v) >= (1/4 + mu) W` for every
nonzero character. A `mu`-spread family that also satisfies the round-5 margins would carry fakes for every `eps` at
`W = O_mu(M log M)`. Turning that into an `eps`-uniform barrier needs one diagonal step (`C_M(eps) -> 0`), which is
true but not stated. For row-permuted Paley bases, Lemma 7.1 shows every character with `|supp c| <= (1/4 - mu) M`
is spread.

**Not proved.**
* Simultaneous fakes for the standard product form at fixed `eps`. The union bound diverges, e.g. at
  `c = (B, 1, ..., 1)`.
* Simultaneous fakes for homogeneous instances (no `Log i`).
* Representations by numbers of higher degree.
* Anything `eps`-uniform.
* Compatibility of the double core with the prime-level parity law, or with additive data.

**Minor corrections.**
* At `M = 8` the census misses positive-dimensional families with `s = M`. The `M = 12` counts are complete.
* The lead's claim "support `>= (1 - 1/t) W` for `t` independent Hadamard bases" is refuted by pair characters for
  `t >= 3`. For `t = 2` it is refuted only by coincidences at `M = 8, 12`.
* "Plotkin-optimal" and "best known fake" rest on the unrefereed `exact.md`.

The `1/eps` is intrinsic to union bounds over representations: there are Bell-number many groupings,
`e^((M/2) log M (1 - o(1)))`.

### 1.2 `weakest.md`: which Diophantine hypothesis suffices

**Referee verdict (`referee_weakest.md`): the headline does not survive (two fatal issues).** Every conditional
theorem of §2 rests on a false hypothesis. The unconditional results survive.

*The hypothesis.* `SLW(k, g, kappa)` concerns `k` forms with independent integer coefficient rows, in logarithms of
multiplicatively independent `alpha_j in Q(i)` with `|alpha_j| = 1`:

```text
max_i |Lambda_i|  >=  c^(n+1) prod_j (1 + B_j)^(-kappa) prod_(j in U) A(alpha_j)^(-g).
```

Here `U` is the union of the supports. Each number is charged once, whatever its coefficients. `SLW^log` asks this
only for systems with `k >= max(k_0, theta log(n+2))`.

**Corrected statement.**

*Conditional part: false or vacuous.*
1. `SLW(k, g, kappa)` is **false** for every `k >= 1`, `kappa >= 0`, `c in (0,1]` and `g < 1` (Prop R1).
   * *Counterexamples.* Polynomial templates `beta_i(t) = prod_l (t/b_l + i)^(sigma_il)` with
     `sum_l sigma_il b_l^(2j-1) = 0` for `j <= s`. They have private primes, `max|Lambda_i| ~ t^(-(2s+1))` and union
     height `~ t^(2D)`.
   * At `k = 5` the hypothesis fails with `{0, +-1}` coefficients for every `g < 1/4`, and at `k = 4` for every
     `g < 3/10`.
   * There are finite, exactly computed violations of `SLW(5, 1/4, 1)` with `c = 1/2`.
   * `SLW^log(theta, g, kappa)` is false for every `theta > 0` and `g < 1/2` (Prop R2).
   * Consequences:
     * Theorem 2.3 (`M <= 458`), Theorem 2.4 (`M = O((log R)^(1-eta))`) and Prop 2.7 are vacuous.
     * Prop 2.5 is true but not sharp.
     * The Dirichlet heuristic for `SLW` is invalid, because one-parameter families beat it.
2. *Repair with coefficient-weighted heights* `prod_j Mah(alpha_j)^(-g B_j)`. The proofs then go through:
   * `SLW_w(k, g, kappa)` for some `k >= 6` and `g < 1/4` implies the uniform theorem with `M <= (5.5^(k-1) + 1)/2`
     for `R >= R_0(C)`;
   * `SLW_w^log(theta, g, kappa)` with `k_0 >= 6`, `theta < 1/log 5.5` and `g < 1/4` implies `M = O((log R)^(1-eta))`.

   Status of `SLW_w`:
   * it holds at `g = 1/2`, which is Liouville;
   * it is false for `k <= 5` at `g < 1/4` (the same templates);
   * it is open for `k >= 6`. It is new, and no heuristic supports it.

   It concerns 7-point configurations, like `V(7, 1/4)`. The circle-restricted variant `SLW_circ` is the
   conclusion restated with a power margin.

*Unconditional part: correct.*

3. **Rank reduction (Thm 1.3).** Let a normalised configuration lie on an arc of length `<= C sqrt R`, let its
   exponent lattice have rank `k`, and let `log N' >= 18(2 log^+ C + log 2)`. Then

   ```text
   M <= ((1 + 2W/rho)^k + 1)/2 <= (5.5^k + 1)/2,      rho = W/2 - 2 log C - log 2.
   ```

   * The uniform theorem is equivalent to bounded rank.
   * Rank `<= theta log W'` with `theta < 1/log 5.5` would give `M = O((log R)^(1-eta))`.
   * Checked exactly on 73,728 windows of actual circles.
   * This is a reduction only. No rank bound is known, and the fakes have full rank `M - 1`.
4. **Support lemma (Lemma 3.2).** On the skew double core every nonzero rational character is nonzero on at least
   `(q+1)/2` cuts. The true minimum is `(q+3)/2` for `q = 7, 11`, not the `(q+7)/2` the write-up states. The columns
   `K e_a` attain `(q+3)/2` for every `q <= 67` tested.
5. **LW blindness (Thm 3.3, Lemma 3.4).**
   * Lang's product-form LW, with `eps <= 1`, any constant and any `b_0` (homogeneous forms included), holds
     simultaneously for every cluster-determined character form, written in the basis `alpha_T`, of a
     random-position angle model on the skew double core. The probability is `>= 1 - e^(-eps q tau/8)`.
   * On the product form this improves on `lw.md`, but only for this restricted family of forms.
   * The fakes' existence at `W = O(M log M)` with point-level residues uses the unrefereed `exact.md` Thm 5.1.
6. **Other routes.**
   * abc over `Q` or `Q(i)` is vacuous on squarefree profiles (pair, triangle, Ptolemy and `n`-term relations), and
     all fakes are squarefree.
   * Bombieri-Lang is vacuous, because the auxiliary varieties are rational.
   * Theorem 5.2: if `QS(m, delta, f, H_0)` holds for some `m >= 2` and `0 < delta < (m-1)/2`, then
     `M <= f(omega(N')) + m - 1`. So `f = o(s)` gives growth (definition in §2.3).
   * Prop 5.3: that hypothesis is false for `m <= 3`.
   * One instance `V(m_0, eps)` with `eps < (m_0 - 6)/2` gives the uniform theorem. The write-up's `(m-6)/4` is
     sufficient, not the threshold.
   * The integrality set `S` of the chord data contains every inert prime `<= M - 2`, so Lang-Vojta would need
     `S`-uniformity.
   * Prop 6.1(c), "no Vojta-type hypothesis gives growth without the uniform bound", is an unproved meta-claim.

### 1.3 `stability.md`: structure of near-extremal clusters

**Referee verdict (`referee_stability.md`): survives as a restricted theorem plus a barrier.** No mathematical error
was found in Sections 1-7. The checks were exhaustive for Theorem 3.1 on 435,196 windows (`N <= 5*10^6`, `C <= 8`),
plus about 245,000 further windows and subsets.

There are two major corrections:
* *Attribution.* The crude gcd lemma, the crude smoothness theorem and its growth corollaries follow in two lines
  from the Cilleruelo-Cordoba theorem (Proc. AMS 115 (1992) 899-905) after dividing by a shared Gaussian factor. The
  CC statement is quoted from memory; check it against the paper. What is new is the `F*` second-order gain and what
  is built on it.
* *Section 8.6.* Its pointer "five simultaneous forms suffice" rests on the false `SLW` (§1.2) and is withdrawn.

Notation: points `z_1, ..., z_M` with chords `<= C sqrt R`, normalised as in `two_thirds` (`G, N', W'`).

**Corrected statement (all unconditional).**
1. **gcd lemma (Prop 2.2).** For every subset `B` with `|B| = s >= 2`,

   ```text
   log Norm gcd_B <= W/s + (4(s-1)/s) log C_B - 8F*(s)/s^2 < W/s - 3 log s + K_0 + 4 log^+ C,      K_0 < 13.9092.
   ```

   Equality holds iff `s = 2` and `|beta|^2 = 2`. Without the `F*` gain the bound follows from Cilleruelo-Cordoba:
   * `W/s + (4(s-1)/s) log C - (2(s-1)/s) log 2` for even `s`;
   * `(W + 4s log C - 2s log 2)/(s+1)` for odd `s`, which is a better exponent.
2. **Smoothness (Thm 3.1, Cor 3.2).** For `M >= 4`, every prime `q | N` satisfies
   `log q < 2W/M - 3 log(M/2) + K_0 + 4 log^+ C`. The crude form `q < N^(2/M) C_+^4` (also from CC plus a halving
   step) gives:
   * `M < 2W/(log P_+(N) - 4 log^+ C)`;
   * `M = o(log R/loglog R)` whenever `log P_+(N)/loglog R -> infinity`, with `C` fixed;
   * `M < max(4, 2/lambda)` if `P_+(N) >= N^lambda C_+^4`.

   These are **restricted** theorems. Every model in `A_1^>(5M/4)` with `W = O(M log M)` has
   `log P_+ <= 2W/M + O(1) = O(log M)`; in `construct.md`, `log p = (3 + o(1)) log M`.
3. **Near-extremal clusters.** Assume `M >= eps log R/loglog R` and `R >= R_0(eps, C)`. Then:
   * `P_+(N) <= M^(4/eps - 3 + o(1))`, and `p^e || N'` implies `p^e < C_+^8 (log R)^(8/eps)`;
   * `(M/4)(1 - o(1)) < omega(N') <= Omega(N') <= (2/eps + o(1)) M`;
   * for `C <= 1`: `omega(N') > M/4`, more than `M/2` distinct cut partitions, and any two rows separated by more
     than `M/4` cuts;
   * the primes below `M^(1-delta)` carry an `O(M^(-delta)/log M)` fraction of `W'`.
4. **Profile and collision stability (Props 4.1, 5.1, Thm 5.2).** Put `A = W'/(M log M)`.
   * The mean pair margin is `Theta(log M)`:
     `(M/(M-1))(log M - 7.9252) + log 2 - 2 log^+ C < mbar <= W'/(2(M-1))`.
   * The geometric mean of `|beta|^2` lies in `[M^(1-o(1)), C^2 M^(1/eps + o(1))]`.
   * The chord, collision and level-profile excesses are each `<= (A - 3)/3 + O(1/log M)` times `F*(M)`.
   * `sum over split p <= M with p not dividing N' of log p/(p-1)` is `<= ((A-3) log M + K')/2`.
5. **Genuine input (Prop 6.1; Matveev 2000, which should be quoted precisely).** In the near-extremal regime, every
   zero-sum `c` with `||c||_1 <= exp((log R)^(eta/2))` and `v = sum c_x a_x != 0` has
   `|supp v| > (1 - eta) loglog R/logloglog R`.
6. **Barrier (Prop 7.1, Cor 7.2).**
   * Items 1-4 hold in the angle-model class `A_1^>(Q_0)` with `Q_0 >= 5M/4`, with `F_free` in place of `F*`. In
     that class Thm 5.2 reads `A >= 4 - o(1)`.
   * The refereed `construct.md` models lie in this class (`W <= (27 + o(1)) M log M`, renormalised after row
     deletion when `M != q + 1`), and their profiles satisfy item 5.
   * So `A_1(Q_0)` together with items 1-5 cannot prove growth.
   * The class does not contain the prime-level parity law, residues at all moduli, or additive relations.
7. **Section 8 routes.** S-unit counting, small-`S` relations (which are parallel chords), abc and smooth numbers near
   squares give no growth. GRH counting is heuristic and reaches only the Liouville scale.

The minor corrections are:
* the dropped `C <= 1` and `log^+ C` qualifiers;
* "integrality of reduced chords is equivalent to genuineness", where only one direction is shown;
* the profile-level-only fake check;
* the gaps in the transfer: `v_q(beta)` in models, and renormalisation after row deletion;
* the unproved "complete content of recursion" remarks;
* the `lw.md` labelling;
* the abc coprimality reduction.

## 2. Growth-rate map

### 2.1 Dictionary

* For a near-extremal cluster, `log M = (1 + o(1)) loglog R`. So `W = c M log M` corresponds to
  `M = (2/c + o(1)) log R/loglog R`.
* Actual clusters satisfy `W >= 3 M log M - O(M)`, hence `M <= (2/3 + eps) log R/loglog R` (`two_thirds/proof.md`).
* **A growth improvement means `W/(M log M) -> infinity` along near-extremal actual clusters.** A loglog gain,
  `W >> M log^2 M/loglog M`, would give `M = O(log R logloglog R/(loglog R)^2)`.
* *Capped.* An input class is capped when it contains fakes with `W <= c M log M`. Then no argument using only those
  inputs proves `M <= (2/c - eta) log R/loglog R`, and a fortiori none proves `o(log R/loglog R)`. Every cap below is
  at the rate `log R/loglog R`; only the constant `c` varies.

### 2.2 Inputs proved capped at `log R/loglog R`

The classes are nested from top to bottom, except that LW rows are additions to `A_0`. Each row lists the inputs, a
witness family of fakes that satisfies them, its size, and its status.

| # | inputs | witness | `W` of the witness | status |
| --- | --- | --- | --- | --- |
| 1 | `A_0`: exponent profile, pair chords, cut/Plotkin, all monomial (Liouville) gaps of `outside.md` Lemma 0.1, rational characters (equivalent to the monomial gaps for large primes, `outside.md` Thm A.1), GAP/Fourier, energy and incidence arguments as reduced in `outside.md`, Haar-generic angles | skew double core (`lw.md` Thm 6.1 parts 1-2) | `(8 + o(1)) M log M + 8M log(1/C)` | proved (round 8; it uses `exact.md` Thm 4.1(iii), which `referee_lw.md` re-derived) |
| 2 | `A_1(Q_0)`: `A_0` plus point-level residues at moduli `<= Q_0` (small-prime collisions, gcd and cofactor data), with each character's actual collision modulus | `construct.md` fully flipped Paley cores | `(27 + o(1)) M log M + 9M log(1/C)` at `Q_0 = 2M`; constant `9(2 + max(K,1))` at `Q_0 = M^K` | proved (round 5). `lower.md`: `(16 + 32A)` at `M^A`. `exact.md` double core: `12` (draft) |
| 3 | P3-char: `A_1` with residues at moduli `<= 3M^A` that factor through prime residues with the genuine norm, so the prime-level parity law holds, imposed on characters | `sharp.md` parity-class designs | `(b max(A,2) + o(1)) M log M`, `b >= 33` | proved (round 6) |
| 4 | `A_1^>(5M/4)` plus everything in `stability.md` §2-5: sub-cluster gcds, smoothness `P_+ <= N^(2/M) C_+^4`, prime counts, Plotkin balance, margin bounds, collision, cofactor and level stability | `construct.md` models (transfer, Prop 7.1) | `(27 + o(1)) M log M` | proved (round 8); the renormalisation after row deletion still has to be written out |
| 5 | proved linear forms in logarithms (Baker, Matveev): they exclude characters supported on `<= (1 - eta) log M/loglog M` primes (`stability.md` Prop 6.1). Baker-type bounds (product of heights) are too weak for Theorem L (`flipped.md`) | `construct.md` (supports `>= 2 sqrt((b-1)M) - 1`), double core (supports `>= M/2`) | as in rows 1-2 | proved |
| 6 | `A_0` plus **any single LW inequality**, product or max form, any `eps >= 0`, any representation by numbers of `Q(i)^*` | skew double core (`lw.md` Thm B with Thm 6.1 parts 1-2; `weakest.md` Thm 3.3) | `(8 + o(1)) M log M`; `(12 + o(1)) M log M` with residues (draft `exact.md` Thm 5.1) | proved relative to `A_0` |
| 7 | `A_0` plus `LW_max(eps, C)` at a **fixed** `eps`, for all forms whose list contains `i`, simultaneously | skew double core (`lw.md` Thm 6.1) | `(max(8, 4/eps) + o(1)) M log M`; `2/eps + o(1)` by the referee's refinement | proved; the constant degrades as `eps -> 0` |
| 7' | `A_1(2M)` plus product-form LW (any constant, `eps <= 1`) for all cluster-determined character forms in the basis `alpha_T`, simultaneously | skew double core in random position (`weakest.md` Lemma 3.4) | `O(M log M)` | LW part proved; the existence of the fakes (gaps, residues) rests on `exact.md` Thm 5.1 (draft) |

**Certificate-level cap (draft).** `round7/p3lower.md` Thm A treats every inequality derived from the pair
inequalities together with the cut/Plotkin/PSD structure of margins, with real weights. Such an inequality cannot
exceed `pi lambda_max(D)`, where `D` is the pair-demand matrix. For the P3-all members of `sharp.md` Thm 6.4,
`lambda_max(D) <= 36 M log M + 111 M`. So no pair-cut certificate proves growth in P3-all, even though those members
have `W` of order `M log^2 M/loglog M`. This is unrefereed.

**Scope of every cap.** The witnesses are angle models. None of them is known to satisfy:
* the prime-level parity law imposed on every monomial (P3-all);
* residues at moduli growing faster than a power of `M`;
* additive (archimedean) cluster data such as `Log((z_x + z_y)/(2 z_u)) = O(Delta)`.

The construct and lower fakes violate the parity law (`referee_construct.md`). The skew double core has never
been given P3 data.

### 2.3 Hypotheses that give growth, and restricted growth theorems

**(a) Conditional theorems that give growth or better.** In every row the implication is proved and refereed, except
the last row, which is a draft. Each hypothesis carries information that the angle-model inputs of §2.2 do not
imply (otherwise the fakes would contradict the conclusion), so no row conflicts with §2.2.

| hypothesis (precise form in the source) | conclusion | status of the hypothesis | source |
| --- | --- | --- | --- |
| Vojta's Main Conjecture (Vojta, LNM 1239, 1987) with `D = 0`, `A = H` and one fixed `eps < 1/2`, for a smooth model of the blow-up of `X_7 = {x_1^2+y_1^2 = ... = x_7^2+y_7^2}` in `P^13` along its diagonal line | uniform: `M(C) <= max(M_0, c_0 C^22 + 1)` | standard conjecture, open | `growth/conditional.md` |
| the same for one instance `V(m_0, eps)`, any fixed `m_0 >= 7`, `eps < (m_0 - 6)/2` | uniform | standard, open | `weakest.md` Prop 6.1(a,b), threshold corrected by the referee |
| Vojta (`D = 0`) on `Mbar_{0,8}`, or the arithmetic lemma (L_8) | uniform | standard; (L_8) is open and Vojta-type | `growth/ptolemy.md` |
| `SLW_w(k, g, kappa)`: as `SLW` in §1.2, with `prod_(j in U) A(alpha_j)^(-g)` replaced by `prod_j Mah(alpha_j)^(-g B_j)`; `k >= 6`, `g < 1/4` | uniform: `M <= (5.5^(k-1) + 1)/2` | new; false for `k <= 5`; open for `k >= 6`; no heuristic support | `referee_weakest.md` §3.1 |
| `SLW_w^log(theta, g, kappa)`: `SLW_w` for all systems with `k >= max(k_0, theta log(n+2))` forms; `k_0 >= 6`, `theta < 1/log 5.5`, `g < 1/4` | `M = O((log R)^(1-eta))`, `eta = 1 - theta log 5.5` | new, open | same |
| `QS(m, delta, f, H_0)` for some `m >= 4` and `1 < delta < (m-1)/2` (definition below) | `f` bounded: uniform; `f(s) = o(s)`: `M = o(log R/loglog R)`; `f(s) <= s^(1-eta)`: `M <= (log R)^(1-eta+o(1))` | non-standard. False for `m <= 3`. The known bounds (Evertse 1996; quoted from memory) are `exp(O(s))`, which is useless. For `m >= 4`, `delta in (1, (m-1)/2)` it is neither refuted nor supported | `weakest.md` Thm 5.2, Prop 5.3, Remark 5.4 |
| rank of the exponent lattice `<= theta log W'` with `theta < 1/log 5.5` (resp. bounded rank) | `M = O((log R)^(1-eta))` (resp. uniform; the converse also holds) | open; the reduction itself is unconditional. The fakes have full rank `M - 1`, so the hypothesis needs non-local input | `weakest.md` Thm 1.3, Cor 1.4 |
| *(draft)* `K`-uniform pair margins and the Poisson-type factorial-moment Hypothesis `H(eta, delta)` on collision counts | `M <= (K/(4(1-eta)) + o(1)) log R logloglog R/(loglog R)^2` | neither hypothesis is known for actual circles; the round-6 design violates `H` | `round7/p3lower.md` Thm C (unrefereed) |

*Definition of `QS`* (`weakest.md` §5.1).
* Work over `K = Q(i)` with normalised absolute values: `|x|_inf = |x|^2`, and `|x|_w = N(w)^(-v_w(x))`.
* For a finite set `P` of split primes, with `s = |P|`, let `S_P` be `inf` together with the places above `P`.
* The system `Sigma_m(P)`: at `inf` the forms `X_1` and `X_i - X_1` (`i = 2..m`); at each finite `w in S_P` the forms
  `X_1, ..., X_m`.
* `QS(m, delta, f, H_0)` asserts: for every `P`, the solutions `x in K^m` of
  `prod_(v in S_P) prod_i |L_(i,v)(x)|_v/||x||_v <= H(x)^(-m-delta)` with `H(x) >= H_0` lie in at most `f(s)` proper
  linear subspaces.

The classical Subspace Theorem gives finitely many subspaces for each fixed `S`. What is needed here is sublinear
dependence on `s = |S|`. `SLW_circ` (`referee_weakest.md` §3.2) is omitted from the table because it restates the
conclusion.

**(b) Restricted unconditional growth theorems.** Each input is violated by the fakes, so none of these extends to
general circles.

| restricted input | conclusion | why the fakes escape | source |
| --- | --- | --- | --- |
| `log P_+(N)/loglog R -> infinity`, `C` fixed | `M < 2W/(log P_+(N) - 4 log^+ C) = o(log R/loglog R)` | `log P_+ = O(log M)` | `stability.md` Cor 3.2; Cilleruelo-Cordoba 1992 |
| `P_+(N) >= N^lambda C_+^4` | `M < max(4, 2/lambda)` | same | same |
| pure Hadamard-core profile, any type | uniform: `M <= 24` for `C <= 1`, `M < max(28, C^(2+o(1)))` | fakes are flipped cores | `growth/walsh-extraction.md` Thm B |
| repeated-Walsh profiles, flips on any row subset | `M = o(log R/loglog R)` | the fakes are Paley, not Walsh | same, Thm C |
| one split prime `p` with bounded `v_p(L_0)` (resp. `L_0` free of split primes below `M/8`) | `M = O(sqrt(log R))` (resp. `M <= (sqrt 24 + o(1)) sqrt(log R/loglog R)`) | one colliding pair raises `v_p(L_0)` at negligible cost | `growth/cotangent.md` Thm A |
| abc over `Q`, and `rad(N) <= N^(1/4 - eta)` | at most one point per arc `C sqrt R`, for `R >= R_0(C, eta)` | fakes are squarefree | `growth/conditional.md` Cor 8.2 |

### 2.4 Hypotheses and methods that give no growth

| hypothesis or method | why it gives no growth | source |
| --- | --- | --- |
| `SLW(k, g, kappa)` with union heights, and `SLW^log` | false for genuine numbers: for every `k` at every `g < 1`, and for `SLW^log` at every `g < 1/2`. The implications are vacuous | `referee_weakest.md` Props R1-R2 |
| Lang-Waldschmidt, single forms (product or max form, any `eps >= 0`, `Q(i)`-representations) | sees a character only through `w(supp v)`; the skew double core has `w(supp v) >= W/4 + (1/2) log P - 1/16` for every character | `lw.md` Thms A-B, `weakest.md` Thm 3.3 |
| `LW_max(eps)` at a fixed `eps`, all forms with `i` listed | fakes at `(max(8, 2/eps) + o(1)) M log M` satisfy it | `lw.md` Thm C |
| Baker/Matveev (proved bounds) | reach supports `<= (1 - eta) log M/loglog M` only; the fakes' supports are `>= c sqrt M` or `>= M/2` | `stability.md` Prop 6.1; `flipped.md` |
| abc over `Q` or `Q(i)` | vacuous on squarefree profiles, and every fake is squarefree. Two-level profiles force `rad(N) >= N^(1/2 - o(1))` | `weakest.md` Props 4.1-4.2 |
| Bombieri-Lang | vacuous: every auxiliary variety (including the chord-data variety) is rational | `weakest.md` Prop 6.2 |
| Lang-Vojta for the chord data | the integrality set contains every inert prime `<= M - 2`, so `S`-uniformity would be needed, which is non-standard | `weakest.md` Props 6.3-6.4 |
| Subspace Theorem with `S = {inf} u primes(N)` and known `S`-dependence | `exp(O(s))` subspaces give only `M <= exp(O(s))`, weaker than `M <= 4 d(N)`. `QS(m <= 3, delta, o(s))` is false | `weakest.md` §5 |
| `S`-unit equation counts (Beukers-Schlickewei `2^(8(r+1))`; Erdos-Stewart-Tijdeman) | bounds of size `2^(Theta(M))` against `M^3` triangle identities | `stability.md` Prop 8.2 |
| single-place auxiliary polynomials (Liouville at `inf`), Ru-Vojta with `S = {inf}`, forced Gaussian divisibility | `tau*(M) = 2 - 1/ceil(M/2) < 2` and `beta_M < 2`; divisibility is exactly critical. Nothing at the `sqrt R` scale | `threshold/` |
| smooth numbers near squares | the cluster is the original problem rotated by `conj(u_1)`; no reduction | `stability.md` §8.3 |
| counting under GRH for Hecke `L`-functions of `Q(i)` (heuristic) | square-root cancellation reaches the Liouville scale even for one form | `stability.md` §8.4 |
| unconditional structure: rank, smoothness, Plotkin balance, collision stability | shared by the fakes (full rank, smooth norms, transfer Prop 7.1) | `weakest.md` Thm 1.3; `stability.md` §7 |

### 2.5 Undecided

1. **Full Lang-Waldschmidt.** The open cases are:
   * `eps -> 0` together with `M`;
   * simultaneous fakes for the standard product form at fixed `eps`;
   * homogeneous instances;
   * representations by numbers of higher degree.

   Two routes would give an `eps`-uniform barrier:
   * (S_mu): a `mu`-spread family satisfying the round-5 margins (Prop 7.2), e.g. `t >= 3` random row permutations of
     a Paley core. Lemma 7.1 handles `|supp c| <= (1/4 - mu) M`; the rest is a two-basis spark problem.
   * (N): a character-level union bound for the double core, weighted by coefficient size (`referee_lw.md` §4.2).

   *Heuristic:* the full conjecture gives no growth either.
2. **New simultaneous hypotheses.** `SLW_w(k >= 6, g < 1/4)` and `QS(m >= 4, delta in (1, (m-1)/2), o(s))`. Both
   would give growth or more; both are non-standard and unsupported. A single rank-6 circle family on arcs `C sqrt R`
   would refute `SLW_w(6)` at the boundary.
3. **P3-all at polynomial moduli** (`sharp.md` Def 1.1 with (e-all); question (P3-all\*) of `sharp.md` §6.7). Is
   `W_min^all = O(M log M)`?
   * The proved bracket at `X = 3M` is `3 M log M - O(M) <= W_min^all <= (66 + o(1)) M log^2 M/loglog M`.
   * Uniform random designs cannot beat the upper end (Prop 6.5, with the per-unit height condition).
   * See §3.
4. **The skew double core against residue data.** Its point-level residue version rests on the unrefereed
   `exact.md` Thm 5.1. Nobody has given it P3-char or P3-all data.
5. **Additive (archimedean) data.** These are not modelled by angle models. Their natural uses reduce to
   single-place auxiliary polynomials, which are barred, but no proved cap covers them.
6. **Prime-weight concentration.** Every fake puts its primes in one window `[P, P e^(1/(4r))]`. Whether
   near-extremal actual clusters must do so is open (`stability.md` §9.3). If they had to, the problem would become
   one about angles of Gaussian primes in very short intervals.

## 3. The most promising open route to an unconditional growth improvement

**Route: lower bounds for P3-all at polynomial moduli.**

*The class.* P3-all_X is `sharp.md` Def 1.1: conditions (a)-(d) with (e-char) and (e-all). For every odd `q <= X`,
each prime `p_j` gets a residue `r_j(q)` mod `q^(A_q)` with `Norm r_j(q) = p_j`. Every monomial `v != 0` must then
satisfy

```text
|sin(<v, phi> - s pi/4)|  >=  nu_s m_(v,s) e^(-w(v)/2),
```

where `m_(v,s)` is the collision modulus of `prod_j (r_j/conj r_j)^(v_j)` against `i^s`.

**Target (T).** For some fixed `A >= 1` and `X = 3M^A`, every member of P3-all_X(M, C) has `W >= f(M) M log M` with
`f(M) -> infinity`.

**What (T) gives.** Normalised actual clusters are members of P3-all_X for every `X`, with their genuine residues
(`sharp.md` Prop 1.2). The parity law that these residues obey is `sharp.md` Prop 2.1; `referee_construct.md` proved
it and checked it on 211,136 actual instances. So (T) gives, unconditionally,

```text
M <= (2/f(M) + o(1)) log R/loglog R = o(log R/loglog R).
```

If `f` matches the best known members (`f ~ c log M/loglog M`), this is `M = O(log R logloglog R/(loglog R)^2)`: a
loglog gain, which is exactly what this round asks for.

**Why this route rather than another.**
1. *It needs no conjecture.* Every input is a theorem about actual circles.
2. *It is the only local class not known to be capped.* `sharp.md` §6.7 already singles out (P3-all\*) as the one
   local input class whose growth rate is undecided. No known `O(M log M)` fake satisfies it:
   * `construct` and `lower` violate the parity law;
   * `sharp.md`'s parity-class designs reach P3-char at `O(M log M)`, but P3-all only at `O(M log^2 M/loglog M)`
     (Thm 6.4);
   * uniform random designs cannot do better (`sharp.md` Prop 6.5, as repaired by its referees);
   * the skew double core has never been given P3 data.
3. *The alternatives are worse.*
   * Every conditional route rests either on a standard conjecture whose effect is unknown or at best uniform (Vojta;
     full LW with `eps -> 0`, which heuristically gives nothing), or on a new hypothesis that nothing supports
     (`SLW_w`, `QS`).
   * Every unconditional structural fact found so far (rank, smoothness, Plotkin, collision stability, Matveev) is
     shared by the fakes.
4. *Round 8 narrows it.*
   * By `stability.md` Cor 3.6, a near-extremal actual cluster has every prime `<= M^(4/eps - 3 + o(1))`. So (T) is
     only needed for profiles on primes of size `M^(O(1))`, with moduli `X = M^A` of the same order.
   * By `construct.md` and `stability.md` Prop 7.1, the point-level class `A_1` is capped even with all of round 8's
     structure added. A proof of (T) must use the parity information on non-pair monomials.

**Constraints a proof must respect** (from the unrefereed round-7 drafts):
* *No pair-cut certificate proves (T).* For the `sharp.md` Thm 6.4 members, `lambda_max(D) <= 36 M log M + 111 M`
  (`round7/p3lower.md` A.5). A proof must therefore use one of two things:
  * the conditions (e-all) on monomials outside the pair lattice;
  * the discreteness of prime weights, i.e. few columns, each of the weight of one prime (`p3lower.md` Question Q).
* *Genuine residues are not a shortcut.*
  * Impose genuine residues on every monomial (class P3-G) over prime sets of polynomial size.
  * Then the certified closure pins every fake angle to the genuine one, so P3-G there coincides with the actual
    problem (`round7/p3G.md` Thms 3.1 and 4.1: computer-verified on specific prime sets, proved in a random-tightness
    model).
  * So (T) must be proved for P3-all, whose residues are constrained only by norms and parity.

**Decisive first test: the skew double core.**
* Its margins are exactly uniform. Every pair of rows is separated by exactly `M` of the `2q` types (`weakest.md`
  Lemma 3.1; exact check for all `q <= 83` in `readme_checks/dcore_pair_separation.py`). With all `2q` primes in one
  window `[P', P' e^(1/(8q))]`, every pair margin `mu_xy = d_xy - W/2` lies in `[log P' - 1/8, log P' + 1/4]`.
* The pair case of (e-char) with `eta = 1`, together with `|sin x| <= |x|`, gives `mu_xy >= D_xy - c_C`. Here
  `D_xy = 2 log m_xy` is the unit-1 pair demand and `c_C = 2 log C - log 2`. The mean margin is `<= W/(2(M-1))`
  (Plotkin). Hence on this profile, for every choice of residue data,

  ```text
  W >= (2(M-1)/K)(D_max - c_C),      K = max mu / mean mu = 1 + o(1)
  ```

  (`p3lower.md` Prop B.1; the proof is the two lines above).
* So the double core can be an `O(M log M)` fake in P3-char, respectively P3-all, **only if** it carries residue data
  of that class with maximal pair demand `D_max = O(log M)`. This condition is necessary, not sufficient: (e-char)
  for the other characters, and (e-all) for every monomial, must still hold.
* Random residues give, heuristically (Poisson), `D_max ~ 4 log^2 M/loglog M`. With that, `p3lower.md` Thm C (draft)
  gives `W >~ 8 M log^2 M/loglog M`. Designed residues can do better: `sharp.md`'s parity-class designs keep pair
  demands `O(log M)` in P3-char on Paley profiles with `b >= 33` copies per label. The double core has only two
  columns per label.

For each `M` the first step is a finite, exact computation: minimise `D_max` over the residue data of the class on
the double core. That gives evidence; the asymptotic statement still needs a proof.
* *If `D_max >> log M` is forced*, the double core is not an `O(M log M)` fake in that class. It then becomes the
  model case for (T): uniform margins, where pair demands alone already force growth. The next step is non-uniform
  margins. By `p3lower.md` Thm A, pair data alone cannot handle that regime, so there the (e-all) conditions or the
  discreteness of the weights must be used.
* *If low-demand data exist*, the remaining conditions have to be checked. For the Paley profile of `sharp.md` §4,
  Prop 6.6 there gives a sufficient condition: pair demands `O(log M)`, plus `m_v <= M^(beta ||v||_1)` for every `v`
  outside the character lattice, followed by the completion sum. The double-core analogue is not written. If the
  remaining conditions also hold, P3-all joins the capped list and every local class is capped. An unconditional
  growth proof would then have to use the actual angles of Gaussian primes of norm `M^(O(1))` through simultaneous
  approximation (`stability.md` §8.6). That is beyond current technology: counting under GRH reaches only the
  Liouville scale.

## 4. Files and caveats

| write-up | referee report | verdict | scripts |
| --- | --- | --- | --- |
| `lw.md` | `referee_lw.md` | survives as a barrier; scope restricted to the max form with `i` listed, and to the angle-model class | `lw_checks/`, `referee_lw_checks/` |
| `weakest.md` | `referee_weakest.md` | headline fails: `SLW` and `SLW^log` are false. Rank reduction, support lemma, LW blindness and the abc/Subspace/Vojta/Bombieri-Lang propositions survive | `weakest_checks/`, `referee_weakest_checks/` |
| `stability.md` | `referee_stability.md` | survives as a restricted theorem plus a barrier. The crude forms are attributed to Cilleruelo-Cordoba, and §8.6's pointer to `SLW` is withdrawn | `stability_checks/`, `referee_stability_checks/` |

Every script runs from its own directory with `python3 <script>`, and its output is stored next to it. The README
adds one check, `readme_checks/dcore_pair_separation.py` (exact; it verifies that every row pair of the skew double
core is separated by exactly `M` types for `q <= 83`).

**Caveats.**
* *Literature quoted from memory.* The network was blocked, so the following were quoted from memory:
  * Lang 1978 (pp. 212-217) and Waldschmidt 2000/2004 (Conj. 2.5); the product form was confirmed by two
    secondary sources;
  * Cilleruelo-Cordoba 1992 (statement also seen in a search snippet);
  * Matveev 2000 (only `c(n) <= c_0^n` is used);
  * Evertse 1996, Beukers-Schlickewei 1996, Erdos-Stewart-Tijdeman 1988, Rosser-Schoenfeld.

  None of these affects a proof except Cilleruelo-Cordoba, which is used only for attribution.
* *Unrefereed inputs.* The following results above rest on unrefereed drafts:
  * the point-level residue versions of the double-core fakes (`exact.md` Thm 5.1);
  * the "Plotkin-optimal / best known fake" description of the double core;
  * everything attributed to round 7 (`p3lower.md`, `p3G.md`).
* *Scope of every barrier.* Each barrier is relative to the angle-model input class it names. None covers the
  prime-level parity law on every monomial, residues at super-polynomial moduli, or additive data.
