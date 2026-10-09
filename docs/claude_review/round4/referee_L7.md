# Referee report on `round4/L7.md`: (L_7), residues, product criterion, structure of balanced seven-curves

Referee checks are in `scratchpad/round4/referee_L7_checks/` (Section 9). Conventions are the same as in
L7.md: *exact* means Python integers and `Fraction`; *NUMERICAL* means floating point and is never used
as a proof.

---------------------------------------------------------------------------------------------

## 0. Verdict

**The result survives in a corrected form, with no fatal errors.** Every rigorous statement in L7.md
is true, after the small repairs listed below. Three framing claims are wrong or overstated, and the
headline item (1) is not new.

| item | status after refereeing |
|---|---|
| Theorem 1 (`delta_p = 0`, `det(P_i,P_j) = ±prod N_T`) | **correct, not new**: it is `ptolemy.md` Lemma 5.2 together with that file's §5.2 remark (lines 492–494), which already states `delta_p = 0` and gives the same median reason. A three-line direct proof is in §1.2 below. Independent exact check: 0 mismatches at 9548 primes (families A and B of §1.3). |
| Corollary 1 (`(L_7) ⇒ (L_7)*`) | correct (one direction only); minor repair needed for the top block |
| Theorem 2(a) | correct |
| Theorem 2(b) | correct, but **a special case of 2(a)**. Remark 2 ("genuinely weaker than (a); no individual curve of class `∝ beta` is needed") is **false**. A product chart always produces Q-rational, approximately balanced curves as in (a). The proof also wrongly calls `mu` a flat morphism (repairable). |
| Theorem 2(c) | a sketch, correctly labelled as one |
| Proposition 3 | correct; all numbers re-derived exactly. Its conclusion ("the missing ingredient is a product structure") is **misdirected**; see §3 and §6. |
| Theorem 4 (i)–(iii), Lemma 5 | correct for rational curves of class exactly `beta` (essentially item 112 with `PGL_2` in place of `Aut(E)`); (iv) needs the extra hypothesis "the 56 contacts are distinct" in its statement |
| "the L6 mechanism cannot refute `(L_7)`" | **heuristic** (parameter counts), and only for class exactly `beta`. The summary states it as a consequence. |
| k = 6 numerics, 0/844 at k = 7, real walks, bookkeeping 7 / 8 / −24 / 31 | reproduced |

**Additional finding (NUMERICAL plus standard deformation theory).** L7.md left open whether the
"twisted cubics meeting 13 planes" of Prop. 3 exist and cover. They do, numerically:
* I found a nondegenerate curve of that class (`cub13`);
* the family of such curves is locally 3-dimensional and dominates `Mbar_{0,7}`;
* at least 8 distinct nondegenerate members pass through one general point.

Since `26·5040·beta = 70·sum_sigma sigma_*cub13 + 21·sum_sigma sigma_*f` exactly, smoothing trees of
free rational curves then gives **free rational curves of class `131040·beta` over `C` that cover
`Mbar_{0,7}`**. This rests on the numerical existence of the covering `cub13` family. So the obstruction to
refuting `(L_7)` through Theorem 2(a) is purely arithmetic: one needs Zariski density of the
**Q-rational** members. Neither curve classes nor a product structure is the missing piece. This
agrees with L7.md §6.3 ("This is the precise gap"), but contradicts the summary's "the missing
piece is building a product structure".

---------------------------------------------------------------------------------------------

## 1. Theorem 1 and Corollary 1

### 1.1 Correctness

I re-derived the proof. With `w_0, w_1, w_2 = (1,0), (0,1), (1,1)`, these three ends reduce to three
distinct points of `P^1(F_p)` for every `p`. So the standard vertex is their median. It lies on the
geodesic `(1,2)`, which is inside the subtree spanned by the ends `1..6`, and this gives
`L_{r_p} = b_p`. The proof in L7.md is correct.

### 1.2 A shorter, self-contained proof (referee)

Write `a_j = x_j/y_j` (`j >= 1`). Then `det(w_i,w_j) = y_i y_j (a_i - a_j)` and
`P_j = sgn(y_j)·(D a_j, 1)`. Hence

```text
v_p det(P_i,P_j) = v_p(D) + v_p(a_i - a_j),     v_p(D) = max_j v_p(y_j) = -min(0, min_j v_p(a_j)).
```

In the cluster picture with the anchor at `∞`, the sum of `l_T(p)` over the visible `T` containing
`i` and `j` telescopes along the chain of discs between `disc(a_i,a_j)` and the smallest disc
containing all `a_k`:

```text
sum_{T ∋ i,j, |T| <= 5} l_T(p) = v_p(a_i - a_j) - min_{k<l} v_p(a_k - a_l).
```

So `delta_p = min_{k<l} v_p(a_k - a_l) - min(0, min_j v_p(a_j))`. With `a_1 = 0` and `a_2 = 1`:

* `min_j v_p(a_j) <= v_p(a_2) = 0`;
* `min_{k<l} v_p(a_k - a_l) = min_j v_p(a_j)`, by `a_1 = 0` and the ultrametric inequality.

Hence `delta_p = 0`. The formula also shows that the normalization is needed: without it,
`delta_p > 0` occurs.

### 1.3 Exact checks (`thm1_cluster.py`, `adm_recount.py`)

`thm1_cluster.py` is independent of the Gromov-product min-formula used in `anchor_check.py`.

| family | prime checks | mismatches | `delta_p != 0` |
|---|---|---|---|
| (A) 300 random normalized configurations, coordinates up to 200 | 5183 | 0 | 0 |
| (B) 300 configurations with deep nested `p`-adic clusters, `p ∈ {2,3,5,7}`, normalized (124/300 trees non-admissible) | 4365 | 0 | 0 |
| (C) the same kind of data **without** the median normalization | 4620 | 0 non-constant offsets | **101 primes with `delta_p > 0`** |

Other checks:

* `anchor_check.py` reproduces as stated: 5208 primes, 0 mismatches, 568 non-admissible trees.
* `adm_recount.py` recounts with an exact brute-force two-chain test, not the greedy `chains_ok`, and
  also gets 568.
* Defect (minor): `anchor_check.py`'s tally "primes with top-block length > 0: 0" is **vacuous**.
  For `T = [6]` the complement `{0}` contains no pair, so its `edges()` returns 0 by construction.
  The real test of `delta_p = 0` is the mismatch count, which is valid.

### 1.4 Novelty: Theorem 1 is not new

`growth/ptolemy.md` already contains this.

* **Lemma 5.2** gives `v_p det(P_i,P_j) = delta_p + sum_{T ∋ i,j} l_T(p)`.
* **§5.2 (lines 492–494)** says: "No extra top block appears (`delta_p = 0`). The reason is that
  the standard vertex is the median of the three anchors `infinity, 0, 1`, and that median lies on
  the `infinity`-ray on the near side of `b_p`."
* L8 Lemma B.5 uses the same normalization `(1,0), (0,1), (1,1)`.

The argument is general, even though ptolemy states it while discussing item 109. L7.md Remark 1
presents the vanishing of the top block as a sharpening of ptolemy §5.2; in fact it is a
restatement. What is new is only the packaging, namely Corollary 1, the residue-free
reformulation.

### 1.5 Corollary 1

The proof is correct:

* `l_T(p) <= v_p det(P_i,P_j)` for `i, j ∈ T`;
* `H_Sigma = delta(Phi)` once the non-admissible `p`-parts are moved into the `t_ij`;
* the bookkeeping of `eps`, `eta_0` and `c` checks.

Minor issues:

* (i) The proof sets `n_[6] := 1`. Under the "all blocks balanced" reading of `(L_7)`, which §1
  mentions, the top block must be inflated by `diag(n,1)` with a fresh prime `n ≈ e^w`, and the
  singletons must be fresh primes. This is harmless, since `Y_j = ±1`.
* (ii) Only `(L_7) ⇒ (L_7)*` is proved. The summary's "So `(L_7)` is purely a statement about how
  the finite boundary heights … are distributed" claims more than that one-way implication.

---------------------------------------------------------------------------------------------

## 2. Theorem 2

### 2.1 (a): correct

I checked each step:

* the binary forms `F_S` are coprime by (iii);
* bad primes come from resultants and spreading out over Knudsen's smooth model;
* `p`-adic trees are locally constant at `p ∈ Sigma`;
* `|F_S(u,v)| ≍ X^{e_S}` in a real cone;
* at good primes at most one `F_S(u,v)` is divisible, which gives a single-edge tree;
* anchoring by Theorem 1 puts the bounded `Sigma`-parts into `t_ij`;
* `Z` is avoided because `f(P^1) ⊄ Z`.

Hypothesis (ii) (birationality) is only used to get infinitely many distinct `Phi`. Since
`log n_S -> ∞`, that also follows without (ii).

### 2.2 (b) is a special case of (a): Remark 2 is false

**Claim.** If `Psi : (P^1_Q)^R ⇢ Mbar_{0,7}` satisfies the hypotheses of 2(b), then the hypotheses of
2(a) hold.

*Proof.*

1. **Choose weights.** Choose integers `a_k >= 1` and `q` with
   `|(sum a_k gamma_k).D_S - q| <= eps q`; this is the first step of L7's own proof.
2. **The curves.** Let `G = prod_k Rat_{a_k}`, the space of maps
   `g = (g_1, …, g_R) : P^1 -> (P^1)^R` with `deg g_k = a_k`. It is a Q-rational variety, so
   `G(Q)` is Zariski dense.
3. **Avoiding the bad locus.** Let `V ⊂ (P^1)^R` be the union of the indeterminacy locus of `Psi`
   and the sets `{F_S = F_{S'} = 0}`, `S != S'`. By hypothesis `V` has codimension at least 2.
   Evaluation `G × P^1 -> (P^1)^R` is smooth, so the `g` whose image meets `V` lie in a proper
   closed subset of `G`.
4. **Class.** For `g` outside it, `Psi∘g : P^1 -> Mbar_{0,7}` is a morphism and
   `deg (Psi∘g)^*D_S = sum_k a_k deg_k F_S = (sum_k a_k gamma_k).D_S`. So its class is
   `eps`-close to `q beta`.
5. **Conditions (i) and (iii).** No point of `P^1` maps into two boundary divisors, and a general
   `g` does not map into the boundary. Pass to the normalization of the image if birationality is
   wanted; this divides the class by the degree.
6. **Density.** If the images of the `g ∈ G(Q)` all lay in a proper closed `W`, then
   `{g : g(P^1) ⊂ Psi^{-1}W}` would be a closed set containing the dense `G(Q)`, hence all of `G`,
   and `ev` would not be dominant. So the images are Zariski dense. `[]`

Consequences:

* **Remark 2 is false.** "(b) is genuinely weaker than (a). No individual curve of class `∝ beta` is
  needed: mixtures of unbalanced families suffice." The hypothesis of (b) implies that of (a): a
  product structure **automatically** supplies the Q-rational, approximately balanced curves that
  (a) asks for. The geometric-sieve proof of (b) is correct but unnecessary.
* For the same reason, the summary's "(b) … No single balanced curve is needed" and the verdict's
  "or a product structure as in Theorem 2(b)" do not describe a separate route. A product
  structure is one particular way of producing the dense Q-rational curves of (a).

### 2.3 Repairs inside the proof of (b) (minor)

* (i) **`mu` is not a morphism.** The multiplication map
  `((u_1:v_1),…,(u_a:v_a)) ↦ (prod u : prod v)` is undefined where some `u_i = 0` and some
  `v_j = 0` (`i ≠ j`), a codimension-2 locus. So "`mu` is flat and surjective" is false as stated.
  * The needed conclusion still holds: no prime divisor of `(P^1)^{R'}` maps into
    `D_S ∩ D_{S'}`. Reason: a prime divisor `E'` not contained in the indeterminacy has
    `mu(E')` either dense or dense in a divisor `E`, by the dimension of the `(R'-R)`-dimensional
    fibres; and `Psi` is defined at a general point of `E`.
  * Alternatively, use weighted boxes `|u_k|, |v_k| ≈ r^{a_k}` and drop `mu` altogether.
* (ii) **Small primes.** For finitely many small `p`, `Y(F_p)` can be all of `F_p^{2R'}`. Such `p`
  must be put into `Sigma_0` and handled by the `p`-adic local condition. The local density
  `1 - O(p^{-2})` is positive only for `p` large.

### 2.4 (c)

(c) is correctly labelled as a sketch. The changes from L8 B.5 are reduced fibres, Galois orbits of
contact points, and archimedean recurrence away from `D(C)`. They look routine, but they are not
written out.

---------------------------------------------------------------------------------------------

## 3. Proposition 3 and the cubic class

### 3.1 Exact verification (`classes7.py` re-run, `cubic_check.py`, `movable_test.py`)

**Divisor identities.**

* `K = -(1/3) sum_{|S|=2} D_S` (coefficient `s(7-s)/6 - 2`).
* `sum psi_i = K + 2B`.
* `(B+8K).gamma = triple - (5/3) pair`.
* On `beta`: `-K.beta = 7`, `B.beta = 56`, `psi_i.beta = 15`, `(K+2B).beta = 105`,
  `(B+8K).beta = 0`. All correct.

**(i) Covering families have `-K.gamma >= 2`.** This holds by Kollár II.3.11: the general member is
free, `g^*T_X` is nef, and it contains `T_{P^1} = O(2)`. Correct.

**(ii) and (iii)** follow from the vanishing of the linear forms `B + 8K` and `7 psi_i + 15K` on
`beta`, together with integrality of `psi_i.gamma`. Correct.

**The cubic class `cub13`.** It is `e = 3` and meets the 13 planes `123,124,125,126,134,135,136,145,246,256,346,356,456`.

* It satisfies Keel, has every `D >= 0`, `pair = 6`, `triple = 13`, `-K = 2`, `(B+8K) = 3`.
* Its `psi`-degrees are `(7,5,5,5,5,6,3)`. I computed them two ways, with two choices of the
  auxiliary pair in the `psi` formula.
* Every 4-set contains at most 3 met planes.
* The decomposition `beta = (35/13) sym(cub13) + (21/26) sym(f)` checks on both invariant
  coordinates.
* In integers: `26·5040·beta = 70 sum_{sigma ∈ S_7} sigma_*cub13 + 21 sum_sigma sigma_*f`.

**Necessary conditions for a covering class, none violated** (`movable_test.py`).

* `cub13` pairs nonnegatively with all `7 × 15` pulled-back Keel–Vermeire divisors (via the
  pushforwards to `Mbar_{0,6}`).
* Its 21 pushforwards to `Mbar_{0,5}` have self-intersection between 5 and 15, so `>= 0`.

**Sanity checks of the Keel–Vermeire formula used.**

* `beta_6.KV = 5`.
* The forgetful fibre has `KV` value 1.
* The value is independent of the label assignment within each pairing.

### 3.2 The cubic family exists and covers (NUMERICAL; `cub13_frame*.py`, `cub13_monodromy.py`, `verify_through_point.py`)

**The frame.** In the frame `1 -> 0`, `2 -> ∞`, `7 -> 1` the movers have degrees `(2,2,2,3)`:

```text
x3 = (s-z34)(s-z35)/((s-p3)(s-p36)),   x4 = (s-z4)(s-z34)/((s-p45)(s-p46)),
x5 = k5 (s-z5)(s-z35)/((s-p45)(s-p56)),   x6 = 1 + c s(s-1)/((s-p36)(s-p46)(s-p56)).
```

The contacts `o_123, o_124, o_125` are at `s = 0, 1, ∞`. The system has 15 unknowns: 4 normalization
equations, and 8 equations for the four non-anchor triple points
`{3,4,5}, {3,4,6}, {3,5,6}, {4,5,6}`.

**The direct searches.** These are recorded as failures:

* Plain Newton (`cub13_frame.py`) converges only to cancellation-degenerate solutions, with
  `x_m ≡ 1` and `c = 0`.
* The `P^4` formulations (`cubic13_exist.py`, `cubic13_barrier.py`) converge only to degenerate
  solutions or not at all.
* With barrier unknowns for the collapsing differences (`cub13_frame_b.py`), one start in 400 gives
  a **nondegenerate** solution.

**The nondegenerate solution.**

* It has exactly 19 contact points, pairwise distinct (minimum distance 0.17).
* Each point lies on exactly one boundary divisor, and these are exactly the 17 divisors predicted
  by the class, with `D_16` met at three points. Since `B.cub13 = 19`, every contact is
  transversal.
* The 12×15 Jacobian has singular values in `[0.046, 5.8]`. So the family is locally smooth of
  dimension 3.

**Covering.**

* Fix a point `X0` on this curve. The 16×16 system "class `cub13` through `X0`" is nonsingular
  there (`sv_min = 1.9e-3`), so the evaluation map is étale at it and the family dominates
  `Mbar_{0,7}`.
* Monodromy loops (with some path jumping, filtered out afterwards) produced **8 distinct
  nondegenerate curves of class `cub13` through the same general point `X0`**. Each was
  re-verified: residual `<= 9e-15`, `sv_min >= 6.6e-4`, and 19 distinct transversal contacts on the
  predicted divisors.
* So the curve through a general point is not unique, and Q-rationality of members is not
  automatic.

### 3.3 Consequence for the "missing ingredient"

Here are the inputs:

* the forgetful fibres `f_i` form covering families, hence free ones, over Q;
* `cub13` forms a covering family, hence free (NUMERICAL above, plus Kollár II.3.11);
* `beta` lies in the cone of their `S_7`-orbits, exactly.

Then a connected tree of `70·5040` cubic translates and `21·5040` fibre translates, glued at general
points, is a tree of free rational curves. By the standard smoothing of trees of free rational
curves (Kollár, *Rational Curves on Algebraic Varieties*, §II.7; Debarre, *Higher-Dimensional
Algebraic Geometry*, §4.2), such a tree deforms to a free rational curve. Its class is
`131040·beta`, and a general deformation avoids every fixed codimension-2 set (Kollár II.3.7). So:

> **Over `C`, free rational curves of class `q·beta` (`q = 131040`) cover `Mbar_{0,7}`**,
> conditional only on the numerically established `cub13` family.

Hence the hypotheses of Theorem 2(a) fail, if at all, only because of the arithmetic requirement:
a Zariski-dense set of **Q-rational** members. Classes are not the problem, and neither is the
geometry over `C`. Prop. 3's sentence "the real obstruction is the product structure, not the
numerics of classes", and the summary's "the missing piece is building a product structure", should
be replaced by "the missing piece is Q-density of rational curves of class `≈ q beta`". L7.md §6.3
already says this. A product structure is one sufficient way to obtain them (§2.2).

It also follows that L7's negative `k = 7` numerics, which concern class exactly `beta`, are not
evidence about the geometry relevant to Theorem 2(a). Any `q` will do there.

---------------------------------------------------------------------------------------------

## 4. Theorem 4 and Lemma 5

**(i)** `pi_j^*D'_{S'} = D_{S'} + D_{S' ∪ j}`. Both are stable, so `r_j [C'].D' = 2` and
`r_j ∈ {1,2}`. Constancy is excluded because the class of a forgetful fibre is not `beta`. Correct.

**(ii)** I re-derived the argument; it is the deck-transformation argument of item 112 §3
(`seven_row_quadratic_quotient_obstruction.md`), with these steps:

* `r_{jj'} | 4`, by pushing forward to `Mbar_{0,5}`; each `D''` pulls back to 4 stable divisors.
* `r_{jj'} = 4`, because `K_j ≠ K_{j'}`.
* `<tau_j, tau_{j'}> ⊂ Aut(C(t)/L)` has order at most 4, so it is the Klein four-group.
* There is `F_2`-independence for at most 4 labels, using 3 anchors outside them.

For rational curves the conclusion `s <= 2` follows because finite subgroups of `PGL_2(C)` have
2-rank at most 2. Correct. It is new as a statement for genus 0: the notes have only the genus-one
version.

**(iii)** The Klein normal form is correct:

* `V_4 ⊂ PGL_2` is unique up to conjugacy, and its normalizer `S_4` permutes the involutions;
* the fixed fields are `C(t^2)`, `C(t + 1/t)` and `C(t^2 + t^{-2})`;
* the degrees `4, 4, 2, 2` follow from `deg x_l = 8` (item 112 (4)).

**(iv) (minor gap).** The unramifiedness of the 25 fibres, and "a equals the cluster value at
exactly one preimage", need `f^*D_{S'} ≠ f^*D_{S' ∪ a}`, that is, that no contact lies on
`D_{S'} ∩ D_{S'∪a}`. That intersection is a nonempty codimension-2 stratum, since the splits are
nested. The proof uses this ("excluded when every contact lies on a single divisor"), but the
statement of Theorem 4, whose standing hypotheses are only birationality, meeting `M_{0,7}`, and
class `beta`, omits it. Add "with 56 distinct contact points" to (iv).

**Lemma 5** is correct as a necessary condition. The kernel of the `25 x 18` matrix contains
`(N, D)`, for one of the `2^25` choices of `eps`. The phrase "8 conditions" is a generic-rank count.
In the summary, "in every case the 7th point lies in the kernel" needs two things: relabel so that
`r_7 = 1` (possible, since at most two labels have `r = 2`), and assume distinct contacts.

**Parameter counts and Corollary 6.** The table rows are consistent with each other. In particular
`k = 6`: one involution `4 + 2 - 1 = 5`, and Klein via (iv) and Lemma 5 `5 + 1 - 6 = 0`. They are
nonetheless heuristic. Two points matter:

* (i) L7.md §5 labels Corollary 6 "(heuristic)". The summary and §0.4 instead draw the conclusion
  "so the Galois–Boolean mechanism that refuted `(L_6)` cannot by itself refute `(L_7)`" (§0.4 also
  has "cannot by itself be dense at k = 7" in bold). That conclusion follows only from the
  parameter counts. It should read "is not expected to".
* (ii) Theorem 4 and Corollary 6 concern class **exactly** `beta`. Theorem 2(a) allows `q·beta`
  and approximately balanced classes. For class `q beta` one gets `r_j | 2q`, and nothing in L7.md
  constrains forgetful symmetries there. §3.3 shows that the relevant curves over `C` are of class
  `q beta` with large `q`. So the mechanism is excluded heuristically only in its literal class-`beta`
  form.

---------------------------------------------------------------------------------------------

## 5. Numerics and bookkeeping (reproduced)

* **`klein6_verify.py`.** Seed 3 has 25 fibres, each split once. All 10 `klein6t` solutions are
  nondegenerate.
* **`gen6.py`.** The residual is `9.3e-15` and the 17×28 Jacobian has rank 17. I also printed its
  singular values: the smallest is 0.252 and the largest 11.48, a robust gap with
  `h ∈ {1e-5, 1e-6, 1e-7}`. So the local dimension is 8 (NUMERICAL).
* **`klein3.py`.** The three logs have 283 + 284 + 277 = 844 runs, all failures, consistent with
  "0/844". The Klein ansatz's fixed anchor-fibre pattern (four `V_4`-orbits per anchor value) is
  consistent with the split counts. For example, `c = d = x` occurs at exactly 4 points, one
  `V_4`-orbit, and `a = b = c = d = x` occurs at exactly one point of that orbit. I did not prove
  that this pattern is the only one possible.
* **`realwalk.py 5/6/7`.** Closed walks are found after 11, 26 and 64 search nodes.
* **BMY.** I checked the universal-family numbers by hand: `e(S) = 4 + 56 = 60`,
  `K_S^2 = 8 - 56 = -48`, `K_S.sigma_i = 13`, `(K_S + sum sigma_i)^2 = -48 + 182 - 105 = 29`, and
  `3 e(S - sum sigma) = 3·46 = 138`. I did not check the orbifold expression
  `88 - 168/n + 224/n^2`. It is a "no obstruction found" remark and nothing depends on it.
* **`classes7.py`.** I re-ran it and checked the counts by hand:
  * `sum(|T|-1) = 15 + 40 + 45 + 24 = 124`;
  * `56 + 75 - 124 = 7` and `57 + 80 - 129 = 8`;
  * `56 - 80 = -24`;
  * the defect `10 + 15 + 6 = 31`.

  The interpretation "31 transitively implied congruences" is heuristic, which is acceptable.

---------------------------------------------------------------------------------------------

## 6. Hidden assumptions and quantifier slips (summary)

1. Theorem 4(iv) and Lemma 5 in the summary omit "56 distinct contacts" (minor).
2. Theorem 2(b) proof: `mu` is a rational map, not a flat morphism; small primes must be put into
   `Sigma_0` (minor; repairable, and (b) follows from (a) anyway).
3. Corollary 1 needs the top-block inflation under the all-blocks-balanced reading (minor).
4. "So `(L_7)` is purely a statement about the finite boundary heights": only one direction is
   proved (minor).
5. "Cannot by itself refute `(L_7)`": a heuristic conclusion stated as a consequence, and limited
   to class `beta` (major, framing).
6. "(b) genuinely weaker than (a)" and "the missing piece is a product structure": false, see §2.2
   and §3.3 (major, framing).
7. "Theorem 1 … new sharpening": already in `ptolemy.md` §5.2 (major, novelty).

None of these affects the truth of a stated theorem.

---------------------------------------------------------------------------------------------

## 7. Novelty

* **Not new:**
  * Theorem 1 (`ptolemy.md` Lemma 5.2 + §5.2, lines 492–494; L8 B.5 normalization);
  * the vanishing of `-K` on non-pair divisors and `K.beta = -7` (`ptolemy.md` §4, L8 table);
  * Theorem 4(i), the divisibility `r_S | 2^{|S|}` (item 112 (7));
  * the deck-involution argument of Theorem 4(ii) (item 112 §3, there for genus one).
* **New** (grep of the 1170 notes and the `claude_review` write-ups finds no prior statement):
  * Corollary 1 as a residue-free reformulation;
  * Theorem 2(a)/(b) as refutation criteria;
  * Proposition 3 and the `cub13` decomposition of `beta`;
  * Theorem 4(ii)–(iv) for genus 0, including the Klein normal form;
  * Lemma 5;
  * the `k = 6` genus-zero Klein solutions. Genus-zero balanced six-curves were already found
    numerically in `ptolemy.md` §5.4 and L8 §6.2; the Klein symmetric ones and the rank-8 check
    are new.

---------------------------------------------------------------------------------------------

## 8. Corrected statement

> **`(L_7)` remains undecided.** Proved (L7.md, with the repairs above):
>
> 1. *(Known: `ptolemy.md` Lemma 5.2 and §5.2.)* With the anchor and two labels at `∞, 0, 1`, every
>    `Phi ∈ M_{0,7}(Q)` has `Y_j = ±1` and `det(P_i,P_j) = ±prod_{T ∋ i,j} N_T(Phi)`. Hence
>    `(L_7) ⇒ (L_7)*` (Corollary 1, new formulation, one direction).
> 2. Theorem 2(a): a Zariski-dense set of Q-rational curves of class approximately `∝ beta`, with no
>    point on two boundary divisors, refutes `(L_7)`. Theorem 2(b), a dominant product chart with
>    `beta` in the cone of its ruling classes, is a valid criterion but a special case of (a).
>    Theorem 2(c) is a sketch.
> 3. Proposition 3 gives necessary linear conditions on product charts. They do not obstruct:
>    NUMERICALLY the class `cub13` is represented by a covering 3-dimensional family, with at least
>    8 members through a general point. Together with forgetful fibres and smoothing of trees of
>    free curves, this gives free rational curves of class `131040·beta` covering `Mbar_{0,7}`
>    over `C`. **The only missing ingredient for a refutation via (a) is Zariski density of
>    Q-rational members.**
> 4. For rational curves of class exactly `beta`: `r_j ∈ {1,2}`, and at most two labels have
>    `r_j = 2`. With two, the Klein normal form holds. With one and 56 distinct contacts, the
>    six-point image has class `beta_6` and the 7th coordinate solves the `25 x 18` system of
>    Lemma 5. Parameter counts (heuristic, class `beta` only) suggest that this
>    involution-based mechanism gives families of dimension at most 2 (genus 0) or 1 (genus one
>    with three involutions). That is too small to be dense, but it says nothing about classes
>    `q beta`, `q >= 2`.

---------------------------------------------------------------------------------------------

## 9. Files (`scratchpad/round4/referee_L7_checks/`)

| file | content | type |
|---|---|---|
| `thm1_cluster.py` | Theorem 1 via the cluster picture (independent of the min-formula): families A, B (deep clusters, non-admissible trees) and C (no normalization, `delta_p > 0` appears) | exact |
| `adm_recount.py` | recount of non-admissible trees in `anchor_check.py` with an exact two-chain test (568, agrees) | exact |
| `cubic_check.py` | `cub13`: planes, `psi`-degrees (two ways), linear tests, cone decomposition | exact |
| `movable_test.py` | `cub13` against all pulled-back Keel–Vermeire divisors; `Mbar_{0,5}` self-intersections | exact |
| `cub_frames.py` | `D`-vector of `cub13`; coordinate degrees in all 35 anchor frames (best `(2,2,2,3)`) | exact |
| `cub13lib.py`, `cub13_frame.py`, `cub13_frame_b.py` | the `(2,2,2,3)` frame; plain Newton (degenerate only); barrier Newton (nondegenerate solution) | NUMERICAL |
| `cub13_monodromy.py`, `verify_through_point.py`, `inspect_mono.py` | curves through a fixed general point: monodromy, re-verification (8 distinct nondegenerate) | NUMERICAL |
| `cubic13_exist.py`, `cubic13_barrier.py`, `cub13_count.py` | failed formulations (`P^4` model; fixed-`x3` count), kept for the record | NUMERICAL |
| `*.npy`, `monodromy_7.log` | saved solutions and logs | — |

Run each script from `referee_L7_checks/` with `python3`. The exact scripts take under 30 s, and
`thm1_cluster.py` about 22 s. `cub13_frame_b.py` and the monodromy runs take minutes.
