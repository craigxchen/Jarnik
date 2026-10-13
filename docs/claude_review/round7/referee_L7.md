# Referee report on `round7/L7.md` (second attempt at (L_7))

Referee scripts: `round7/referee_L7_checks/` (library `reflib.py`, tests `t1`–`t10`, logs `*.log`).
All of them were written from the definitions in `L7.md` alone. They do **not** import any author code
(`L7_checks/r8/*`). Exact arithmetic throughout: Python `Fraction`, `sympy` rationals, and
`flint.fmpz_poly` over `Z[s]`.

Classes of statements used below: **Proof** (checked line by line), **Data** (exact finite
computation), **Evidence** (exact computation that supports, but does not prove, a statement),
**Heuristic**.

---------------------------------------------------------------------------------------------

## 0. Verdict

**The computational core survives. The headline barrier is overstated.**

* **Confirmed by independent exact computation:**
  * Theorem A (the proof is correct);
  * Theorem B (the class table, with one gap in the proof that I close here);
  * Proposition 5.4, `lambda_A(beta) <= -23` and `lambda_A(mu) <= -10`.

  For Proposition 5.4 my own code gives:
  * 1260 distinct classes `sigma_* gamma_A`;
  * `sum comp = -9576`, average `-38/5`, range `[-12, 2]`;
  * `comp(mu) = -10`.

  The validation against Lemma 4.1 has 0 mismatches over 8820 (label, class) pairs. The same
  per-class values come out for three independent members of `V_A`.
* **Not proved as stated (major, scope).** The Verdict, Section 3.3, the Consequences after
  Proposition 5.4, and Section 9.3 say that `V_A` "cannot be the last family of **any** addition or
  section construction of approximately balanced curves", and that r7a Theorem 3 "does not apply to
  `V_A`, in any ordering". Corollary 4.5, which is the only tool, proves this only for curves `C'`
  that
  * (i) avoid `Sing(F_A)`, and
  * (ii) are **moving** sections, or sections meeting every fibre in the same component as the
    base section.

  Hypothesis (i) is **not automatic for good curves**. I exhibit a point `x_m in Sing(F_A)` that
  lies on exactly one boundary divisor, `D_{O b1 b2}` (Section 6; Data plus Lemma 5.3). For curves
  through `Sing(F_A)`, Lemma 4.2 bounds `Sigma^2` only by `comp(C')`, which is larger than
  `lambda_A(C')`. The dropped hypothesis matters for r7a's (G1) as well.
* **Proposition 5.2 is not proved (major for Section 9).** Its proof applies Lemma 4.2 to images of
  general ruling lines without checking that they avoid `Sing(F_{V_k})`. It also uses
  `lambda_{V_k}(gamma_k) = 2`. That equality fails for birational product charts whose leaves pass
  through a singular point; a 2-dimensional analogue gives `lambda = 1`.
  * The Section 9.1 table describes round 4 Theorem 2(b) as "birational product chart".
  * Theorem 2(b) actually allows any dominant `(P^1)^R ⇢ X`.
  * For `R > 4`, or for non-birational charts, `lambda_V` is not even defined. So "no candidate"
    under 2(b) is not supported.
* **Minor issues.**
  * r7a Theorem 3 is mischaracterised.
  * The parenthesis "(as every density argument needs)" in Corollary 4.5 is false.
  * Hypothesis (H) is never stated formally, and "exactly when" in Proposition 4.4 is overstated.
  * The model surface in Proposition 4.4 is the wrong one.
  * Theorem B does not address whether the parameter space is irreducible (closed here).
  * The novelty of Theorem A is not checked.
  * The conditional statement of Section 9.1 skips the twist and field-of-moduli issue.
* **Correctly labelled as evidence, and reproduced exactly:**
  * equality `lambda_A(beta) = -23`;
  * the list of `mu_T` differences;
  * the other event sets (`-31`, `-32`).

  Reproducing them exactly, over `Z` where the author used mod `p`, does not make them proofs: the
  decomposition `beta = sum_T mu_T - 32 mu` has a negative coefficient.
* **Growth and actual circles.** L7.md claims nothing for actual circles, and I agree that nothing
  changes. The proved bound `M <= (2/3 + eps) log R/loglog R` and the conditional Vojta theorem are
  untouched. `(L_7)` remains undecided.

`result_survives = true`, with the corrected statement of Section 9.

---------------------------------------------------------------------------------------------

## 1. What was recomputed independently

| claim in L7.md | referee check | result |
|---|---|---|
| §3.1: (A1); exactly three coaxial triples; three points on five planes; eight cube points; `dim Lambda = 4` | `t1_arrangement.py` | all confirmed |
| Thm B(2): member class table; 18 simple contacts; good | `t2_member.py`, members through `x = (6,-6,-3,-5,2,0,1)`, `(3,-2,4,0,2,-5,-6)`, `(1,7,-4,2,-9,5,3)`; the curve passes through `x` (all 35 cross-ratios at `s = 0` agree) | table reproduced exactly 3/3; `pair 6`, `triple 12`, `(B+8K) = 2` |
| `lambda_A(gamma_A) = 2` | Lemma 5.3 count on the leaf itself | `comp = 2` (3/3) |
| `comp(mu) = -10` | `t4_mu.py`, 4 random filtered integer Möbius curves, class `mu` verified | `-10` (4/4) |
| Prop 5.4: 1260 classes, `sum = -9576` | `t3_bound.py` (my own enumeration of `S_7`, dedup by class vector) | 1260 classes, `sum = -9576`, average `-38/5`, min `-12`, max `2`; histogram identical to the author's `lamA_bound_part0_of3` partial counts |
| validation by Lemma 4.1 | the same run: Lemma 5.3 count for the 7 forgetful foliations on every one of the 1260 test curves, against the closed formula `-psi_i + sum_j D_ij` | 0 mismatches / 8820 |
| equality evidence: a second member gives identical values (author: 60 classes) | `t3_bound.py` with two further members, **all 1260 classes** | per-class identical for all three members |
| `mu_T` differences, `3(-10) + 7 = -23` | `t8_muT.py`, 2 instances per `T` | list reproduced exactly; sum of differences `7` |
| other event sets `-31`, `-32` (author: mod p) | `t10_variants.py`, exact over `Z` | `-23, -31, -32, -32`; `comp(mu) = -10, -14, -16, -16` |
| `lambda_{f_i}(tau gamma_A) in {-2,-3,-4}`, `lambda_{f_i}(beta) = -9` | `t9_misc.py` | confirmed |

*How `comp` is computed in my code.*
* The tangent to the leaf at `R = span(p,q)` is obtained from first-order deformations of
  `R` inside the web quadric `B`: `B(p,w1) = B(q,w2) = B(p,w2) + B(q,w1) = 0`, with `w1, w2` in a
  fixed complement. This is not the author's ruling parametrisation.
* Orders at finite `t` come from one `gcd` and one `lcm` of the 35 reduced fractions `W_k/g_k^2`.
* The order at infinity is `2 max(deg N_k, deg D_k) - deg W_k`.
* `comp` does not change when `v` is multiplied by a rational function, so degenerate frames
  along the curve do no harm.

---------------------------------------------------------------------------------------------

## 2. Theorem A (Grassmannian model): proof correct; novelty unchecked (minor)

I re-derived steps 1–5.
* Rank 2 of the kernel uses (A1) correctly: a rank-1 `L = p mu^T` puts `p` on the six planes
  whose `v_a` are not killed by `mu`.
* Step 4 also needs `R ⊄ H_a`. This is implicit in `Phi_A(R)` being defined.
* The field-of-definition statement is correct, because signed maximal minors are polynomial in
  the entries.

The dictionary is also correct:
* a coaxial triple gives the Schubert divisor `sigma_1(l)`;
* the nine coaxial pair divisors are contracted (lines meeting `l` give a triple collision);
* `rho(X) - rho(G) = 42 - 1 = 41`.

Data: `t2_member.py` reconstructs `R_x` from `x` and recovers `x` exactly (all 35 cross-ratios).

*Novelty (minor).* The map is the restriction of the 7 linear forms `l_a` to 2-planes
`U ⊂ Q^4`. This is a point of `G(2,7)` (the 2-planes inside the fixed 4-space spanned by the
`l_a`), taken modulo the torus. So Theorem A is a special case of the Gelfand–MacPherson / Gale
correspondence `G(2,7)//T ~ (P^1)^7//PGL_2`, restricted to a sub-Grassmannian `G(2,4) ⊂ G(2,7)`.
* That this restriction is birational is easy and may be folklore.
* L7.md labels it "new" without any literature search.
* A quick web search found nothing either way.
* It should be called "new to this project; literature not checked".

## 3. Theorem B (the family `V_A`): confirmed; one gap closed

* **Class.** Data, `t2_member.py`, own code: the table of Theorem B(2) holds for members through
  three different rational points. All 18 contacts are simple, single-divisor and at finite
  parameters.
* **Uniqueness of the member through a general point.** The proof is correct. The closure step
  works because `Psi_A` is continuous at `x in U_B`, so `R_x in rho`.
* **Gap (minor), now closed.** Part 2 deduces "the class is constant on the family" from an open
  condition that holds at sample points. That inference needs the parameter space
  `{(Q, rho)}` to be **irreducible**: a double cover of `P(Lambda)` that could split into one
  component per ruling. L7.md does not address this.
  * Data (`t9_misc.py`): `det(sum c_i Q_i) = 4900 · (linear form) · (cubic form)`, squarefree.
  * So the branch divisor is reduced and nonempty, and the double cover is irreducible.
  * (Sampling `R_x` at random points would in any case land in both rulings, if there were two
    components.)
* **"First known rich addable Q-family."** This is correct within the project. Through any `Q`-point
  of `U_B` the member is a good `Q`-curve of class `gamma_A`, with `(B+8K) = +2`, so rich good
  `Q`-curves are Zariski dense. This is a genuine positive result.

---------------------------------------------------------------------------------------------

## 4. Section 4 (the corrected addability criterion)

### 4.1 Lemma 4.2: correct as stated; the "lambda(gamma_V) = 2" remark is not general

**Lemma 4.2.** Re-derived; correct under its hypotheses: `C ∩ Sing(F_V) = ∅`, `iota` defined along
`Sigma`, fibres mapped into leaves. Along `Sigma` the vertical tangent maps into `f^*F_V`, which is
saturated, and the map is nonzero.

**What happens without the hypothesis `C ∩ Sing = ∅`.** This is the point that matters below.
* The vertical tangent `T_{Y/P^1}|_Sigma` maps into the **saturation** of `f^*F_V` in `f^*T_X`.
* The degree of that saturation is exactly the count `comp(f, v)` of Lemma 5.3. It equals
  `lambda_V(C) + (excess at the points of C ∩ Sing(F_V))`.
* So in general one only has `Sigma^2 <= comp(C)`.

**Section 4.1's claim `lambda_V(gamma_V) = 2`.** It is justified by "a general member is a smooth
leaf avoiding `Sing(F_V)`". That is an assumption, not a fact about families with one member
through a general point.
* Example: the lines through a point `A` of `P^2`. This is a family with one member through a
  general point, but every leaf passes through the radial singularity `A`. Here `T_F = O(1)`, so
  `lambda(gamma) = 1`.
* For `V_A` the value 2 is correct (Data: `comp = 2` on the leaf, and `comp` is an upper bound,
  while `deg f^*F >= deg T_{P^1} = 2` on a leaf).
* For forgetful fibres the value 2 is also correct.

### 4.2 Corollary 4.3 and r7a Theorem 3: mischaracterised (minor)

Corollary 4.3 is correct: (G1) together with `C ∩ Sing(F_V) = ∅` gives `lambda_V(C) >= 2 delta >= 0`.

But L7.md's account of what this does to r7a is off:
* §0 says r7a Theorem 3 "is not correct" and that "its proof uses that [(G1) is open]".
* §4 is headed "r7a Theorem 3 is not correct".

**r7a Theorem 3** (`L7_prev_attempt1.md` §4.3) **includes** the hypothesis "the curves of
`F_{j-1}` satisfying (G1)–(G3) … have Zariski-dense union". Its proof uses that hypothesis, not
openness.
* What Corollary 4.3 refutes is r7a's *Remark 1*: "(G1)–(G3) fail only when a curve meets a fixed
  codimension-2 locus or is itself special".
* Corrected: r7a Theorem 3 is a conditional statement whose density hypothesis is not known to hold
  for any family. By Corollary 4.3, it can hold at a late step only through curves that meet
  `Sing(F_V)`.
* That loophole is real for `V_A` (Section 6 below). So "r7a Theorem 3 does not apply to `V_A`, in
  any ordering" (end of Corollary 4.5) is **not proved**. What is proved: no curve avoiding
  `Sing(F_A)` with `lambda_A(C) < 0` satisfies (G1).
* For the forgetful families the loophole is absent. `Sing(F_i)` is the set of nodes of fibres,
  which lies in `∂^{(2)}`, and good curves avoid it.

### 4.3 Proposition 4.4: correct in substance; three wording issues (minor)

The Zariski-lemma argument is correct. So is the base-point-freeness of `|Sigma + mF|` for
`m >= max(1, -Sigma^2)`. I also checked:
* `R^1 p_* O_Y = 0` for tree fibres;
* the restriction sequence;
* `(Sigma + mF)^2 = Sigma^2 + 2m > 0`, so the general member is irreducible.

Issues:
* **Wrong model surface.** "Resolve the normalisation of the closure" of the *swept surface* need
  not fibre over `C`, because members through different points of `C` may meet. The right object
  is a resolution of the incidence surface `{(c, y) : y in M_c}`.
* **(H) is never stated as a single hypothesis.** §5.1 lists three bullets. For the conclusion
  "the image is a good curve" one also needs `iota` to be birational onto its image along `sigma'`
  (r7a's (G3)), which is not among them.
* **"Adding V to C is possible exactly when m >= -lambda_V(C)" is overstated.**
  * Necessity holds only for sections meeting every fibre in the same component as `Sigma`.
  * Sufficiency is conditional on (H).

### 4.4 Corollary 4.5: the parenthesis is false; the headline drops hypotheses (major)

The inequality `lambda_V(C') >= sigma'^2 >= 0` is correct **for a moving section `sigma'` with
`C' ∩ Sing(F_V) = ∅`**.

* **The parenthesis "(as every density argument needs)" is false.**
  * Every curve `C'` that avoids `Sing(F_V)` and is transverse to `F_V` is the base section of its
    own leaf surface, with `Sigma^2 = lambda_V(C') < 0`.
  * So "`C'` is a section of a `V`-ruled surface" is no restriction at all.
  * A density argument can let the base curve vary over a dense family and take **rigid** sections.
    Corollary 4.5 says nothing about those.
* **What survives without "moving".** For sections in the same component (Proposition 4.4 type),
  one has `sigma' ≡ Sigma + mF` with `m >= -lambda_V(C)`. Hence
  `sigma'^2 = lambda_V(C) + 2m >= m >= 0`, and Lemma 4.2 again gives `lambda_V(C') >= 0`.
  So **contact-preserving additions** are covered without the word "moving". Sections that change
  component are covered only if they move.
* **The headline drops both hypotheses.** L7.md §0, §3.3, the Consequences after Proposition 5.4,
  and §9.3 drop both "avoids `Sing(F_A)`" and "moving / same component" when they state that
  `V_A` "cannot be the last family of **any** addition or section construction".
  * The moving / same-component hypothesis can simply be restored.
  * The `Sing` hypothesis is not automatic (Section 6).

---------------------------------------------------------------------------------------------

## 5. Lemma 5.3 and Proposition 5.4: correct, and independently confirmed

**Lemma 5.3.** Re-derived; correct.
* `f^*F` is a line bundle. Since `f(P^1) ⊄ Sing(F)`, `f^*F -> f^*T_X` is injective.
* The orders as `T_X`-sections are at least the orders as `F`-sections, with strict inequality
  exactly at points of `Sing(F)`.
* Summing over `P^1` (including `t = infinity`) gives `deg f^*F <= comp(f,v)`, independently of the
  rational section `v` chosen.
* The local immersion claim holds. It also follows from Keel's closed embedding
  `X -> prod Mbar_{0,4}`, so the "at most one boundary divisor" hypothesis is not even needed for
  step 1.

**Proposition 5.4.** Every input checks out.
* *Decomposition.* `beta = (35/12) avg_{S_7} sigma_*gamma_A + (1/12) mu`. The pair and triple
  degrees check, and the `D_S` span `N^1(X)`.
* *Stabiliser and representatives.*
  * The stabiliser of `gamma_A` has order 4: `b1 <-> b2` and `(a1 c1)(a2 c2)`. Data: exactly 1260
    distinct classes.
  * Each relabelled curve `sigma ∘ f_x` is good of class `sigma_* gamma_A`, because `sigma` is an
    automorphism permuting the `D_S`.
  * Its image is not in `Sing(F_A)`. My code requires the 3x4 tangent system to have generic rank 3,
    and all 35 `W_k` were nonzero in every case.
* *Arithmetic.* `(35/12)(-38/5) + (1/12)(-10) = -133/6 - 5/6 = -23`.

So `lambda_A(beta) <= -23` and `lambda_A(mu) <= -10` are **proved**: Lemma 5.3 plus exact Data,
confirmed by an independent implementation (Section 1).

**Equality** (`= -23`) is evidence, correctly labelled. I add:
* per-class identity of all 1260 values for three different members;
* exact (over `Z`) reproduction of the `mu_T` table.

Both are consistent with the test curves avoiding `Sing(F_A)`. Neither proves it.

**Section 5.3, items 3–4 (local analysis).** These are correctly marked as relying on an assumption
("`Psi_A` contracts no interior divisor") and on evidence (`d_R`). They are not used in the proof of
Proposition 5.4.
* One small slip: "the six event divisors have `c = +1`: from `2 = -2·2 + 6c`" assumes that the six
  values are equal. They are not related by the order-4 stabiliser.
* The local radial computation gives `+1` anyway, so the conclusion stands.

---------------------------------------------------------------------------------------------

## 6. `Sing(F_A)` meets the boundary away from `∂^{(2)}`: the dropped hypothesis is substantive

L7.md never determines `Sing(F_A)` on `X`, and its headline statements silently assume that good
curves avoid it. They need not.

### 6.1 The non-unique-quadric locus (Data, `t5_degen.py`)

The lines `R` with `rank(Lambda|_R) <= 2` (more than one web quadric through `R`) were computed
exactly, on random 2-plane sections of a chart of `G(2,4)`. Each section gives 10 rational points.
Classifying them, the locus is

```text
6 sigma_2(P_e)  +  sigma_2(v_13)  +  sigma_11(H_a1)  +  sigma_11(H_c1)  +  sigma_11(Pi),
class 7 sigma_2 + 3 sigma_11.
```

* This matches the Thom–Porteous class `c_1^2 - c_2` of `Sym^2 S^*`.
* Here `v_13 = l_1 ∩ l_3`, and `Pi` is the plane spanned by `l_{a1c2}` and `l_{a2c1}`.
* All components lie in the boundary. So **in the interior `M_{0,7}`, non-uniqueness of the quadric
  causes no singularity**.

### 6.2 A singular point on a single boundary divisor (Data)

The web contains rank-2 quadrics `P1 ∪ P2` coming from partitions of `E` into two triples
(`t6_sing_probe.py`). Take

```text
P1 = < P_{a1b1c1}, P_{a1b2c2}, P_{a2b2c1} >,   P2 = < P_{a1b1c2}, P_{a1b2c1}, P_{a2b1c1} >,   m = P1 ∩ P2.
```

**Where `m` sits.**
* `m` meets the axis `l_2` (labels `O, b1, b2` collide) and nothing else.
* `m` lies in none of the seven planes.
* So its point `x_m` lies on **`D_{O b1 b2}` only**, not in `∂^{(2)}`.
* The web quadric through `m` is unique and equals `P1 ∪ P2`. Since `m = P1 ∩ P2`, it has no
  "ruling through `m`".

**Test 1 (Lemma 5.3, in `X`; `t7_sing_m.py`, two seeds).**
* Pencils of lines `(y, Pi')` with `y in m` and `Pi' ⊃ m` are Möbius curves through `x_m`.
* Each is a good curve whose class is exactly that of a generic pencil (three coaxial triples and
  twelve non-coaxial pairs, once each). The contact at `x_m` is a simple contact with
  `D_{O b1 b2}` alone.
* Results:
  * pencils through `m`: `comp_A = -1` in 21 of 21 cases;
  * generic pencils: `comp_A = -2` in 56 of 56 cases.
* Lemma 5.3 gives `comp >= lambda(class)` with equality iff the curve avoids `Sing(F_A)`, and the
  classes are equal. So **every tested good curve through `x_m` meets `Sing(F_A)`, with excess 1**.

**Test 2 (direction field in `G`; `t11_dirlimit.py`, control `t11_control.py`).**
* Along six random pencils through `m`, the limit as `s -> 0` of the tangent line of the leaf
  (the Plücker point of `span(pi, pi')`) gives **six distinct limits**.
* Control: at a generic line, and at a generic line meeting `l_2`, five pencils give one limit.
* So the line field of `F_G` has no continuous extension at `m`, i.e. `m in Sing(F_G)`.

I did not prove separately that `Phi_A` is a local isomorphism at `x_m`; in the local coordinates of
`D_{O b1 b2}` this is clear. Together with Test 1 the conclusion is
**`x_m in Sing(F_A) \ ∂^{(2)}` (Data)**.

### 6.3 Consequences

* For a good curve `C'` through `x_m`, or through any other point of `Sing(F_A) \ ∂^{(2)}`,
  Lemma 4.2 gives only `sigma'^2 <= comp(C') = lambda_A(C') + excess`.
* So Corollary 4.5 does not exclude a moving `V_A`-section near `L beta` that passes through
  `Sing(F_A)` often enough to make up the deficit `23 L`.
* The same loophole affects r7a's (G1) (Section 4.2).
* For forgetful fibres there is no loophole, since `Sing(F_i) ⊂ ∂^{(2)}`.

*Heuristic only (not proved): the loophole looks unexploitable.* Rational curves of class `L beta`
have expected dimension `-K.(L beta) + 1 = 7L + 1`. Each passage through an isolated singular point
costs 3 conditions and gives excess about 1. Even if `Sing(F_A)` had surface components inside
boundary divisors, at 1 condition per passage, the total excess would be `≲ 7L < 23L`, unless
passages of excess `>= 4` were cheap. This does not replace a proof: constructed families need not
have the expected dimension, and the full `Sing(F_A)` is not known.

**Required correction.** Either keep the hypothesis "`C'` avoids `Sing(F_A)`" in every statement
derived from Corollary 4.5 (§0 items 3–4, §3.3, the Consequences after Proposition 5.4, §9.1,
§9.3), or prove a bound on the total excess of approximately balanced good curves.

---------------------------------------------------------------------------------------------

## 7. Proposition 5.2 (product charts): proof gap, and the scope of Theorem 2(b) is misstated

**What the proof assumes.** The proof uses two facts:
* (a) `lambda_{V_k}(gamma_k) = 2`;
* (b) `lambda_{V_k}(gamma_j) >= 0` for `j != k`, by Lemma 4.2 applied to the image of a general
  `t_j`-line as a section of a `P^1 x P^1`.

Both need the curves involved to avoid `Sing(F_{V_k})`: the leaves for (a), the images of
`t_j`-lines for (b). Nothing in the hypotheses ensures this.

**Correct general formula.** Let `Z` be a common resolution, `p: Z -> (P^1)^4`, `q: Z -> X`. Then

```text
lambda_{V_k}(gamma_j) = 2·[j = k] + sum_{E'} b_{E'} (E'.l_j),
```

where `E'` runs over the divisors of `(P^1)^4` contracted by `Psi`, and `b_{E'}` is the order along
`E'` of the generator of the `t_k`-foliation, read as a section of `q^*F_X`.
* `b_{E'} >= 0` when `F_X` is regular at the generic point of `Psi(E')`.
* `b_{E'}` can be negative when `Psi(E') ⊂ Sing(F_{V_k})`.

**2-dimensional analogue (Proof).** Take `X = P^2`. Let `Psi(s,t)` be the intersection of the line
through `A` with slope `s` and the line through `B` with slope `t`. This is birational.
* The `t`-lines go to the lines through `A`. Their foliation is radial at `A`, with `T_F = O(1)`.
* The curve `{t = slope(BA)}` is contracted to `A`. Every leaf meets it.
* Hence `lambda_1(gamma_1) = 1`, not 2. In the formula above, `b = -1` (blow-up computation).
* With `beta := gamma_1`, `c_1 = 1`, `c_2 = 0`, the inequality `lambda_1(beta) >= 2c_1` fails.

Transporting this to `Mbar_{0,7}` would need a chart with `beta` in its cone, which I did not
construct. So the inequality of Proposition 5.2 is **unproved**, not refuted. The proof is
complete under the extra hypothesis that `Psi` contracts no divisor `E'` with `E'.l_j > 0` (for
some `j`) into `Sing(F_{V_k})`.

**Scope (major for §9).** Round 4 Theorem 2(b) is stated for **dominant**
`Psi: (P^1)^R ⇢ X`, any `R`, with no birationality.
* L7.md §9.1 calls it "birational product chart".
* §5.2 and Q1 then say a positive answer to Q1 "would show that no birational product chart
  refutes `(L_7)`".
* For `R > 4`, or for non-birational `Psi`, the ruling families have several members through a
  general point, and `lambda_V` is not defined.
* So the §9.1 entry "no candidate: forgetful `-9`, Kapranov `-26`, `V_A <= -23`" covers only the
  birational `R = 4` case, and that case only modulo the gap above.

---------------------------------------------------------------------------------------------

## 8. Other points

1. **Proposition 5.1** (sufficiency) is conditional on (H), and is clearly labelled as such. The
   induction is correct as written. It also needs r7a's (G3) (birationality onto the image) for
   goodness; see Section 4.3.
2. **Section 6 (Möbius curves).** Re-derived the formulas.
   * Contributions: `(B+8K) = 6 tau_3 + 11 tau_4 + 15 tau_5 - 70` (pure `mu`: `-70`; a 3-, 4- or
     5-event contributes `+6`, `+11`, `+15`).
   * Expected dimension `15 - tau_3 - 2 tau_4 - 3 tau_5`.
   * The only rich types of dimension `>= 3` are `(12,0)`, `(10,1)`, `(8,2)`, with no 5-events.
     Data: `t12_mobius_types.py`, an exhaustive search over `(tau_3, tau_4, tau_5)`.
   * `pair = 6`, `-K = 2`.
   * The table "16 of 23 types realised" and the monodromy degrees are NUMERICAL and labelled as
     such. Not rechecked.
3. **Proposition 7.1.** The proof is correct (`6 <= 9`, `11 <= 12`, `15 <= 15`; pair zeros `6D`;
   `-(5/3)(-9 n_5) = 15 n_5`). "Free-linear" is a definition, not a theorem about all
   constructions. The sentence "every rich family is cut out by nonlinear conditions in the graph
   presentation" is a paraphrase of the contrapositive, and should not be read more broadly.
4. **Section 8 (class `beta`).** The following counts are correct:
   * `d_a = 3`: each pair has 1 pure contact and 5 triple events, so `d_a + d_b = 6`;
   * `5 d_a + D = 36`;
   * expected dimension `49 - 6 - 35 = 8 = -K.beta + 1`;
   * the sextic surface in `P^6`.

   The continuation results are NUMERICAL, with honest stop reasons.
5. **Section 9.1, conditional statement.** "`M_k` rationally connected + Colliot-Thélène ⇒
   `M_k(Q)` dense ⇒ `(L_7)` false" skips steps:
   * a smooth proper model is needed;
   * a `Q`-point of the moduli space gives a genus-0 curve over `Q` that may be a pointless conic
     when `q` is even (field of moduli versus field of definition);
   * the goodness conditions must be shown to be open and nonempty on `M_k`.

   It is labelled conditional. The chain of implications should be written out, or marked as a
   sketch.
6. **Section 4.3 example.** The degrees `D_{ijk}.e = -1` and `D_{ij}.e = 1` and the class
   `[C] + m f_i + e_{ijk}` are correct.
7. **Attribution inside the project.**
   * "First rich addable Q-family": correct relative to round 4 and r7a.
   * The claim that cub13 "has at least 8 members through a general point" is cited from r7a. Not
     rechecked.

---------------------------------------------------------------------------------------------

## 9. Corrected statement

`(L_7)` is not decided. What this round establishes is the following.

1. **Theorem A (proved).** For seven planes in `P^3` over `Q` satisfying (A1) and (A2),
   `R -> (R ∩ H_a)_a` is a `Q`-birational map `G(2,4) ⇢ Mbar_{0,7}`, with inverse given by the
   kernel of a `7 x 8` linear system. It is the restriction of the Gale / Gelfand–MacPherson
   correspondence to a sub-Grassmannian; the literature has not been checked.
2. **Theorem B (proved; class by exact Data, irreducibility of the parameter space added by the
   referee).** For the star arrangement and the six cube points `E`, the rulings of the web
   quadrics form an irreducible 3-dimensional `Q`-family `V_A`.
   * Exactly one member passes through a general point, computed over its residue field.
   * The general member is good, with class `gamma_A`: `pair 6`, `triple 12`, `(B+8K) = +2`.
   * Rich good `Q`-curves are Zariski dense in `Mbar_{0,7}`.
3. **Corrections to r7a (proved).**
   * For a family `V` with one member through a general point and a good curve `C` **avoiding
     `Sing(F_V)`**, r7a's (G1) forces `lambda_V(C) >= 0` (Corollary 4.3). So r7a's Remark ("(G1)–(G3)
     fail only on a codimension-2 locus") is false.
   * r7a Theorem 3 survives only as a conditional statement whose density hypothesis can hold at
     late steps only through curves meeting `Sing(F_V)`.
4. **The invariant (proved).** Let `C'` be a good curve, the image of a section `sigma'` of a
   `V`-ruled surface along which `iota` is defined, with **`C' ∩ Sing(F_V) = ∅`**. If `sigma'`
   either moves or meets every fibre in the same component as the base section, then
   `lambda_V(C') >= 0`.
5. **Proposition 5.4 (proved; exact Data confirmed independently).**
   `lambda_A(beta) <= -23` and `lambda_A(mu) <= -10`.
   * Equality `= -23` is evidence: three independent members give identical per-class values,
     and the `mu_T` table is reproduced exactly.
   * The other event sets give `-31` and `-32`. This is evidence: the decomposition has a negative
     coefficient, and the classes of those families are unchecked.
6. **Consequence (proved).** No approximately balanced good curve **that avoids `Sing(F_A)`** is a
   moving (or same-component) section of a `V_A`-ruled surface. The same holds for every
   `S_7`-translate of `V_A` and for the forgetful families. For the forgetful families the
   avoidance is automatic for good curves.

   For `V_A` it is **not automatic**: `Sing(F_A)` contains the point `x_m` on `D_{O b1 b2}` alone
   (Section 6, Data). Whether passages through `Sing(F_A)` can make up the deficit `23L` is open.
   A dimension heuristic says no, but this is not a proof.
7. **Not proved.**
   * "`V_A` cannot be the last family of **any** addition or section construction". What is proved
     is the restricted statement 6.
   * "r7a Theorem 3 does not apply to `V_A` in any ordering". What is proved is the restriction to
     curves avoiding `Sing(F_A)`.
   * Proposition 5.2 as stated. Its proof needs the leaves and the images of general ruling lines to
     avoid `Sing(F_{V_k})`.
   * Anything about round 4 Theorem 2(b) beyond birational `R = 4` charts.

**Meaning.** Unchanged from L7.md §9.3, and I agree with it.
* Nothing here changes the uniform theorem or the growth rate.
* The proved bound `M <= (2/3 + eps) log R/loglog R` and the conditional Vojta theorem are
  untouched.
* No statement of L7.md involves inputs from actual circles. The objects are curves on
  `Mbar_{0,7}`, linked to circles only through the missing lemma `(L_7)` of the Ptolemy route.
* The result is a partial barrier, and its scope is smaller than the headline says: the curves
  that avoid `Sing(F_A)`, and moving or same-component sections.

---------------------------------------------------------------------------------------------

## 10. Files (`referee_L7_checks/`)

| file | content |
|---|---|
| `reflib.py` | the referee's exact library: arrangement, web, kernel line `R_x`, member through `x`, contact class and goodness, `F_A` tangent by first-order deformation in the web quadric, the Lemma 5.3 count, forgetful tangents |
| `t1_arrangement.py`, `t1.log` | (A1), coaxial triples, five-plane points, cube points, `dim Lambda = 4` |
| `t2_member.py`, `t2.log` | Theorem B class for three members; `comp = 2` on leaves |
| `t3_bound.py`, `t3_first.log`, `t3_second.log`, `t3_third.log`, `t3_vals_*.pkl` | Proposition 5.4: 1260 classes, sum `-9576` for three members; forgetful validation (0 mismatches over 8820) |
| `t4_mu.py`, `t4.log` | `comp(mu) = -10` |
| `t5_degen.py`, `t5.log` | the non-unique-quadric locus |
| `t6_sing_probe.py`, `t6.log` | the lines `m = P1 ∩ P2` for the 10 partitions of `E` |
| `t7_sing_m.py`, `t7.log`, `t7b.log` | pencils through `m`: `comp -1` (21/21) vs generic `-2` (56/56), same class |
| `t11_dirlimit.py`, `t11.log`, `t11_control.py`, `t11_control.log` | discontinuity of the leaf direction at `m` (6 distinct limits; control: 1) |
| `t8_muT.py`, `t8.log` | the `mu_T` table, sum `-23` |
| `t9_misc.py`, `t9.log` | `det` of the web quadric squarefree (irreducible double cover); `lambda_{f_i}` values |
| `t10_variants.py`, `t10.log` | other event sets `-23`, `-31`, `-32`, `-32`, exact |
| `t12_mobius_types.py`, `t12.log` | rich Möbius event types |
| `NOTES.md` | working notes |
