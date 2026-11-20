# Referee report rays1: Part X.2, first half

Scope: Part X.2 of `astra_free_X/answer.md`, from "## X.2." (line 2348) through the paragraph ending "This proves the assertions. ∎" (line 2393). That covers the statement, the coarse profile, the class encoding, the local construction with the positional digit-sum covering (X.1), the index bound (X.2), the angular completion and the count of words. The three paragraphs after line 2393 (the deliberate failure at l\* < 8P, the auxiliary heights (X.3) on the 5 mod 8 version, and the limit of the order argument) belong to block rays2. They are commented on only where item (4) requires it. I also read `x_finite_work.md`, `x_digit_check.py`, `final10.txt`, README_IX §§0, 1.3, 1.4, 1.7, 1.8, sharp.md §§1.1–2 and construct.md §1.1.

Referee scripts (my own code, exact integer arithmetic, nothing of Astra's imported) are in `ref10_rays1/`:
`digit_ref.py`, `gauss_ref.py`, `fraction_ref.py`. Outputs are quoted in §6.

## 0. Verdict

1. **The digit-sum covering (X.1) is correct.** I checked the proof line by line and found no error. Three of its inequalities are tight or need a short justification that she does not write out: the digit bound 16H+1 ≤ B−1, the closed interval [a_t/B, (a_t+2)/B], and the rounding loss below 4g. My own implementation passed every check: 320 blocks of 512 or 1,152 slots at r = 2, 3, 5, 7, 13, both with one-hot rows and with random detecting rows, 2,631,680 exact digit-interval checks, and 6,985,602 tested-word non-divisibility checks at H = 2, 3 (0 failures). I also reran her script, which printed her reported counts: 200 blocks, 1,152 slots, 921,600 intervals, 32,438 words. **Proved.**
2. **The index bound (X.2) and the collision bound are correct.** For every tested word, gcd(Σ v_j L_j, d_a) divides G = lcm(1, …, ⌈U⌉). Hence a collision of (∏ t_j^{v_j})^n with i^s at l^a forces d_a | 4nG and l^a ≤ 6nG, at every level and every cutoff. I checked this with actual Gaussian residues mod l^a for l = 3, 5, 7, 13 (8,068 tested words × 12 powers × 4 units, 0 failures). **Proved.**
3. **The object lies in F^char_X for every X, and it satisfies (G) and the stated restricted e-all at every X.** It does so with W/(m log m) → 12g/ρ = 1536 H²/ρ and m/s → ρ/(256 H²). **Proved**, with the corrections in §1.
4. **Single labels: every power of every label passes e-all at all four units and every cutoff.** **Proved.** This already holds with one class (b = 1) and needs no growing partition.
5. **Corrections.** (i) The words not tested are exactly those whose class-sum vector is zero. Every word of odd mass, and every word with nonzero total coefficient sum, is tested. (ii) "Counterexample to a growth argument requiring this fraction of the tests" holds only for arguments whose set of words is fixed independently of the residue data. An argument that picks its words from the data, as X.4 does (two labels with equal residues), is not blocked. (iii) The constant 1536H²/ρ grows with H, so no single fake covers words of unbounded length. (iv) The threshold q_0 is astronomically large (log q must exceed e^U with U ≥ 3·10⁹). This is a pure existence statement. **CORRECT WITH CORRECTIONS** overall; no step is false.

## 1. Corrected statement of X.2, first half (quoted in full; referee wording)

> **Proposition X.2 (restricted ray fake).** Fix C > 0, an integer H ≥ 2 and ρ ∈ (0.99, 1). Put B = 64H², g = 2B = 128H², d = 3. Fix any function b(q) → ∞ with U(q) := 100·g·B^{T(q)} = o(log log q), where T(q) = b(q) (one-hot rows) or T(q) = ⌈2H log₂(2b(q))⌉ + 2 (random detecting rows). Then there is q_0 = q_0(C, H, ρ, b(·)) such that for every prime q ≡ 3 mod 4, q ≥ q_0, there is a tuple (p, e, a, φ, θ_0, δ, k, r) with the following properties.
>
> * Profile: the g-fold replicated Paley profile of IX.3/IX.5 (column-zero deleted, a_ij = 1 + χ(j − i), e_j = 2), with s = g(q − 1) physical primes p_j ∈ [P, 2P], P = q³, in one sector and one radial bin, and m = (ρ/2 + o(1))q marks. W = Σ 2 log p_j satisfies W/(m log m) → 12g/ρ = 1536H²/ρ. The limit does not depend on C; the threshold does.
> * The tuple is in F^char_X(m, C) of sharp.md Def. 1.1 for **every** X ≥ 3, with the same data. In particular the same-level clause holds at every parent prime.
> * Every v ≠ 0 satisfies both Lemma 0.1 gaps (G).
> * The labels are partitioned into b = b(q) classes, balanced or with a bounded number of reserved singletons. For every v ∈ Z^s with 1 ≤ ‖v‖₁ ≤ H whose class-sum vector w(v) = (Σ_{cl(j)=κ} v_j)_κ is nonzero, every n ∈ Z∖{0}, every s ∈ Z/4 and every X, inequality (1.2) holds for nv. In particular this holds for every n e_j, for every e_j ± e_k with cl(j) ≠ cl(k), for every e_j + e_k, and for every word of odd mass ≤ H.
> * The words of mass h ≤ H that are not covered are exactly those with w(v) = 0. They make up a fraction O_H(1/b) + O_H(1/s) of all words of mass h, and none of them when h is odd. Among differences the uncovered fraction is (s/b − 1)/(s − 1) for balanced classes.
> * The angles φ_j can be chosen off the finite set {±arg π_j + kπ/2}, so they are nongenuine. Consequently, by IX.7, no label's residue family is cofinitely the reduction of one element of Q(i)*.
>
> The proposition does **not** assert e-all for the words with w(v) = 0, so the tuple is not claimed to lie in F^all_X. Any argument that uses only P3-char, (G) and e-all on a *data-independent* set of words of mass ≤ H, contained for some balanced partition in the set above, cannot prove W/(m log m) → ∞. It can prove at most a constant ≤ 1536H²/ρ.

Labels: **Proved** (all bullets).

## 2. The object, precisely (task item 1)

* **Profile.** Paley matrix K_ij = χ_q(j − i) on F_q, q ≡ 3 mod 4. Column 0 is deleted. Each of the q − 1 remaining columns is replicated g times; I call these g copies a *block* J. The entries are a_xj = 1 + K_xJ for every j ∈ J, with e_j = 2. The marks are a ρ-density sample of {x : K_x0 = 1}. The exponent row of a mark is constant on every block. Labels are physical primes p_j ∈ [P_0, P_0 e^σ] ⊂ [P, 2P], P = q³, in a sector of width b\*/g, b\* = min(1, C)/1000. There are s = g(q − 1) of them. This is IX.5 "details of the replicated profile", refereed Proved (README_IX §1.7).
* **Classes.** The b classes are an arbitrary balanced partition of the s labels (by index). They are unrelated to blocks, and a block may contain any mixture of classes. The detecting rows c_1, …, c_T ∈ {0,1}^b satisfy: every nonzero w ∈ Z^b with ‖w‖₁ ≤ H has some t with c_t·w ≠ 0. The random-row existence proof is correct: P(c·w = 0) ≤ 1/2 by exposing one used coordinate, there are at most 2(2b)^H vectors, and T ≥ H log₂(2b) + 2 suffices. *Simplification (Proved):* the b one-hot rows detect every nonzero w trivially. The random rows only improve the dependence of U on b, and that is irrelevant because b → ∞ arbitrarily slowly.
* **Coarse data (unchanged from IX.3/IX.5).** The block angles Φ_J = Σ_{j∈J} φ_j, the offsets δ, the windings k, θ_0 and the global witnesses r_j^old = i^{u_j} π_j ∏ γ_l^{V_jl} on the replicated matrices. Only Φ enters the realisation equation, because ⟨2a_x − e, φ⟩ = Σ_J (2a_xJ − 2)Φ_J. The coarse noncharacter gaps come from the Haar completion on the nullspace (IX.5), which needs d > 2; here d = 3.
* **Exterior residues (new, genuinely local).** At each odd l ∤ N and each level a, T_{l^a} is cyclic of order d_a = (l − χ₄(l)) l^{a−1} (sharp.md Prop. 2.1). Write t_j = r_j/conj r_j = g^{L_j}. The constraints are:
  * parity: L_j ≡ [(p_j/l) = −1] (mod 2). This is exactly the norm condition of Def. 1.1(d), by Prop. 2.1(4).
  * block sums: Σ_{j∈J} L_j ≡ Σ_{j∈J} L_j^old (mod d_a) for every block J.
  * digits: for each prime r | d_a with r^e ‖ d_a and r^e ≥ U, the residue L_j mod D = r^{k_r} lies in a cell of (X.1) determined by the class of j.

  The residues r_j are then any norm-p_j units with these ratios: a genuine factor times a norm-one u with u² equal to the required quotient. The data at different l are chosen **independently**. At each l the levels are made compatible by König's lemma on the finitely branching tree of admissible finite data. This step is correct, because reduction preserves the parities, the block sums and the gcd property. The block products ∏_{j∈J} t_j remain reductions of the **global** IX.5 block owners. Only the individual labels are local.
* **Parent primes.** At p_l the data are kept exactly as in IX.3/IX.5: reductions of the global witnesses with the singular factor γ_l deleted. Note that (1.2) runs only over q ∤ N, so parent primes play no role in e-all. Only the same-level clause of e-char involves them, and it is unchanged.
* **P3-char at every X (Proved).** At exterior q every zero-sum mark character is ∏_x ρ_x^{c_x} = i^{−c·k} ∏_J (∏_{j∈J} t_j)^{Σ_x c_x a_xJ}. This depends only on the block products, which are unchanged; δ and k are unchanged too. So (1.1) holds verbatim as in IX.3/IX.5, with the same uniform-in-X bound. Parent-prime data are unchanged. (a), (b) and (d) are immediate.
* **Fine angles.** Independently for each block, φ|_J is drawn from normalised Haar measure on the coset {Σ_{j∈J} φ_j = Φ_J}. This coset is a connected (g − 1)-torus. Real lifts are taken with exact block sums.
* **Weight.** W = log N' = Σ_j 2 log p_j = 6g(q − 1) log q (1 + o(1)) and m log m = (ρ/2) q log q (1 + o(1)), so W/(m log m) → 12g/ρ. Nothing depending on the residue data or on any lift height enters W. **Confirmed (Proved):** the weight charged is the physical parent weight.

## 3. Single labels and the general index argument (task item 2)

Let Λ = Σ v_j L_j, so that ∏ t_j^{v_j} = g^Λ.

1. **Collision ↔ index.** i ∈ T_{l^a} has order 4. If (g^Λ)^n = i^s mod l^a, then 4nΛ ≡ 0 mod d_a, so d_a/gcd(d_a, Λ) divides 4n and d_a | 4n·gcd(Λ, d_a). Prime powers a ≥ 2 are included, since d_a contains l^{a−1}. The 2-part and the unit are handled together by the factor 4. **Proved;** I also checked it with actual residues (§6.2).
2. **(X.2): gcd(Λ, d_a) | G for every tested word.** For each r: if e < k_r, the r-part is at most r^e ≤ r^{k_r−1} < U. If e ≥ k_r, the digit lemma gives D = r^{k_r} ∤ Λ, so the r-part is at most r^{k_r−1} < U. Every prime-power part of the gcd is < U, so the gcd divides G = lcm(1, …, ⌈U⌉). The property is intrinsic (it says the order of ∏ t_j^{v_j} in T_{l^a} is at least d_a/G) and passes to reductions. **Proved.**
3. **Collision product.** l^a ≤ (l/(l − χ)) d_a ≤ (5/4)·4nG ≤ 6nG. So m_(nv,s) ≤ lcm(1, …, ⌊6nG⌋) = e^{O(nG)}, uniformly in X and maximised over units separately. **Proved.**
4. **Angles.** Every tested v has 0 < ‖v‖₁ ≤ H < g, so it is nonconstant on some block J. A coefficient c ≠ 0 constant on J would cost mass ≥ g. On the Haar coset of J, n⟨v, φ⟩ is uniform mod 2π conditional on all other blocks. The forbidden set for the four units therefore has measure O(e^{O(nG)} P^{−n‖v‖₁/2}), using w(nv) ≥ n‖v‖₁ log P. The union bound over the ≤ (2s)^h words of mass h ≤ H and all n ≥ 1 is
   Σ_h (2gq)^h Σ_n e^{O(nG)} q^{−3nh/2} = O_H(q^{−1/2+o(1)}) = o(1),
   because G = e^{(1+o(1))U} = (log q)^{o(1)}. The mixed ordinary gaps (G) cost Σ_h (2gq)^h P^{−h/2} = o(1). Block-constant words are coarse and were already handled. So a set of positive measure succeeds, and removing a null set gives nongenuine angles. **Proved.**
5. **Singletons specifically.** For v = n e_j the class sum is n ≠ 0, so items 2–4 apply. With b = 1 (one row c = (1)) all words with Σ v_j ≠ 0 are tested, so the single-label claim needs no growing partition. The "local order checks" in the referee's question are exactly item 2 read as order(t_j^{n}) ≥ d_a/(G·gcd(n, d_a)) at every l, a. **Proved.**
6. **Compatibility with the residues forced by P3-char.** The mark ratios Ω_i/Ω_{i'} (Ω_i the IX.5 owner ratios) are products of block products, which are preserved. Parities are the physical ones. **Proved.**

## 4. The digit-sum covering (task item 3)

**Lemma (referee statement of (X.1)).** Let H ≥ 1, B = 64H², g ≥ B − 1 (here g = 2B), and rows c_t ∈ {0,1}^b (t = 1..T). Let D be a positive integer with D ≥ 100·g·B^T. For a class κ, call an integer L *κ-admissible* if L mod D lies in some cell [Dz/B^T, D(z+1)/B^T] with z = Σ_t z_t B^{T−t}, z_t ∈ {a_t, a_t + 1}, a_t = 16H·c_t[κ]. When D is even, additionally require a prescribed parity.

* **(A) Attainment.** For any g slots with any classes, any parities, and any S whose parity agrees with the sum of the parities when D is even, there are admissible L_j with Σ L_j ≡ S (mod D). The same holds mod r^e when D = r^{k} | r^e: correct one slot by a multiple of D.
* **(B) Digit intervals.** For admissible L of class κ and each t, B^{t−1}L/D mod 1 ∈ [a_t/B, (a_t + 2)/B].
* **(C) Covering.** If ‖v‖₁ ≤ H and some row has u = Σ_j v_j c_t[cl(j)] ≠ 0, then D ∤ Σ v_j L_j.

**Line-by-line check.**

* (A) In digit t the g slots contribute a digit sum ranging over the full interval [Σ a_t(j), Σ a_t(j) + g]. That is g + 1 ≥ B consecutive integers, whatever the classes are. Choosing digit sums from the least significant digit upward, with carries, gives Σ z ≡ τ (mod B^T) for every τ. With τ = ⌊S B^T/D⌋ − ⌊g/2⌋, the translate S\* of S satisfies S\* − DZ/B^T ∈ [⌊g/2⌋, ⌊g/2⌋ + 1)·D/B^T. This holds because DZ/B^T changes by D when Z changes by B^T. Note that DZ/B^T need not be an integer, since D and B^T are coprime for odd r; the argument is in exact rationals and is unaffected. The attainable sums form all integers (all integers of the right parity, when D is even) between Σ⌈left⌉ and Σ⌊right⌋. The loss from rounding is < 2 per endpoint, so < 4g in total. Since D/B^T ≥ 100g, both margins ⌊g/2⌋·100g − 4g and (g − ⌊g/2⌋ − 1)·100g − 4g are positive. **Correct.**
* (B) B^{t−1}·z/B^T = Σ_u z_u B^{t−1−u}. The digits with u < t are integers. The digit u = t gives z_t/B. The tail Σ_{u>t} z_u B^{t−1−u} plus the cell length B^{t−1−T} is ≤ 1/B, *because every z_u ≤ 16H + 1 ≤ B − 1* (this is where B = 64H² ≥ 16H + 2 is used; Astra's sentence "the entire tail of the cell occupies at most one extra digit unit" is this computation). So the representative lies in [z_t/B, (z_t + 1)/B] ⊂ [a_t/B, (a_t + 2)/B]. The interval is closed, and the right endpoint can be attained in the limit. **Correct.**
* (C) Σ_j v_j B^{t−1}L_j/D ≡ u/(4H) + E (mod 1) with |E| ≤ ‖v‖₁·2/B ≤ 2H/B = 1/(32H). Since 1 ≤ |u| ≤ ‖w(v)‖₁ ≤ H, the number u/(4H) is at distance ≥ 1/(4H) from Z, so the sum is nonzero mod 1. Hence D ∤ B^{t−1}Λ and so D ∤ Λ. **Correct.** The requirement is only B > 8H². The choice 64H² is not optimal, and g could be reduced to about B as well. This would lower the constant 1536H², but not its order H².
* **What it is used for.** At every exterior l and every prime r | d_a with r^e ≥ U, each label's log is confined mod r^{k_r} to class-dependent cells, while the block sums remain free. The point is that block sums (P3-char) and the word covering (e-all) are compatible: the covering uses only the residues mod D in cells of length ≥ 100g, and g slots with g ≥ B can absorb any prescribed sum.

**Observation (Computed).** In my greedy realisation of (A), most slots sit at their left endpoint. Slots of the same class and the same cell then get **identical** logs mod D. In the test, 67,706 untested words (w(v) = 0) had D | Λ. So same-class differences can collide at every digit component. This is consistent with X.2, which makes no claim about them.

## 5. Freedoms that genuine circles do not have (task item 4)

1. **Independent local residues at each exterior prime.** At each exterior prime the residues are chosen with no global object. For a single label, IX.7 (contrapositive) shows there is no cofinite rational lift. Its hypotheses, (G) and the unit −1 singleton rays at every cutoff, are satisfied, and the angles are nongenuine. **Proved.** Only the block products are global.
2. **Collision-poor data, i.e. failure of the product formula.** For a genuine circle and a word v, the four integers D_s = Y_v, X_v ± Y_v, X_v have |D_s| = ν_s^{−1} e^{w(v)/2}|sin(⟨v,φ⟩ − sπ/4)|. As X → ∞ the collision product m_(v,s) equals the part of |D_s| prime to 2N. In X.2, by contrast, m_(nv,s) ≤ e^{O(nG)} = p^{o(n)} for every tested word, at every cutoff. The local orders are near-maximal simultaneously at all l, which no global element can achieve. This is the decisive freedom.
3. **Angles decoupled from residues.** The fine angles are Haar on block fibres and are tied to the residues only through the inequalities.
4. **Weight.** W is the physical log N'. Finite-cutoff rational lifts r = πγ/conj γ exist at each X (README_IX §1.8, finite-cutoff lift proposition, height e^{O(X)}), but no lift height enters W. **Confirmed.** The quantitative lower bound (X.3) for those lift heights is in block rays2 and was not checked here.

## 6. Computations (Computed; referee's own code)

### 6.1 `digit_ref.py`

Four parameter sets were run: (H, b, rows) = (2, 3, one-hot), (2, 4, random T = 14), (3, 3, one-hot), (3, 3, random T = 18). In each set, detection was verified exhaustively on all 24, 40, 62 and 62 nonzero class vectors respectively. For r ∈ {2, 3, 5, 7, 13}, e ∈ {k_r, k_r + 2}, with 4 trials × 2 blocks each, slots were given random or constant classes and random parities, and random targets were used. For each block the run checked:

* the exact location claim S\* − DZ/B^T ∈ [⌊g/2⌋, ⌊g/2⌋ + 1)D/B^T and that S\* is attainable;
* the parities;
* the lift to r^e;
* (B) for every slot and every digit, in exact integer form a_t D ≤ B·(B^{t−1}L mod D) ≤ (a_t + 2)D;
* (C) for every signed word of mass ≤ H on a label subset spanning both blocks.

Output:
```
blocks built 320 digit-interval checks 2631680 tested-word checks 6985602 untested words with D | Lambda 67706
ALL DIGIT CHECKS PASSED
```
Astra's `x_digit_check.py`, run unmodified (not imported), printed `digit blocks 200 slots per block 1152 signed word tests 32438 individual digit intervals 921600`, matching her report. Her script tests one block, one-hot rows on 4 colours, no lift to r^e and no CRT. Mine adds two blocks, random rows, the lift, and (in 6.2) the CRT.

### 6.2 `gauss_ref.py`

The check uses actual Gaussian residues for l = 3, 5, 7, 13, at level a = k_l + 1, so that the l-component carries digits (modulus about 10¹³; U = 8.6·10¹¹ for H = 2, b = 3). There are 1,024 physical primes in 2 blocks of 512. The run:

1. built a generator t_0 = z_0/conj z_0 of T_{l^a} and checked its order;
2. ran the digit construction at the l-component and block sums plus parity at the other components, combined by CRT;
3. formed r_j = sqrt(p_j/N(z_0)^{L_j}) z_0^{L_j} and checked Norm r_j ≡ p_j and r_j/conj r_j = t_0^{L_j};
4. checked that the opposite parity gives a nonresidue, so that no norm-p_j factor exists for it;
5. checked that the block products of the new t_j equal those of the old data.

Then, for every word of mass ≤ 2 on 34 labels, every n ≤ 12 and every unit, it computed the exact collision level a' by Gaussian arithmetic and checked d_{a'} | 4n·gcd(Λ, d_{a'}). For tested words it also checked v_l(gcd(Λ, d_a)) ≤ k_l − 1 and l^{a'} ≤ 6n·gcd(Λ, d_a). Results: 2,024 / 2,026 / 2,008 / 2,010 tested words, 0 failures, with maximal tested collision levels 8, 6, 4, 4. `ALL GAUSSIAN LOCAL CHECKS PASSED`.

### 6.3 `fraction_ref.py`

This counts, exactly and for balanced classes, the fraction of all words of mass h with zero class-sum vector:

| s | b | h = 1 | h = 2 | h = 3 |
|---|---|---|---|---|
| 24 | 2 | 0 | 0.229 | 0 |
| 24 | 4 | 0 | 0.104 | 0 |
| 24 | 6 | 0 | 0.0625 | 0 |
| 48 | 8 | 0 | 0.0521 | 0 |

Odd masses are always fully covered, because ‖w‖₁ ≡ h (mod 2). At h = 2 the uncovered words are exactly the same-class differences, a fraction (s/b − 1)/(s − 1) of the differences. This matches her 1/b + o(1).

## 7. Quantifiers (task item 5)

The order of choice is:

1. C, H and ρ are fixed first. Then g = 128H², B = 64H², d = 3, and the limiting weight constant 12g/ρ are determined.
2. A schedule b(q) → ∞ with U = o(log log q) is chosen. It does not depend on the prime set.
3. For each large q, in this order: the Paley marks; the physical primes and their assignment; the coarse data of IX.3/IX.5; the partition into b(q) classes and the rows; the local data at all exterior l (independently, each a compatible family); and finally the fine angles, from a set of positive measure.

No exceptional set depends on X: all bounds are uniform in X. The threshold q_0 depends on C, H, ρ and b(·). It is not effective in any practical sense, because e^{O(G)} = q^{o(1)} requires log q ≫ e^{U}, and U ≥ 100·512·256² ≈ 3.4·10⁹ already at H = 2, b = 2. "Passes" means inequality (1.2) for nv at all four units, all n ≠ 0 and all X. No standing caution of the record is violated:

* the record's Subspace-theorem caution is not touched, because no Subspace theorem is used;
* no Gaussian gcd or denominator is discarded, because W is physical and no lift enters;
* the common phase θ_0 is kept, with the block sums fixed exactly;
* there is no template-dependent threshold used as if it were uniform.

## 8. Verdict table

| claim | verdict |
|---|---|
| Statement of X.2 (first and second bullets; third bullet as corrected in §1) | CORRECT WITH CORRECTIONS (meaning of "fraction", the H² constant, data-independence of the word set) |
| Coarse profile and preservation of P3-char at every X, including parent-prime levels | CORRECT (Proved) |
| Class encoding (random detecting rows) | CORRECT; one-hot rows also suffice |
| Digit-sum covering (X.1): attainment, digit intervals, covering | CORRECT (Proved and Computed) |
| Index bound (X.2) and collision bound l^a ≤ 6nG at every level/cutoff | CORRECT (Proved and Computed) |
| Every power of every single label passes e-all at every unit and every cutoff | CORRECT (Proved; holds already with b = 1) |
| Angular completion, union bound, nongenuine angles, no cofinite rational lift via IX.7 | CORRECT (Proved) |
| W = O_{C,H}(m log m), m ≍_H s | CORRECT, constant 1536H²/ρ + o(1), m/s → ρ/(256H²) |
| Fraction 1 − o(1) of words of each mass; ≥ (1 − 1/b − o(1)) C(s,2) differences | CORRECT; sharper: every odd-mass word is covered |
| Deliberate failure, (X.3) auxiliary heights, limit of order argument | NOT CHECKED (block rays2) |
