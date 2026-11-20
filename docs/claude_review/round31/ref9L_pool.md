# Referee report, round 31, Part IX, block "pool"

Scope: Astra's Part IX, §IX.0, from "PROVED (main theorem)" to the end of "PROVED (primes below a polynomial cutoff)"
(`astra_lemma_IX/answer.md` lines 1939-1971). That covers the main theorem (IX.0), the integer triangle, the Gaussian
triangle, the probability argument, the corollaries and the polynomial-cutoff version. The two bounded-exponent
paragraphs (lines 1973-1988) and the same-norm P2 comparison (line 1990 on) belong to other blocks. I read the twin
version of the bounded-exponent paragraph only because it uses the same union bound (§3.6). Scripts are my own and
are in `round31/ref9L_pool/`: `twin_sim.py`, `brute_test.py`, `circle_sim.py`, `gauss_triangle_check.py` and
`union_bound.py`. I did not import Astra's `ix_checks.py`.

## 0. Verdict

**The block is correct.** I checked every step by hand and recomputed every numerical inequality. Nothing is false.
Nothing needs a correction of substance. I make four precision remarks: the implied constant is absolute, the
dependence on C enters only through s_0(C), coverage on circles is as stated in §4, and the novelty is narrower than
the summary suggests.

| Item | Label | Verdict |
|---|---|---|
| (IX.0), twin, all intervals and all residue choices at once | Proved | CORRECT (verified statement in §1) |
| Integer triangle (IX.0a) and the deterministic heavy block | Proved | CORRECT |
| Gaussian triangle, z_j - z_i = 2iηH_ij t_ij with t_ij in Z \ {0} | Proved | CORRECT (also checked exactly on 300 random triples) |
| Probability and union bound (4^s and 7^s states, (2Y)^3, binom(K-s+k,k)^{-1}) | Proved | CORRECT (exponent recomputed) |
| Corollaries P = A s log s (O(log s)) and P = s^{1+ε} (O(1)) | Proved | CORRECT (needs A ≥ e^40, or 2e^40 for split pools) |
| Primes below X ≤ s^A, no lower bound on the primes | Proved | CORRECT |
| Simulation, s = 6..10 | Computed | B_max = 3 or 4 and M_max ≤ 3, far below every scale; does not test asymptotics |
| Novelty | — | the twin statement is new. The circle statement and the exposure lemma were already in the parallel thread's unrefereed notes. The triangle (three marks, deterministic heavy block) is a genuine simplification of the quartet proof. |

**Significance.** The block proves a rate o(ω(q)) on a random family of parents: B, M = O(log s) = O(log log q) at
log q ~ s log s. Before this, the twin had no such statement at all, and the circle version existed only in unrefereed
notes. The family has density zero, and the theorem says nothing about any individual modulus or norm. Like every
random-parent theorem in the record, it leaves the worst-case growth rate for general R unchanged.

## 1. Statement as verified (Proved)

**Theorem (IX.0, verified form).** Fix C > 0. There is s_0(C) such that the following holds for all s ≥ s_0(C) and
all P, K. Let 𝒫 be any set of K primes in [P, 2P] with K ≥ e^40 s. Put J = log(K/s) (so J ≥ 40) and
b = max(3, 1 + ⌈100 log(2P)/J⌉). Let S be a uniformly random s-element subset of 𝒫 and q = ∏_{p∈S} p. Outside an
event of probability at most exp(-sJ/20), the following hold.

1. *(twin)* For every real interval I of length C√q, and for every choice of two (possibly equal, possibly non-unit)
   residues r_p, r'_p mod p at each p | q, at most b - 1 integers y ∈ I satisfy y mod p ∈ {r_p, r'_p} for every p | q.
2. *(circle)* If 𝒫 consists of primes ≡ 1 mod 4, outside a second event of the same probability bound, every arc of
   length C q^{1/4} on x² + y² = q contains at most 4(b - 1) lattice points.

Order of quantifiers: C is fixed first. s_0 depends only on C, and **not** on P, K or the pool. The residues, the
interval and the arc are chosen after q. The bounds b - 1 and 4(b - 1) are **absolute**: the constant in Astra's
"≪_C" may be dropped, since C enters only through s_0(C). One explicit choice is
s_0(C) = max{(3.2 + 9 log C_0)/2.24, log C/3.7, 1}, with C_0 = max(2, C, C²); the "large s" in the Gaussian
heavy-block step needs only a similar bound. Two remarks:

* *Uniformity in C (Proved, a free extension).* The events in the union bound grow with C, and the C-dependence
  costs 9 log C in the exponent. With s_0 absolute, the conclusion therefore holds simultaneously for **all**
  C ≤ exp(sJ/400). The exponent becomes ≤ -0.084 sJ (`union_bound.py`). The heavy block survives because
  log C/2 ≤ sJ/800 ≤ s(log P/20 - (log 2)/5) when log P ≥ 40. This covers intervals up to q^{1/2}(K/s)^{s/400}.
  At those lengths the larger sieve already gives the same order, so the extension is cosmetic.
* *Simultaneity of the two assertions:* their union has probability ≤ 2 exp(-sJ/20), as stated.

## 2. The integer triangle (Proved, every step checked)

Setup. Let y_1 < ... < y_b be distinct admissible integers in I. Put h_ij = ∏{p : y_i ≡ y_j mod p}. Then h_ij | y_j - y_i,
so h_ij ≤ |y_j - y_i| ≤ C√q, and t_ij := (y_j - y_i)/h_ij is a nonzero integer. Equal residues r_p = r'_p and
non-unit residues need no special treatment, because only congruences between marks are used.

* **Bit separation.** A prime with k of the b marks in one class separates k(b - k) ≤ b²/4 pairs. Hence
  Σ_{i<j} log(q/h_ij) ≤ W b²/4, with W = log q, and the mean over the binom(b,2) pairs is ≤ W/2 + W/(2(b-1)).
  So mean log h_ij ≥ W/2 - W/(2(b-1)), and mean log|t_ij| ≤ log C + W/(2(b-1)). CORRECT.
* **Good triple.** Each pair lies in exactly b - 2 triples, and (b-2)·binom(b,2)/binom(b,3) = 3. So the mean over
  triples of log|t_12 t_23 t_13| is 3 times the pair mean, and some triple has
  log|t_12 t_23 t_13| ≤ 3 log C + 3W/(2(b-1)). CORRECT.
* **(IX.0a).** Let G be the product of the primes at which the three marks are congruent. At any other p the three
  marks occupy exactly two classes. This is the place where "at most two residues" is used. So exactly one of the
  three pairs agrees at p. Hence h_12 = G A_12, h_23 = G A_23 and h_13 = G A_13, with q = G A_12 A_23 A_13 and the
  four factors pairwise coprime. Dividing (y_2 - y_1) + (y_3 - y_2) = y_3 - y_1 by G gives
  t_12 A_12 + t_23 A_23 = t_13 A_13. CORRECT. The identity and the factorisation were also asserted exactly by
  `twin_sim.py` on every admissible triple it found.
* **Heavy block on every triple.** At the primes separating marks 1 and 2, mark 3 agrees with mark 1 (block A_13) or
  with mark 2 (block A_23). Their product q/h_12 has log ≥ W/2 - log C. So one of A_13, A_23 has
  log ≥ W/4 - (log C)/2, and that block contains at least (W/4 - (log C)/2)/log(2P) primes. For "≥ s/5" one needs
  W/4 - (log C)/2 ≥ (s/5) log(2P). Using W ≥ s log P, it suffices that s((log P)/20 - (log 2)/5) ≥ (log C)/2. This
  holds when P > 16 and s ≥ s_0(C). Here P > 16 is automatic, since K ≤ #primes in [P, 2P] < P gives
  P > e^40 s. CORRECT. The bound is deterministic and holds on every triple, so no independence between the block
  and the cofactors is needed, as Astra says.

## 3. The Gaussian triangle (Proved)

* Since q is a squarefree product of split primes, unique factorisation in Z[i] gives z = η ∏_p g_p with g_p ∈ {π_p,
  π̄_p} and η a unit. The representation is unique once the generators π_p are fixed (for example by a fixed rule on
  a + bi). Keeping one unit class loses a factor of at most 4, hence the bound 4(b - 1). CORRECT.
* For i, j in the same class, let X = ∏_{p separating} g_ip. Then z_i = ηH_ij X and z_j = ηH_ij X̄, so
  z_j - z_i = ηH_ij(X̄ - X) = -2iηH_ij Im X. So t_ij = -Im X is a rational integer, nonzero because z_i ≠ z_j.
  Moreover |z_j - z_i| = 2√h_ij |t_ij|, with h_ij = |H_ij|² = ∏ p. CORRECT. Exact check
  (`gauss_triangle_check.py`): the identity, the integrality of t_ij and the triangle relation
  t_12 A_12 + t_23 A_23 = t_13 A_13 in Z[i] hold on 300 random triples, with 10 primes from [100, 2000] and random
  units.
* Chord ≤ arc ≤ C q^{1/4} gives 2√h_ij ≤ C q^{1/4}, i.e. log h_ij ≤ W/2 + 2 log(C/2). The same bit separation gives
  mean log|t_ij| ≤ log(C/2) + W/(4(b-1)), which is half the twin's coefficient cost. Replacing C by
  C_0 = max(2, C, C²) is a harmless overestimate. The separating mass is ≥ W/2 - 2 log(C/2). CORRECT.
* Dividing by 2iη and by the product of the Gaussian factors common to all three marks gives (IX.0a) in Z[i]. Each
  non-constant rational prime appears in exactly one A_ij, as π_p or π̄_p. The arc centre is never used, only chord
  lengths. CORRECT.

## 4. The probability argument (Proved, the core, checked as a rigorous statement)

*Probability space.* S is a uniform s-subset of 𝒫. Give it an independent uniform ordering. Then (p_1, ..., p_s) is a
uniform injective sequence, and "S is bad" has the same probability in both models. (The record's earlier slip in
Part III, a sorted ordering, does not occur here.)

*Events.* An event E(λ, F, t) is indexed by three things:
- a labelling λ of the slots {1, ..., s}: 4 labels for integers (G, 12, 23, 13), and 7 labels for Gaussian integers
  (G, or one of the 3 blocks together with one of 2 orientations relative to the fixed π_p). **7 = 1 + 3·2 is
  right.** G slots need no orientation, because G cancels;
- a forced block F ∈ {12, 23, 13} with k := |λ^{-1}(F)| ≥ s/5 (at most 3 choices);
- a signed triple t ∈ ([-Y, Y] ∩ Z \ {0})³, with Y = C_0³ exp(3 s log(2P)/(2(b-1))).

E holds when t_12 A_12 + t_23 A_23 = t_13 A_13, where A_X is the product of p_slot (or of its labelled Gaussian
factor) over the slots labelled X.

*Coverage.* Suppose some interval, some residue choice and b admissible marks (or b points in one unit class)
exist. Take the good triple of §2. Its labelling is read off the triple's classes in the sample's own order. Its
heavy block has k ≥ s/5. Its cofactors satisfy |t_ij| ≤ |t_12 t_23 t_13| ≤ Y **individually**, because every
|t| ≥ 1. So the event E(λ, F, t) holds. **Residues and the interval are not enumerated, and need not be:** E is a
statement about the sampled primes alone, and the triple, the labelling and t are existential witnesses. Choosing
the triple after seeing q is therefore harmless. The block structure depends on q, but every possible structure is
in the union. CORRECT.

*Conditional probability.* λ and F are fixed before sampling. Condition on the values in the s - k slots outside F.
The F-slots are then a uniform injective k-tuple from the K - s + k unused primes. Given t, (IX.0a) fixes
A_F = (±t A ± t' A')/t'' as a single number (in Z[i] in the Gaussian case, where the norm |A_F|² is the rational
product). Unique factorisation in Z gives at most one unordered k-set of primes with that product. Hence
P(E | exposed) ≤ k!/(K-s+k)_k = binom(K-s+k, k)^{-1}. Next, (K-s+k)/k ≥ K/s is equivalent to (s-k)(K-s) ≥ 0, and
binom(n,k) ≥ (n/k)^k. So the probability is at most (s/K)^k ≤ exp(-sJ/5). CORRECT.

*Union bound (recomputed, `union_bound.py`).* The log of the union bound is at most
log 3 + s log 4 (resp. 7) + 3 log 2 + 9 log C_0 + 9 s log(2P)/(2(b-1)) - sJ/5. With b - 1 ≥ 100 log(2P)/J the fourth
term is ≤ 0.045 sJ, and with J ≥ 40, s log 7 ≤ 0.0487 sJ. The total is ≤ -0.106 sJ + O(log C_0) in the Gaussian case
and ≤ -0.120 sJ + O(log C_0) for integers. This is ≤ -sJ/20 once s ≥ s_0(C). CORRECT, with a margin of 0.056 sJ.

*Hidden dependencies:* none. Only λ and F must be fixed before the exposure, and they are, as enumeration indices.

*Coverage on circles (scope, not an error).* The theorem concerns points of norm exactly q. It also covers clusters
on n = q|Γ|² all of whose points are divisible by a fixed Gaussian integer Γ: dividing by Γ shortens chords by
|Γ| ≥ |Γ|^{1/2}. It does not cover norms whose Gaussian-gcd-normalised norm is not a sampled q. That excludes the
cap-two norms of P2: the squarefree theorem does not apply to them, and the same-norm comparison relies on the
bounded-exponent paragraph, which is refereed in another block.

## 5. Corollaries and the polynomial cutoff (Proved)

* **P = A s log s.** The number of primes in [P, 2P] is ~ P/log P ~ A s, so K/s → A and A ≥ e^40 (1 + o(1)) is
  needed. For split primes K ~ A s/2, so A ≥ 2e^40. Then J ≈ log A, and b - 1 ≈ 100 log(2P)/log A = O(log s). Also
  log q = s(log s + log log s + O(1)) ~ s log s, so log q/log log q ~ s. CORRECT: the bound is O(log s) where the
  larger sieve gives order s, outside a fraction e^{-Ω(s)}.
* **P = s^{1+ε}.** K ~ s^{1+ε}/((1+ε) log s) (halved for split primes), J = ε log s - log log s - O(1) → ∞, and
  log(2P)/J → (1+ε)/ε. So the bound is b - 1 ≤ 2 + 100(1+ε)/ε + o(1): O_ε(1), for s ≥ s_0(C, ε). CORRECT. On circles
  this is the regime D = 1 + ε of the record's refereed random-parent theorem (finite M_D for every D > 1), with a
  different and explicit constant.
* **Primes below X ≤ s^A, no lower bound.** Since the primes are distinct, q ≥ ∏_{i≤s} p_i ≥ s!, so
  W ≥ log s! = (1 - o(1)) s log s, and W ≤ s log X ≤ A s log s. The heavy block has
  k ≥ (W/4 - (log C)/2)/log X ≥ s/(5A) for s ≥ s_0(A, C). With b - 1 ≥ 100 A log X/J and J ≥ 40A, the exponent per
  sJ/A is ≤ log 7/40 + 9/200 - 1/5 = -0.106. CORRECT: the bound is O_A(1 + log X/J), with exceptional fraction
  exp(-Ω_A(sJ)). Small primes cause no problem in the twin: at p = 2 the two classes impose nothing, but the
  argument never uses that a restriction is informative. "Cutoff of order s log s" means X = A′ s log s with
  log A′ ≥ 80 (A = 2): the bound is O(log s), again with an astronomically large constant.
* *(Outside the block, read for the shared union bound.)* In the twin with bounded exponents 1 ≤ e_p ≤ E, the pair
  divisibilities use the full p^{e_p}. The separated pairs then have t_ij divisible by p^{v_p(r - r')}, which is
  harmless. There are 1 + 3E labels, the block has ≥ s/(5E) slots, and the exponent per sJ/E is
  1/40 + 9/200 - 1/5 = -0.13 < 0. The twin paragraph is correct. I did not referee the Gaussian bounded-exponent
  extraction.

## 6. Simulation (Computed; it does not test asymptotics)

*Method (exact integers).* A finite set occupies ≤ 2 classes mod p if and only if every triple of it has a pair
congruent mod p. For squarefree q, a set is therefore admissible for some residue choice exactly when
q | (b-a)(c-b)(c-a) for all of its triples a < b < c. After translating the minimum to 0, `twin_sim.py` enumerates
all admissible triples (0, a, b) with 0 < a < b ≤ L = ⌊C√q⌋. For each of the 3^s assignments of primes to the
factors a, b, b - a, it enumerates the index-q lattice {q_1 | a, q_2 | b, q_3 | b - a} in the triangle, using a
Gauss-reduced basis. It then runs an exact clique search with all inner triples checked. This gives B_max(q, C)
exactly: the maximum over all residue choices and all intervals. The code was validated against exhaustive brute
force on 120 cases (`brute_test.py`; s = 3, primes < 40, C = 1, 2, 3), with identical triple sets and maxima.
`circle_sim.py` lists all 4·2^s points of norm q exactly, finds the maximal arc count by an angular sliding window,
and checks the pair identity exactly.

*Twin results* (`twin_runs.txt`, `twin_runs2.txt`):

| s | pool | C | samples | B_max | log q/log log q |
|---|---|---|---|---|---|
| 6 | primes in [100,200] | 2 | 50 | 3 (all) | 8.8 |
| 6 | [100,200] | 4 | 50 | 4 (all) | 8.8 |
| 6 | [100,200] | 8 | 30 | 4 (all) | 8.8 |
| 6 | [1000,2000] | 4 | 30 | 4 (all) | 11.6 |
| 8 | [100,200] | 2 / 4 | 30 / 30 | 3 / 4 (all) | 10.9 |
| 8 | [300,600] | 8 | 10 | 4 (all) | 12.6 |
| 10 | [100,200] | 2 | 10 | 3 (all) | 12.8 |
| 10 | [1000,2000] | 2 | 10 | 3 (all) | 17.0 |

The values 3 and 4 are the generic floor, not growth. Every 4-set found was a two-dimensional cube
{0, q_1, q_2, q_1 + q_2} with q_1 q_2 = q (checked exactly for the extremal sets). Such a cube exists whenever q has
a divisor with q_1 + q/q_1 ≤ C√q. A cube of dimension 3 would need two of its generators above q^{2/3}, which is
consistent with IX.1. The triples come from Minkowski: the triangle has area C²q/2 against covolume q.

*Circle results* (`circle_runs.txt`). For split pools in [100, 200], [200, 400] and [1000, 2000], s = 6..14 and
C = 1..8, M_max ∈ {1, 2, 3}, with a single 3 in 185 samples. All same-unit pairs on the extremal arcs satisfied
z_j - z_i = 2iηH t with t ∈ Z \ {0}.

*Reading.* At these sizes the hypotheses K/s ≥ e^40 fail by a factor of order 10^16, and Astra's bound b - 1 would be
several hundred. The data are consistent with the theorem but test nothing asymptotic. **Heuristic:** a first-moment
count with independent residues predicts O(1) genuine marks already at P = A s log s. Concretely, a b-set occupies
≤ 2 classes at each p with probability ≈ 2^{b-1}/p^{b-2}, which gives ≈ C^{b-1} exp(s[(b-1) log 2 - ((b-3)/2) log P])
expected b-sets. So the factor log s in (IX.0) is very probably an artefact of the 4^s / 7^s entropy against the
cofactor cost, not the truth.

## 7. Novelty: what the record had

* **Refereed record** (`expert_note.md` lines 896-902; `round31/README_III.md` §1, §2 and the wording in §3). It has
  random parents on circles only, with split primes from [r^D, 2r^D] and **fixed D > 1**: M_D = 1 for D > 4, ≤ 2 for
  8/3 < D ≤ 4, O(1) for D > 2 (pair first moment), and finite for every D > 1 (quartet exposure). It also has the
  bounded-exponent extension (1 + 3E states), all on the random support. It has **no rate statement** at
  P ≍ s log s, and **no random-parent statement for the twin**.
* **Parallel thread, unrefereed notes** (`round31/astra_free_III/iii_graph_notes.md`). Lines 131-135 contain the
  exposure lemma in exactly Astra's form: three disjoint slot sets, signed coefficients ≤ U, the bound
  4^r(2U)³/binom(K-r+s, s), and a random ordering. Line 155 states the circle estimate
  O(1 + log P/log(K/r)) for K ≥ e^30 r, giving O_C(log r) at P = C r log r. It is proved through quartets, with a
  Jensen bound on the expected matching size. The referee of Part III (`ref3_density.md`) read lines 141-153 only.
  This estimate never entered README_III or the expert note, so it is unrefereed. Astra cites both items.
* **What is new here.** (a) The twin statement (IX.0)(1): the first random-parent theorem for the two-class problem,
  uniform over all residue choices and intervals, with rate O(log s) against the larger sieve's order s. (b) A
  three-mark proof: the signed triangle replaces the quartet's Plücker relation, and the heavy block is
  **deterministic on every triple**, so the Jensen and matching-proportion step of the quartet proof disappears.
  This is a real simplification, and it is what makes the twin case work. (c) A first refereed proof of the
  circle rate at P ≍ s log s. As a statement it is not new to the thread. (d) Bounded exponents and "all profiles"
  were already in the record for circles. In the twin they are new only as instances of (a).
* Bourgain-Rudnick-type almost-all statements concern density-one sets of n. This family has density zero, as in
  Part III.

## 8. What this does not do

It gives no bound for any individual q, no worst-case rate, and no improvement of Theorem 1.1 for general R. It does
not apply to narrow-bin parents or to the cap-two norms of P2 (§4). It has no implication for Conjecture 9 beyond
the random family. Astra states these limits herself (summary item 11), correctly.
