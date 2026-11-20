# Referee report, block "recon": Part X.4 (arbitrary cuts at the critical length)

Material refereed: answer.md lines 2203-2230 (section X.4), with the two limiting examples of X.5 (lines 2243-2245). X.1 (q^2 g^2 | H, and a triple with max edge < sqrt(qg) is degenerate) and X.2 (line lemma for at least five vertices) are taken as given; another referee checks them. All scripts are mine, in `ref10L_recon/` next to this report. Her script was not imported.

## 0. Verdict

- **First statement (B_int <= B_free <= 2 max(4, B_int) when q has a prime p > C^2; no factor 2 at C = 1): Proved.** It is correct as stated. One note on wording: the CRT alignment to prescribed residues may swap the two labels at a prime. The cut (the partition) is kept, and the unordered pair of residues is the prescribed one. If B_int is defined by optimising over the classes, as she defines it, no alignment is needed.
- **Second statement (dyadic parent, deletions O_C(1 + s/B)): Proved.** I checked every step. I found no gap in the three places she flagged: the alignment on common anchors, the small bit classes and the complementary pairs. Two constants are loose but valid. The far-graph degree is at most 16 (she uses 32). The complementary-pair bound 1 + floor(C^2) holds. A full corrected statement with explicit constants is in section 3.
- **Limiting examples: Computed.** Both reproduce exactly. They are the four-point table at q = 17·23·29, and the five-mark table at q = 59·67·71·83, C = 6, with H = -1152 q^2 on the triple (1,4,5).
- **New observations (Proved / Computed).**
  (i) If every prime of q exceeds C^2, then every non-complementary signature pair has **at most one** admissible magnitude. The only freedom in a sign-free table sits on complementary pairs.
  (ii) At q = 59·67·71·83, C = 6, with the classes of her five-mark example, B_free = 8 and B_int = 5. The 8-vertex table lives on the even-weight subgroup of {0,1}^4. Every choice of its complementary magnitudes needs at least 3 deletions. So B_free > B_int happens at a dyadic modulus, and the deletions in X.4 cannot be removed.
- **Scope.** The result shows that the sign-free relaxation is **not** a relaxation at scale s. On a dyadic parent with s >= s_0(C), B_int >= B_free - O(1 + C^2 + s/B_free). It gives **no** upper bound on B_int or B_free. It neither proves nor refutes B = o(s).

## 1. Definitions used

q is squarefree. For each p | q a nonzero δ_p mod p is fixed. A *signature* is σ_i ∈ {0,1}^{primes}. For two signatures, h_ij is the product of the primes where they agree, and D_ij = log(q/h_ij).

A *sign-free table* on B vertices consists of integers d_ij = d_ji with:
- 0 < |d_ij| <= C sqrt q;
- gcd(d_ij, q) = h_ij;
- d_ij ≡ ±δ_p (mod p) at each separating p, the sign free per pair and per prime.

Signatures are automatically distinct, because q | d_ij would force |d_ij| >= q > C sqrt q when q > C^2.

B_free(q, C, δ) is the maximal B. B_int(q, C, δ) is the maximal number of distinct integers in an interval of length C sqrt q that occupy, at each p, at most two residues differing by ±δ_p. Up to translation the classes are {0, ε_p δ_p}, ε_p = ±1.

## 2. First statement: Proved

**Lower bound.** A genuine configuration y gives a table d_ij = y_j - y_i. At a separating p, y_j - y_i ≡ ±δ_p ≢ 0, so gcd(d_ij, q) = h_ij exactly.

**Upper bound.** Let p > C^2. Split the table by the bit at p. Every triple inside one side has p | g, so g > C^2 and max edge <= C sqrt q < sqrt(qg). By X.1 the triple is degenerate. By X.2 a side with at least 5 vertices embeds: there are integers x_i with |x_i - x_j| = |d_ij|. Their spread is at most max |d| <= C sqrt q.

The embedded side is a genuine configuration:
- At each prime p', a same-bit pair has p' | d_ij, so x_i ≡ x_j. Hence each bit class has one residue.
- A cross pair has x_i - x_j = ±d_ij ≡ ±δ_{p'} ≢ 0. So the two residues ρ_0, ρ_1 satisfy ρ_1 - ρ_0 ≡ ±δ_{p'}. This value is automatically the same for all cross pairs.

So the side counts towards B_int(q, C, δ) when B_int optimises over classes. If the prescribed classes are a fixed pair {r_0, r_0 + δ}:
- When ρ_1 - ρ_0 ≡ δ, translate by r_0 - ρ_0 at that prime.
- When ρ_1 - ρ_0 ≡ -δ, translate by r_0 + δ - ρ_0. Then ρ_0 goes to r_0 + δ and ρ_1 goes to r_0. The labels swap and the partition is kept.

Distinct primes allow one global translation t by CRT. Translation does not change the interval length. Sides with fewer than 5 vertices contribute at most 4. Hence B_free <= 2 max(4, B_int).

**At C = 1.** X.1 makes every triangle degenerate, and X.2 then embeds the whole table if B >= 5. So B_free <= max(4, B_int).

**Computed.** stmt1.py brute-forces B_int (all ε) and B_free (maximum clique of the Cayley graph on {0,1}^k of admissible patterns). The data:
- 10 prime sets, k = 3..5 primes from 3 to 29;
- C^2 ∈ {1, 2, 4, 9};
- up to 400 δ-vectors per (q, C).

There was no violation. In every case B_free <= max(4, B_int). The maximum ratio B_free / max(4, B_int) was 1. In these small cases B_free <= 4. The test is therefore weak, which is why section 4 adds a larger modulus.

## 3. Second statement: Proved, with explicit constants

**Theorem (corrected statement, quantifiers explicit).** Let C > 0 and C_0 = max(1, C). Let q be squarefree with s prime factors, all in [P, 2P], where P > max(16, C_0^8). Since s <= π(2P), this holds for all s >= s_0(C). Fix nonzero δ_p. Let (d_ij) be any sign-free table on B vertices. Then one can delete at most

N <= 64 N_u + 32 + 1 + floor(C^2) <= 320 s/B + 97 + floor(C^2) vertices,

where N_u < 5 s/B + 1 is the number of "unbalanced" primes. If B <= 2277, delete everything; N = B is then O(1). After the deletions there are integers y_i such that:
- they lie in an interval of length C sqrt q;
- y_i ≡ ρ_{p, σ_ip} (mod p) at every p | q, with ρ_{p,1} - ρ_{p,0} ≡ ±δ_p;
- |y_i - y_j| = |d_ij| for **every** surviving pair.

The constants are absolute apart from the floor(C^2) term and the threshold s_0(C). Nothing depends on ε, on A, or on the number of primes per column: the cuts are arbitrary. Squarefreeness and P > C^2 are used. Dyadicity is used only in the count of unbalanced primes. With unequal primes the same proof gives Σ_{unbalanced} log p <= 4W/B + 8 log C_0. The number 32 can be replaced by 16.

**Step-by-step check.** Put W = log q and η = 2 log C_0 / W.

**(a) Proved.** h_ij | d_ij ≠ 0 gives h_ij <= C sqrt q, so D_ij >= W/2 - log C_0. Set v_i = (σ'_ip sqrt(log p / W))_p, where σ' = ±1 is the bit. Then v_i·v_j = 1 - 2D_ij/W <= η.

**(b) Proved, with a better constant.** A far pair means D_ij > 5W/8, i.e. v_i·v_j < -1/4. Let S be the sum of the vectors of the k far neighbours of i. Then |S| >= -v_i·S > k/4, so

k^2/16 < |S|^2 <= k + η k(k-1).

Hence k(1 - 16η) < 16(1 - η).
- For η <= 1/256 this gives k <= 16.
- k <= 32 needs only η < 17/512.

Under the hypotheses, W >= s log P >= 8 s log C_0, so η <= 1/(4s). Both bounds hold for s >= 64. Her Δ = 32 is valid. Δ = 16 is valid as well.

**(c) Proved.** In a triple, a nonconstant prime separates exactly two of the three pairs, and a constant prime separates none. So D_12 + D_13 + D_23 = 2(W - log g). If all three edges are good (D <= 5W/8), then log g >= W/16. We need g > C^2, i.e. W > 32 log C_0. This holds. X.1 then gives degeneracy.

**(d) Proved.** Greedy colouring of the far graph uses Δ + 1 colours. The largest colour class is a set of pairwise good marks of size >= B/(Δ+1). It has size >= 2Δ + 5 once B >= (Δ+1)(2Δ+5) = 2277 (for Δ = 32). Every triple in it is all-good, hence degenerate. With size >= 5, X.2 embeds it.

**(e) Proved: alignment on common anchors.** Let i ∉ R, and let G_i be the set of anchors good to i. Then |G_i| >= Δ + 5. All edges of {i} ∪ G_i are good, because anchors are pairwise good. So all its triangles are degenerate, and by X.2 it embeds, say by z.

Two isometric embeddings of a finite set T ⊂ Z with |T| >= 2 differ by a unique isometry of Z. Proof: fix two distinct points a, b of T. Any c ∈ T is determined by |c - a| and |c - b|, because two distinct points of a line have no common equidistant pair other than c itself. Let φ be the isometry carrying z|G_i onto the R-coordinates of G_i, and put y_i = φ(z_i). Then y_i is the unique integer whose distance to each anchor a ∈ G_i is |d_ia|. No choice is involved.

**(f) Proved.** Let (i, j) be good, with i, j ∉ R. They have at least |R| - 2Δ >= 5 common good anchors A_ij. The set {i, j} ∪ A_ij is all-good, so it embeds. Its anchor part aligns with R's coordinates by the unique isometry, which then sends i and j to the unique points with the prescribed anchor distances. These are y_i and y_j. So |y_i - y_j| = |d_ij|. If i ∈ R, the claim is immediate from the definition of y_j.

**(g) Proved: small bit classes.** Let K be a bit class at p with |K| > 2Δ. If the good graph on K were disconnected, take a component X with |X| <= |K|/2 and any x ∈ X. All of K \ X (more than Δ vertices) is far from x. This contradicts the degree bound. Her counting is correct, and it needs only one member of X.

Good same-class edges have p | d_ij = ±(y_i - y_j), so K has a single residue. If both classes are large, a vertex of K_0 has more than Δ vertices in K_1, so it has a good edge into K_1. That edge gives ρ_1 - ρ_0 ≡ ±δ_p. Small classes are handled only through (h).

**(h) Proved: unbalanced primes.** A prime is unbalanced when its minority class has at most 2Δ marks. Then |Σ_i σ'_ip| >= B - 4Δ, and

(B - 4Δ)^2 Σ_unbal log p <= W |Σ v_i|^2 <= B W + 2B(B-1) log C_0.

With B > 8Δ, (B - 4Δ) > B/2, so

N_u log P <= 4W/B + 8 log C_0 <= 4 s log(2P)/B + 8 log C_0.

Using P >= 16 and P > C_0^8 gives N_u < 5s/B + 1.

Delete the minority classes, at most 2Δ N_u vertices. The coordinates y were computed before the deletion, so deletion only restricts the established properties. At every unbalanced prime the survivors are in the majority class. That class has size >= B - 2Δ > 2Δ, so it is large and has one residue. Every balanced prime has both classes large, so the residues there are correct with difference ±δ_p.

**(i) Proved.** Survivors with equal y would have equal residues, hence equal signatures, which is impossible. Let m be the mark with the smallest y. Every surviving mark that is good to m satisfies y_m <= y_j = y_m + |d_jm| <= y_m + C sqrt q. Deleting the at most Δ far neighbours of m confines all survivors to [y_m, y_m + C sqrt q]. Mark m itself may later be deleted; the interval still contains the survivors.

**(j) Proved: all remaining distances, and the complementary pairs.** Put u = y_i - y_j.
- At an agreement prime, p | u and p | d. So p^2 | u^2 - d^2.
- At a separating prime, which is necessarily balanced, u ≡ ±δ_p and d ≡ ±δ_p. So p | u^2 - d^2.

By squarefreeness, q h_ij | u^2 - d^2. Since |u^2 - d^2| <= C^2 q, we get u^2 = d^2 whenever h_ij > C^2. This covers far pairs as well: every pair with h > 1 has h >= P > C^2.

Only complementary pairs (h = 1) remain. Each signature has one complement, so these pairs form a matching. If some unbalanced prime exists, no complementary pair survives. Otherwise, for every complementary pair:
- y_i + y_j ≡ ρ_{p,0} + ρ_{p,1} at each p. This is one residue mod q, and the sums lie in an interval of length 2C sqrt q < q. So the sum takes a single value.
- (y_i - y_j)^2 ≡ δ_p^2 at each p, a single residue mod q, and it lies in [0, C^2 q]. So it takes at most 1 + floor(C^2) values.

Sum and absolute difference determine the unordered pair of coordinates. Coordinates are distinct, so the pair determines the marks. Hence there are at most 1 + floor(C^2) complementary pairs. Deleting one endpoint of each finishes the proof.

**Added observation (Proved).** Suppose every prime of q exceeds C^2, and d, d' are admissible for the same signature pair with h > 1. Then d^2 ≡ d'^2 (mod q) and h | d, d', so qh | d^2 - d'^2. Hence |d| = |d'|. A sign-free table on such a modulus is therefore determined by the signatures and δ, except on complementary pairs.

**Computed.** tables.py and bigger.py found no pattern with h > C^2 that has two admissible magnitudes. They cover the moduli of section 2 and q = 59·67·71·83 at C^2 ∈ {36, 49, 58}.

## 4. Attacks and the limiting examples: Computed

**Four-point table** (examples.py). Take q = 11339 = 17·23·29, signatures 000, 011, 101, 110 (Sidon fails: they sum to 0), and δ = (12, 6, 17). The sides are 17, 29, 17, 29 and the diagonals are 46, 46, with 46 < sqrt(q) ≈ 106.5. All local conditions and gcds hold. All four H values are 0. No line embedding exists. This is consistent with X.2's restriction to at least 5 vertices.

**Five-mark table** (examples.py).
- The marks are 0, 11786, 11859, 23645, 28542.
- At q = 23295029 they occupy two residues at every prime. The signatures are distinct, and the table is genuine at C = 6.
- The complementary pairs are (1,4), with d = 23645, and (2,3), with d = 73. Both have gcd 1.
- Replace d_14 by 73. This preserves every local condition. The maximum is d^2/q = 34.97 < 36.
- The triple (1,4,5) has sides 73, 28542, 4897 and H = -1152 q^2. Every other triple is degenerate.
- Deleting mark 1 or mark 4 (1-indexed) leaves an embeddable table; one deletion is enough. As X.4 predicts, the changed edge is the complementary one.

**B_free > B_int on this modulus** (bigger.py). With the same δ and C^2 = 36: B_free = 8 and B_int = 5. At C^2 = 49 and 58: B_free = 8 and B_int = 6.
- The 8-clique is the even-weight subgroup of {0,1}^4.
- Its four complementary pairs each admit the two magnitudes {73, 23645}. Every other edge has a single magnitude.
- Over random choices of the complementary magnitudes, at least 3 deletions (sometimes 4) were needed to make the table line-embeddable at C^2 = 36. At C^2 = 49 and 58 the minimum was 2.

This is consistent with both statements: 8 <= 2·max(4, 5), and four complementary pairs <= 1 + 36. It shows that some deletion is genuinely necessary.

With 40 random δ-vectors on the same primes, B_free = B_int <= 2 in every case.

**Limitation.** The thresholds of the second statement (B > 2277, P > C_0^8) cannot be reached by brute force. The computations test the first statement, the uniqueness observation and the mechanism, not the asymptotic count.

## 5. What X.4 gives and does not give

- **Gives.** Fix C and a dyadic squarefree parent with s >= s_0(C), and fix δ:
  - B_int(q, C, δ) >= B_free(q, C, δ) - O(1 + C^2 + s/B_free);
  - B_free >= cs implies B_int >= cs - O_{c,C}(1);
  - so B_free = o(s) if and only if B_int = o(s) along any sequence of such parents.

  Option (iii) of the question therefore collapses to option (ii). Freeing the signs cannot change the count at scale s. A sign-free construction of order s would be a genuine construction of order s on the same modulus at the same C. This includes Sidon patterns, since the property passes to subsets. At unbalanced primes the cut may become constant after the deletions.
- **Does not give.** It gives no bound on B, no o(s) for worst-case q, and no example of order s. It says nothing about circles, where X.5 shows the separated-square congruence fails.

## Files

- `ref10L_recon/examples.py`: the two limiting examples.
- `ref10L_recon/stmt1.py`, output `stmt1.out`: brute force of B_int, B_free.
- `ref10L_recon/tables.py`, output `tables.out`: uniqueness check and random tables at small moduli.
- `ref10L_recon/bigger.py`, output `bigger.out`: the q = 59·67·71·83 study.
