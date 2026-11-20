# Referee report on Part IX, section IX.3, second half (block factor2)

Refereed text: `astra_free_IX/answer.md`, lines 2071–2085. That is the nonpair characters and (IX.7), the union bound and offsets, the real correction, the moment construction and the spread version, the paragraph "which stronger axiom fails", and the COMPUTED paragraph. The first half of IX.3 (Paley set-up, the right inverse U, (IX.4)–(IX.6), the same-level deletion and the pair estimate) belongs to another block. I use its statements as given and say where my verdicts depend on them.

Scripts (my own, exact arithmetic; none of Astra's code is imported):
`ref9_factor2/q11_check.py` (output `q11_out.txt`) and `ref9_factor2/paley_checks.py` (output `paley_out.txt`), both under `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/`.

## 0. Verdicts

| item | verdict |
|---|---|
| (1) The three nonpair bounds from VIII.3, used within their proved range (n ≥ 4) | **Proved** (correct) |
| (1) (IX.7) in the three ranges, with a_ρ = (1−ρ)/2 − o(1) | **Proved** (correct). ρ < 1 is used only at n ≈ m. a_ρ → 0 as ρ → 1. |
| (2) Density, avoided intervals, count, two cases, total failure → 0 | **Proved** (correct). Exponents recomputed in §2. |
| (3) Real correction t_j | **Proved** (correct). Exact identity checked. |
| (3) Moment construction for the remaining monomials, zero windings, compatibility with the witnesses | **Proved**, subject to the range-of-d caveat in §3 |
| (3) Spread version: even winding differences | **Proved** and **Computed** |
| (4) "These witnesses fail e-all at a large cutoff" | **Proved**, with small wording corrections (§4). The conclusion is about these witnesses (more generally, about every S-unit witness family). It is not about every witness family for the model. |
| (5) The q = 11 computation | **Computed**. Reproduced when U is taken from FLINT's HNF transform (max entry 9; 330 and 90 identities; (23, 173, i^3)). My own right inverse gives max entry 18 and also passes. |
| (6) What is capped | P3-char at **every** cutoff, with X-independent global witnesses. P3-all (e-all) is not capped at large cutoffs. P3-all\* at polynomial cutoffs is not touched. |

Overall: the second half of IX.3 is correct. Given the first half and the VIII.3/VIII.10 moment construction, it proves this:

> **Theorem (IX.3, corrected statement).** Fix C > 0, ρ ∈ (0.99, 1) and a moment-construction parameter set as in VIII.3 or VIII.10 with prime exponent d (any d admissible there; all such d exceed 3/2). For every prime q ≡ 3 mod 4 with q ≥ q_0(C, ρ, d) there is an angle model with m = (ρ/2 + o(1))q marks and ω(N') = q − 1, N' = ∏ p_j², p_j ∈ [q^d, 2q^d]. It also has global witnesses r_j ∈ Q(i)^*, r_j r̄_j = p_j. The model and witnesses lie in F^char_X(m, C) of `round5/round6_sharp/sharp.md` Def. 1.1 for **every** X ≥ 3, with the same tuple (p, e, a, φ, θ_0, δ, k) and with r_j(q) the reductions of r_j (with the singular factor deleted at parent primes). Here W/(m log m) = 4d/ρ + o(1). Diagonalising ρ → 1 gives m/ω(N') → 1/2 and the stated limsups: 4d with d → 23/12 from VIII.10 alone (P_+ ≤ m^{23/12+o(1)}, W ≤ (23/3 + o(1)) m log m), and 11/6, 22/3 only if IX.4 is accepted. For X ≥ X_0(model), the same model and witnesses are **not** in F^all_X.

The threshold q_0 grows as ρ → 1 (roughly log q_0 ≫ 1/((1−ρ)d)). It does not depend on any radius, since the object is an angle model. No exceptional sets, discarded gcds or drifting phases occur in this block. The common phase is the explicit −(π/q)Σε plus constants, carried exactly.

## 1. Nonpair characters and (IX.7)

**Definitions.** S ⊂ 𝒜 = {i ∈ F_q^* : K_{i0} = 1} is the mark set, with |S| = m = ρ(q−1)/2 + O(q^{2/3}) < q/2. A character is c ∈ Z^S with Σc = 0, c ≠ 0. Then v = K_S^T c restricted to I = F_q^*, n = ‖c‖_1 (even), h = ‖v‖_1, and b_0 = log L_0/log P_0 with L_0 = max p_j ≤ P_0 e^σ. So b_0 ≤ 1 + σ/log P_0 = 1 + o(1). The pairs are exactly n = 2, so the nonpairs are exactly n ≥ 4. This includes 2(e_a − e_b), which has n = 4 and h = 2q + 2.

**The three bounds (Proved; refereed form in `round31/ref8_fakecore.md` §2.5).**
* (B) (n+1)h ≥ q‖c‖_2² + n q/2. Proof: ‖v‖_2² = q‖c‖_2², from KK^T = qI − 𝟏𝟏^T, Σc = 0, and v_0 = Σc_i K_{i0} = 0. Also y² ≤ (n+1)y − n for integers 1 ≤ y = |v_j| ≤ n, and |supp v| ≥ (q+1)/2. The support bound uses only that c is supported on fewer than q/2 marks; when c ≡ 0 mod q one divides out the q-power, which preserves the support. So (B) holds for all characters.
* (A) h ≥ 1.2q for n ≥ 4. This follows from (B) and ‖c‖_2² ≥ n: h ≥ (3q/2)·n/(n+1) ≥ (3q/2)(4/5).
* (C) h ≥ qn/m. This follows from n h ≥ ‖v‖_∞ h ≥ ‖v‖_2² = q‖c‖_2² ≥ q n²/m.

IX.3 uses all three only for n ≥ 4, which is inside their proved range (the "s ≥ 4" condition of VIII.3). Computed: there are no violations among about 1.6·10⁶ characters at q = 43, 67, 103, for S = 𝒜 and for a random S of density 0.9. Every n = 4 vector with entries ±1 and with pattern (2, −1, −1) was tested, plus random ones (`paley_out.txt`). These numerics do not stress the bounds near n = m, where true h is far above them. The proof carries the claim.

**(IX.7) recomputed (Proved).** Put a_ρ := (1−ρ)/2 − η for any fixed η > 0. For q ≥ q_0(ρ, η, d):

* *n ≤ q/10:* h − b_0 n ≥ 1.2q − b_0 q/10 ≥ 1.09q ≥ (1 + a_ρ)q, since a_ρ < 0.005.
* *q/10 < n ≤ m:* by (B) with ‖c‖² ≥ n, h − b_0 n ≥ (3q/2)(1 − 1/(n+1)) − b_0 n. Its n-derivative is (3q/2)/(n+1)² − b_0 < 0 for n ≥ q/10, so the minimum is at n = m:
  3q/2 − 15 − b_0 m ≥ q(1 + (1−ρ)/2) − q^{2/3} − 16.
* *n > m:* by (B) with ‖c‖² ≥ n²/|supp c| ≥ n²/m, h ≥ (n/(n+1))(qn/m + q/2). Hence
  h − b_0 n ≥ n(q/m − b_0) + q/2 − 3q/(2m).
  This is increasing in n because q/m − b_0 = 2/ρ − 1 + o(1) > 0. At n = m it equals 3q/2 − b_0 m − O(1/ρ), the same slack as before.
* *Second inequality:* by (C), b_0 n/h ≤ b_0 m/q = ρ/2 + o(1) ≤ 0.51. So h − b_0 n ≥ 0.49h. This needs only m < 0.51q, not ρ < 1.

**Where ρ < 1 is used, and how a_ρ behaves.** Only in the middle and last ranges at n ≈ m, where the slack above q is exactly (1−ρ)/2·q − O(q^{2/3}). So a_ρ = (1−ρ)/2 − o(1) → 0 as ρ → 1. At ρ = 1 the first inequality degenerates to h − b_0 n ≥ q − o(q), and the union bound below fails, since P_0^{−aq/2} must beat 2^{5q}. Therefore m/ω(N') → 1/2 is reached only along a diagonal ρ ↑ 1, with q_0(ρ) → ∞. Astra's (IX.3) is a limsup over such a diagonal, which is legitimate. The point to note is that b_0 = 1 + o(1) here, because the witnesses' collision cost uses the physical primes themselves. In VIII.3 the auxiliary primes gave b_0 < 0.51. That is why the slack is only (1−ρ)/2 instead of a constant.

## 2. The probability estimate (Proved; exponents recomputed)

* *Density.* Conditional on all other offsets, c_iε_i/2 with c_i ≠ 0 is uniform on an interval of length |c_i|Δ/(100L_0) < 1. Its density modulo one is at most 100L_0/Δ. Convolution with independent terms does not increase this, so the density of c·ε/2 mod 1 is ≤ 1 + 100L_0/Δ = O(1 + L_0/Δ).
* *The avoided intervals are sufficient.*
  * With zero windings, c·δ = (π/2)c·ε. The e-char left side at unit η is |sin((π/2)(c·ε/2 − arg η/π))|, whose zeros mod 1 lie in {0, 1/2}.
  * The ordinary gaps for v(c) are |sin(⟨v,φ⟩ − sπ/4)| with ⟨v,φ⟩ = (π/4)c·ε. Their zeros are also at c·ε/2 ∈ {0, 1/2} mod 1.
  * If y is the distance of c·ε/2 to {0, 1/2} + Z, then each sine is ≥ y.
  * The worst right-hand side is max(2^{−1/2}m_{c,η}, 1)·e^{−w/2} ≤ L_0^{n/2}P_0^{−h/2}. This uses m_{c,η} ≤ √2 L_0^{n/2} (proved in the first half: m divides A_c − ηĀ_c ≠ 0, which is divisible by 1+i, uniformly in X) and w ≥ h log P_0.
  * So radius r = 2L_0^{n/2}P_0^{−h/2} suffices for all four units, both ordinary gaps and every cutoff, with a factor 2 to spare.
* *Count.* The number of v ∈ Z^{q−1} with ‖v‖_1 = h is ≤ 2^{q−1}·binom(h+q−2, q−2) ≤ 2^q·binom(h+q, q) ≤ 2^{h+2q}. v determines c, since qc = (Kv)_S by the frame identity, the zero column contributing nothing.
* *1/Δ.* By full coverage W = 2Σ_j log p_j ≤ 2(q−1)(log P_0 + σ). So 1/Δ = 10e^{W/4}/min(1,C) ≤ 10(P_0e^σ)^{(q−1)/2}/min(1,C).
* *Total.* Write λ = log P_0 ≥ d log q. A single c fails with probability ≤ 4r(1 + 100L_0/Δ) ≤ 16000 min(1,C)^{−1}·L_0(P_0e^σ)^{(q−1)/2}·P_0^{(b_0n−h)/2}.
  * For 1.2q ≤ h ≤ 3q (at most 2q values of h, count ≤ 2^{5q}), the first bound of (IX.7) gives a total ≤ exp(−(a_ρ/2)qλ + 5q log 2 + λ(1+o(1)) + qσ/2 + O_C(1)). Here qσ = O(q^{2−d} log q) = o(qλ). The total tends to 0.
  * For h > 3q, use q < h/3 and the second bound. Each h contributes ≤ exp(−h(0.245λ − λ/6 − (5/3)log 2 − σ/6) + λ + O_C(1)) = exp(−h(0.078λ − 1.16 − o(1)) + λ + O_C(1)). The sum over h > 3q is ≤ exp(−0.23qλ + O(q)) → 0.

So the total failure probability tends to zero as q → ∞, for fixed C, ρ and d. Astra's sketch is correct.

## 3. Real correction, moment construction, spread version

* **t_j (Proved, Computed).** For x ∈ S, Σ_{j∈I} K_{xj}K_{ij} = qδ_{xi} − 2, because the j = 0 term is K_{x0}K_{i0} = 1. Hence 2K_S t = ε − (2/q)(Σε)𝟏, and v·t = c·ε/2 for every character. The point phases are δ_x = (π/2)ε_x ∈ [0, (π/2)Δ] ⊂ [0, Ce^{−W/4}], with common phase θ_0 = −(π/q)Σε + const. Definition 1.1(b) leaves θ_0 free, so "unrestricted common phase" is exactly right. |t_j| ≤ Δ/2. Checked exactly with rational ε at q = 43 and 67.
* **Remaining monomials (Proved, with one caveat).** The VIII.3 noncharacter analysis uses t only through |t_j| ≤ Δ: the deterministic absorption hΔ for h ≤ H, and conditioning on everything except c_0, ζ or X_k otherwise. Every added term has constant K_S-image, so it contributes zero to characters and needs no windings. P3-char asks nothing of noncharacters beyond the ordinary gaps (c). The witnesses enter only characters, and through u only, and u = 0 when the windings are zero. So the claim holds.
  * Caveat (wording): IX.3 says "d > 3/2". The character part needs only d > 3/2. The moment construction needs a d that VIII.3/VIII.10 (or IX.4) proves admissible, i.e. d = 1.99, d = 2 − ξ with 0 < ξ < 1/12, or IX.4's range. The correct quantifier is "d admissible for the moment construction used". Astra's text does say IX.4 supplies the height.
* **Spread version (Proved, Computed).** For x ∈ S:
  (Kb)_x = (1/12)Σ_{j∈I}K_{xj} + (1/3)(τ_1 + Kz)_x = −1/12 + (1/3)(τ_1 + Kz)_x,
  where Σ_{j∈I}K_{xj} = −K_{x0} = −1, and the column-zero terms are constant on S. τ_1 + Kz is constant mod 3 on S (VIII.4). So (Kb)_x − (Kb)_{x'} ∈ Z, and the windings k = 2Kb have even differences.
  * Moreover v·b = c·k/2 ∈ Z for characters. So the background moves neither the e-char left side (which sees only δ = (π/2)ε) nor the two centres mod 1.
  * Exact check at q = 43, 67, with a random T of size 3 solved over F_3: all differences are integers.

## 4. "Which stronger axiom fails" (Proved, with corrections)

Record statement (sharp.md (1.2)): for every v ∈ Z^r∖{0} and s ∈ Z/4,
|sin(⟨v,φ⟩ − sπ/4)| ≥ ν_s m_{(v,s)} e^{−w(v)/2}, where m_{(v,s)} = ∏_{q≤X, q∤N, q odd} q^{a_{q,s}(v)} and a_{q,s}(v) = max{a ≤ A_q : ∏(r_j/r̄_j)^{v_j} ≡ i^s mod q^a}.

The argument, made precise:
1. **Some row has J_j > p_j.** F = I + 2V, so F_jj is odd and p_j | J_j = ∏_l p_l^{|F_jl|}. With distinct primes, J_j = p_j forces F_j = ±e_j. If every row were ±e_j, F would be diagonal with entries ±1, and H = DF would have rows of the same ℓ¹ length as the rows of D. On I, the rows of D have length q+1 (pair distance q+1, column zero equal) and those of H have length 2. Contradiction.
2. Write α = ∏_{F_jl>0}π_l^{F_jl}∏_{F_jl<0}π̄_l^{−F_jl}. Then r_j/r̄_j = (−1)^{u_j}α/ᾱ exactly in Q(i), and |α|² = J_j.
   * α/ᾱ is not a root of unity: π_j or π̄_j divides α but not ᾱ. Hence arg α/π is irrational, and Im α^n ≠ 0.
3. For n even the unit (−1)^{nu_j} is 1. For odd q ∤ N, the collision α^n/ᾱ^n ≡ 1 mod q^a is equivalent to q^a | Im α^n. So at X ≥ |Im α^n|, m_{(ne_j, 0)} equals |Im α^n| with its 2-part and its parent-prime part removed.
   * Primes p_l with F_jl ≠ 0 do not divide Im α^n: exactly one of π_l, π̄_l divides α^n − ᾱ^n.
   * The other parent primes and 2 are frozen on n ∈ 4 + TZ, where T is the order of α/ᾱ in (Z[i]/p^{e+1})^* (respectively mod 2^{g} with g large) and e is the valuation at n = 4. Then the valuation at n equals the valuation at 4.
   * Irrational rotation gives |sin(n arg α)| ≥ 1/2 infinitely often along that progression.
   * So m_{(ne_j,0)} ≥ J_j^{n/2}/(2c(N', j)). e-all with s = 0 requires this to be ≤ p_j^{n/2}, which fails for large n.
4. m_{(v,s)} is nondecreasing in X. So failure at X_0 implies failure for all X ≥ X_0.

**Corrections.**
* "n a multiple of four" can be any even n. This is harmless.
* "≫_{N'}" means a constant depending on (N', j, α), i.e. on the model.
* The bad cutoff in the proof is |Im α^n|, which is astronomically large.

**Scope.** The proved conclusion is that *these witnesses* fail e-all for all X ≥ X_0(model). The same proof shows a little more (Proved): *every* witness family r_j in the S-unit group ⟨i, π_l, π̄_l : p_l | N⟩ with r_j r̄_j = p_j satisfies e-all at every cutoff only if each r_j is an associate of π_j or π̄_j. Write r_j = unit·∏π_l^{x_l}π̄_l^{y_l}. Then J = ∏p_l^{|x_l−y_l|}, and the argument above applies.

It does **not** show that "no witnesses for this model satisfy e-all". Witnesses that vary with the modulus are not ruled out, and neither are global non-S-unit ones. That is the subject of IX.7, not of this paragraph.

**Heuristic.** These witnesses probably fail e-all already at polynomial cutoffs. For a row with J_j ≫ p_j the ratio is "random" mod ℓ, and a singleton collision at a prime ℓ > √(2p_j) breaks (1.2). The expected number of such ℓ ≤ X is about 4 log(log X/log √p_j). Evidence: at q = 11 the singleton test already fails at ℓ = 13 (FLINT U) and ℓ = 17 (my U). Nothing in IX.3 depends on this.

## 5. The q = 11 computation (Computed, reproduced)

Data: q = 11 and marks S = 𝒜 = {2, 6, 7, 8, 10}. All columns have both endpoints. Each mark has an odd exponent exactly at its own label, and the rows of D have length 12. Primes are recovered from her script: the first ten primes ≡ 1 mod 4 from 101 (101, …, 197) on labels 1..10, with π = x + iy, x minimal. Windings are k_x = 2x² + 4.

Results:
* **Right inverse.**
  * FLINT's HNF transform, the same library routine she used, gives max|U| = 9, max|V| = 13 and max|F| = 26. This reproduces her figure.
  * My own unimodular column reduction gives a different integral right inverse: max|U| = 18, max|V| = 24, max|F| = 48.
  * Both satisfy DU = I and DF = H exactly. "Largest entry nine" is a property of that particular U, not of D.
* **Exact identities in Q(i).** The witnesses satisfy N(r_j) = p_j and r_j/r̄_j = (−1)^{u_j}∏γ_l^{F_jl}. The marks satisfy ρ_x/ρ_y = γ_x/γ_y for all 10 pairs. This is stronger than the modular checks.
* **Modular checks.** 330 exterior pair identities (11 primes ℓ ≤ 37, powers 1–3, 10 pairs) and 90 same-level identities (with γ_l deleted, at p_l, p_l², p_l³) all pass, for both choices of U. Owners are never same-level partners.
* **e-all failure.**
  * With FLINT's U, 7 of 10 rows have J_j > p_j. The singleton test fails at (ℓ, p_j, s) = (23, 173, 3), and 23² = 529 > 346 = 2·173. This reproduces her claim.
  * Smaller failures also exist: (13, 101, 0), with 13² > 101. Her script searched only from ℓ = 19.
  * For rows 0–2 the monomial 4e_j already violates e-all once X ≥ |Im α^4|.
  * With my U, 5 of 10 rows have J_j > p_j, and the first failure is (17, 157, 0). So which (ℓ, p) fails depends on U, but the failure itself does not.

## 6. What IX.3 caps, in the record's vocabulary

**Capped (Proved, conditional on the first half of IX.3).**
* P3-char = F^char_X(M, C) of sharp.md Def. 1.1 (a)–(d) with (e-char) and the same-level clause, **at every cutoff X simultaneously**, with X-independent global witnesses. That is, it caps residues at all moduli, not only at moduli ≤ 3M^A with A fixed.
* This fills the P3-char half of the sharp.md table row "any of these, moduli up to N^A: not proved" (§7), at growth scale with an A-independent constant: W ≤ (23/3 + o(1))m log m from VIII.10, or (22/3 + o(1)) if IX.4 holds.
* It supersedes Theorem 5.1's 33 max(A, 2) for large A.
* The §7 obstruction to coherent designs is avoided because the auxiliary directions are the physical primes themselves (r_j ≡ π_j up to squares, F ≡ I mod 2). This matches the sharp.md §0 item 4 remark that only the genuine π_j up to squares reproduce the Legendre pattern. The price is b_0 = 1, paid for by the slack h − n ≥ (1 + a_ρ)q.
* Since F^char_X ⊂ A_X, it also caps (R) and A_X at all X.
* The fakes have m ≍ ω(N') and W ≍ m log m. So the class cannot prove M = o(log R/loglog R), and in particular it cannot prove M = o(ω(N')).

**Left uncapped.**
* e-all (F^all_X) at large cutoffs. These witnesses, and every S-unit witness family other than the genuine primes up to associates, fail it.
* P3-all\* (W_min^all = O(M log M) at moduli ≤ 3M^A). IX.3 proves nothing either way; heuristically its witnesses fail there too.
* The association of prime norms with their genuine Gaussian arguments.
* The other items of README_VIII §3(e) not addressed by IX.3.
