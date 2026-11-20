# lit8L_verify: adversarial check of the three literature readers (round 31)

## Answer first

I found no refutation. Every quote that the divisors, cg and priorart readers marked READ is in the downloaded full text, and it means what they say. There is one small mislocation (item P2 below) and two omissions that the record should use (C7 and C8). The new searches turned up nothing that changes the status of (L_8):
- No construction of unboundedly many divisors in an invertible class at a fixed alpha > 0.
- No earlier derivation of Rudin's conjecture from Vojta, Bombieri-Lang or abc.

Lenstra 1984 and CHN 2008 / ePrint 2004/339 are still unread. The hosts www.ams.org, eprint.iacr.org and pub.math.leidenuniv.nl failed again: curl got "CONNECT tunnel failed, 403" and WebFetch got "getaddrinfo ENOTFOUND". All Lenstra and CHN wording below is SNIPPET.

Texts I re-opened (all already downloaded by the readers, all v1):
- Hales, arXiv:2410.05030v1: `lit8L_divisors/dl/hales/p.txt`
- Cilleruelo-Granville (CG), arXiv:math/0608109v1: `lit8L_cg/dl/cg_v1/v1.txt`. Its md5 is the same as `lit8L_priorart/dl/cg/paper.pdf`.
- Shkredov-Solymosi, arXiv:2005.11559v1: `lit8L_priorart/dl/shso/paper.txt`
- Gonzalez-Jimenez-Xarles, arXiv:1301.5122v1: `lit8L_priorart/dl/gjx/paper.txt`
- Xarles, arXiv:0909.1642v1: `lit8L_priorart/dl/xarles/paper.txt`
- Green-Harper: `lit8L_priorart/dl/gharper/paper.txt`
- Aryan, arXiv:1502.05062v1: `lit8L_priorart/dl/sqmod/paper.txt`

## 1. Check of the readers' READ findings

### Divisors reader (Hales 2410.05030v1)

- D1 (author, title, v1). **Confirmed READ.**
- D2 (nothing proved at or below 1/4). **Confirmed READ.**
  - §4.1: "the number of divisors of N congruent to r mod S is less than a constant c(α) depending only on α whenever S^α > N and α > 1/4".
  - Thm 4.1: "In addition, c(1/3) = 12".
  - §4.1: "a large portion of [4] is dedicated to proving that c(1/3) = 12 and that c(α) = O((α − 1/4)^−2)".
  - Thm 4.2 (CHN): "Given α > 1/4 there is an integer n0 > 0 such that for all integers n > n0 and s > n^α ... with gcd(r, s) = gcd(s, n) = 1".
  - §4.2: "The case is much less clear when α is closer to 1/4, the threshold of both the Lenstra and Coppersmith et al. methods with significant work only invested in the case α = 1/3."
  - Note: Hales's reference [4] is Lenstra's 1983 Leiden report, not the 1984 Math. Comp. paper.
- D3 (open questions are only at alpha >= 1/3, plus whether Lenstra's bound is sharp for alpha < 1/2). **Confirmed READ.**
  - §4.2: "For α ≥ 1/3 It's currently an open problem as to whether 6 is the best possible".
  - §5 (3): "Is Lenstra's upper bound sharp for values of α < 1/2? Specifically, can an integers N, S, r be found with gcd(N, S) = gcd(S, r) = 1 and S^3 > N so that N has at least 7 divisors congruent to r mod S?"
- D6 and D7 (largest published examples are 6, or 7 counting signed divisors; no construction near 1/4). **Confirmed READ.**
  - §4.2: "N(ℓ) has exactly 6 divisors congruent to 1 mod S(ℓ)".
  - §4.2: "α = log_N(S) ≈ 0.3584".
  - §4.2: "7 divisors congruent to 1 mod S(x) among both positive and negative divisors".
  - Hales has no example with alpha < 1/3.
  - The reader's elementary argument is correct. If n = p_1 ... p_t with every p_i ≡ 1 (mod q), all 2^t divisors are ≡ 1 (mod q), and n > q^t, so alpha < 1/t. This is the same mechanism as prod(1 + iq). It gives at most 2^{floor(1/alpha)}, so it is bounded at each fixed alpha.
- D8 ("16 divisors" is the record's own computation). **Confirmed:** no source I read mentions 16.

Lenstra wording (SNIPPET, from two web searches; full text not reached):
- Abstract: for s > n^{1/3} and gcd(r, s) = 1, "there exist at most 11 positive divisors of n that are congruent to r modulo s".
- The "not bounded by 11, and not even by a polynomial function of log(n^2)" remark is about a non-invertible class: divisors of n^2 that are ≡ n mod 2n.

One search summary claimed that Lenstra says the count is unbounded "for smaller s". A second search shows that this sentence is the non-invertible remark. **No source shows unboundedness in an invertible class below 1/4.** Hales's c(1/3) = 12 (combinatorial constant) and Lenstra's "11 positive divisors" (snippet) do not conflict.

### CG reader (math/0608109v1)

The only arXiv version is v1. I did not re-check that v2 returns 404.

- C1, Conjecture 9 and Theorem 11 (§4 p.8). **Confirmed verbatim:** "Conjecture 9. There exists a constant N such that there are no more than N solutions 0 < x1 < x2 < ··· < xN < x1 + b^{1/2} to the equation x_i^2 ≡ a (mod b), for any given a and b. Theorem 11. Conjecture 9 implies Conjecture 1."
- C2, Conjecture 10. **Confirmed verbatim.**
  - It is stated for "integer d ≥ 1", with |Ω(q)| ≤ d for each prime power q. It reads "for any ǫ > 0 there exists a constant N(d, ǫ) such that for any integer b there are no more than N(d, ǫ) integers n, 0 ≤ n < b^{1−ǫ} with n ∈ Ω(b)."
  - The text then says: "In theorem 9 we proved such a result with the exponent '1 − ǫ' replaced by '1/d − ǫ'. We strongly believe Conjecture 10 with '1−ǫ' replaced by '1/d', analogous to Conjecture 9." Poonen's question is the d = 4 case at exponent 1/2.
  - Theorem 9 is stated for intervals "of length n^{α_d(k)}, where α_d(k) = (1−ε_d(k))/d". So Theorem 9 is about arbitrary intervals, while Conjecture 10 is stated for the initial segment [0, b^{1−ǫ}).
  - The word "critical" does not appear near these statements. The reader is right that it is the record's word, not CG's.
- C3, Conjecture 12 and Theorem 12 (§5 p.9). **Confirmed verbatim.**
- C4, Conjectures 13-15 and Theorem 16. **Confirmed verbatim.**
  - Conj 13: "For any α < 1/2 ... #{(a, b), a^2 + b^2 = n, N ≤ |b| < N + n^α} ≤ C_α".
  - Conj 14: "in an arc of length R^{1−ǫ} is bounded uniformly in R".
  - Theorem 16: "Conjecture 13 implies Conjecture 19".
  - Conj 19 is trivial for α ≤ 1/2, "yet is completely open for any α > 1/2".
- C5, no implication from the arc conjectures to Rudin. **Confirmed READ.** In v1, Conjecture 1 is reached only by Theorem 11 (from Conj 9) and Theorem 12 (from Conj 12). Conjecture 2 is reached by Theorem 1 (from Conj 3, as p -> 4) and Theorem 4 (from Conj 5). The second flowchart ends at Conj 19 and Conj 21.
- C6, Bombieri-Lang remarks in §3. **Confirmed verbatim:**
  - "theorems 2, 3 and 4 show that the Bombieri-Lang conjecture implies σ(k) ≪ k^{4/5}"
  - "we obtain σ(k) ≪ k^{5/7}. Anyway this upper bound was improved unconditionally in [3] and [4]."
  - Theorem 2: "Assume the Bombieri-Lang conjecture. Then Σ r^2_{E+E}(n) ≪ |E|^{11/4}".
- C7, Theorem 10 (§4 p.7). **Confirmed verbatim:** "If x1 < ··· < xk are solutions to the equation x_i^2 ≡ a (mod b), then xk − x1 > b^{1/2 − 1/(2ℓ)}, where ℓ is the largest odd integer ≤ k." The reader's k = 2 counterexample is correct as printed: with b = 3, a = 1, the roots 1 and 2 give a difference of 1, which is not greater than 3^0 = 1. I did not re-run the reader's numerics.
- C8, no R^{1/2} uniform arc conjecture in this paper. **Confirmed (NOT FOUND in v1).**
  - The reader omitted one sentence that matters. It sits right after Conjecture 13 (§5 p.11): "It is simple to prove (5.1) for any α ≤ 1/4 (and Conjecture 13 for α ≤ 1/4 with N ≪ n^{1/2−α}), but we cannot prove (5.1) for any α > 1/4."
  - Here n = R^2, so n^{1/4} = R^{1/2}. CG therefore state that the arc problem is open beyond length R^{1/2} even near the axis. The record may cite this.
- C9, TR(1/2+eps) is the d = 2 case of Conjecture 10. **Confirmed.** It needs the evident translation from [0, b^{1−ε}) to arbitrary intervals, which CG do not write down.

### Priorart reader

- P1, CG conditional bounds. **Confirmed** (see C6). "Vojta" does not occur in CG v1. CG §9 Lemma 1, under Baker's form (9.1) of abc, bounds only A + B, the size of the progression. **Confirmed:** "Then (9.1) implies that A + B ≪ exp(O(T^{9τ/(1−τ)}))".
- P2, Shkredov-Solymosi. **Confirmed, with a small location fix.**
  - Conjecture 4 is in §3: "a corollary of the Bombieri–Lang Conjecture by Caporaso, Harris, and Mazur [7] – called the Uniformity Conjecture".
  - The quote "this is slightly stronger than unconditional result from [4] but weaker than the bound on such intersection from [5]" ([4] = BGP, [5] = BZ) is in the paragraph before Theorem 10. That paragraph attributes O(|P|^{2/3}) to "Lemma 8, combining with the classical Burgess' method". The same bound is then restated as eq. (11) of Theorem 10: "In particular, for any arithmetic progression P one has |A ∩ P| ≪ |P|^{2/3}."
  - So "Theorem 10, eq. (11)" is correct, and the sentence comparing the bound with BZ is in the text just before it.
- P3, GJX. **Confirmed:** "We prove the conjecture above for 6 ≤ N ≤ 52". The Bombieri-Lang statements are of finiteness or constant type (Lemma 20(3), Remark 22, "c(k)").
- P4, Xarles. **Confirmed** (closing remarks): "One could even guess if c_{d,k} ...".
- P6, no earlier conditional derivation of Rudin. I ran three searches of my own:
  - "Vojta conjecture implies Rudin's conjecture squares";
  - "Rudin conjecture ... abc OR Vojta OR Bombieri-Lang O(sqrt N)";
  - "squares in arithmetic progression moduli space OR M_{0,n} OR integral points conditional bound Rudin".

  None found a derivation of Rudin from Vojta, Bombieri-Lang, abc or uniformity, and none found any use of Mbar_{0,n}. The nearest hit is Pasten, "Powerful values of polynomials and a conjecture of Vojta" (not opened; snippet says it is about Büchi-type problems). This is NOT FOUND, not a proof of absence.
- P7, inverse sieve. I spot-checked one quote. Green-Harper is present, and its "(p + 1)/2 residue classes modulo p for all sufficiently large primes p" setting is in `gharper/paper.txt`. I did not re-read the other inverse-sieve papers beyond their headers, which match the titles the reader gives. The reader's larger-sieve remark is their own derivation and should stay labelled as such.

## 2. What the record may now say

### (i) Status of LEN(theta) for theta <= 1/4

Basis:
- READ: Hales 2410.05030v1.
- SNIPPET: Lenstra 1984 and CHN 2008.
- NOT FOUND: any invertible-class unboundedness or boundedness result at or below 1/4.

Proposed sentence:
> Lenstra (1984) proved LEN(alpha) for every alpha > 1/4, with c(alpha) = O((alpha−1/4)^{-2}), and Coppersmith–Howgrave-Graham–Nagaraj (2008) obtained O((alpha−1/4)^{-3/2}) for n > n_0(alpha) (as restated in Hales, arXiv:2410.05030v1, §4.1, Thms 4.1–4.2). Hales calls 1/4 "the threshold of both the Lenstra and Coppersmith et al. methods" (§4.2). We found no result, in either direction, for invertible classes at alpha ≤ 1/4; Lenstra's unboundedness remark (divisors of n^2 that are ≡ n mod 2n) concerns a non-invertible class (Lenstra 1984 seen through snippets only). The published examples we read have 6 divisors (Cohen; Hales §4.2), or 7 counting signed divisors, all at alpha > 1/3. The elementary family prod_{i≤4}(1+iq) gives 16 divisors at exponents just below 1/4 (our computation); families of this type give at most 2^{⌊1/alpha⌋} at a fixed alpha. Whether LEN(1/4) or LEN(1/4−kappa) holds is open as far as we could establish.

Replace "no construction ... is known, and the best known gives 16" with the sentence above.

### (ii) Cilleruelo-Granville statements and numbering

Basis: READ, arXiv math/0608109v1, the only arXiv version. The CRM 2007 published version was not read.

Proposed sentences:
> (a) TR(1/2+eps) is the case d = 2 of Cilleruelo–Granville's Conjecture 10 (arXiv math/0608109v1, §4), which they state for initial segments [0, b^{1−ε}) and every ε > 0 (the translation to arbitrary intervals is evident). Their Theorem 9 proves the exponent (1−ε_d(k))/d for arbitrary intervals, and they write that they "strongly believe Conjecture 10 with '1−ǫ' replaced by '1/d'". TR(1/2) is their Conjecture 9, which implies Rudin's conjecture (their Theorem 11).
>
> (b) CG also remark (§3) that Bombieri–Lang implies σ(k) ≪ k^{4/5} (Theorems 2–4) and, by a direct count on genus-2 curves, σ(k) ≪ k^{5/7}, and that both are beaten unconditionally.
>
> (c) CG's arc conjectures (Conjectures 13–15, arcs of length R^{1−ε}) are linked only to Conjecture 19 (Theorem 16), not to Rudin's conjecture. CG note: "It is simple to prove (5.1) for any α ≤ 1/4 ..., but we cannot prove (5.1) for any α > 1/4" (that is, arcs of length R^{1/2} near the axis).
>
> (d) This paper contains no uniform conjecture for arcs of length R^{1/2}; at that scale it has only Theorem 13 (at most k points on arcs of length R^{1/2−1/(4[k/2]+2)}).

Also:
- Drop the word "critical".
- Cite CG Theorem 10 as stated for k ≥ 3, since it fails for k = 2 as printed.
- If the record's uniform bound "U / M(C)" is attributed to CG, it needs a different source, perhaps CG's "Close lattice points" [15], which nobody has read.

### (iii) Novelty of "Vojta on Mbar_{0,8} implies Rudin"

Basis:
- READ: CG, SS, GJX and Xarles in full.
- NOT FOUND: any earlier derivation, in three searches of my own plus the priorart reader's.
- Not opened: BGP 1992, BZ 2002, CHM 1997.

Proposed sentence:
> We found no earlier derivation of Rudin's conjecture (or of σ(k) ≪ k^{1/2+o(1)}) from Vojta's conjecture, Bombieri–Lang, the uniformity conjecture or abc, and no use of Mbar_{0,n} in this problem. Earlier conditional results are weaker than the unconditional bound σ(k) ≪ k^{3/5+o(1)} (Bombieri–Zannier): Bombieri–Lang gives σ(k) ≪ k^{5/7} (Cilleruelo–Granville, §3), and the Caporaso–Harris–Mazur uniformity conjecture gives |A ∩ P| ≪ |P|^{2/3} (Shkredov–Solymosi, arXiv:2005.11559v1, Thm 10 eq. (11), which they call "weaker than" Bombieri–Zannier's bound). The implications to Rudin's conjecture that are in the literature come from analytic or combinatorial conjectures: Rudin's Λ(p) conjecture (CG Thm 1, giving k^{1/2+ε}), and CG Conjectures 9 and 12 (Thms 11, 12). Searched: CG, Shkredov–Solymosi, González-Jiménez–Xarles, Xarles, García-Fritz, plus arXiv and web searches; not opened: Bombieri–Granville–Pintz, Bombieri–Zannier, Caporaso–Harris–Mazur, Pasten.

### (iv) Inverse-sieve literature

Basis:
- READ by the priorart reader. I spot-checked only the Green-Harper quote.
- NOT FOUND: any result about two classes modulo each prime power dividing a single q, or about short arcs.

Proposed sentence:
> The inverse large-sieve literature we read (Green–Harper; Walsh; Croot–Mao–Pohoata–Yip; and others listed in the referee report) concerns sets in [N] with few residue classes modulo many primes; we found nothing on sets occupying ≤ 2 classes modulo each prime power dividing a single modulus q in an interval of length q^{1/2}. (Our own adaptation of the larger sieve gives boundedness for intervals of length q^{1/2−δ} and fails at exponent 1/2.)

## 3. Reader claims I could not confirm

1. All Lenstra 1984 wording: the "11 positive divisors" statement, the non-invertible n^2 example, and the absence of any remark on invertible classes at alpha <= 1/4. This is SNIPPET only. ams.org and the Leiden hosts were unreachable again.
2. CHN 2008 / ePrint 2004/339 remarks at alpha <= 1/4. SNIPPET; I saw only Hales's restatement. eprint.iacr.org and ams.org were unreachable.
3. Cohen 1982-83. Not opened; known only through Hales §4.2.
4. That arXiv math/0608109v2 returns 404. I did not re-check it. The CRM 2007 published numbering is unread.
5. The reader's numeric check of CG Theorem 10 for k = 3..6 (script `lit8L_cg/scripts/t10.py`). Not re-run. I checked the k = 2 counterexample by hand.
6. The arXiv search-listing results (divisors reader D5, priorart searches). Not re-run; I ran my own web searches instead.
7. Inverse-sieve papers other than Green-Harper. Not re-read beyond their headers. The Helfgott-Venkatesh, BGP, BZ, CHM and Pasten texts remain unread by anyone.
8. The "16 divisors at alpha = 0.2317 and 0.2411" computation (judge §1.3). This is the record's own computation and was not re-run here.

To finish items 1 and 2, the session network must allow www.ams.org, eprint.iacr.org or pub.math.leidenuniv.nl.
