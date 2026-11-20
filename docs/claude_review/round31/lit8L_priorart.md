# lit8L_priorart: literature status of "V_8 implies Rudin's conjecture"

Downloads are in `round31/lit8L_priorart/dl/<short-name>/` (pdf and pdftotext output). Every arXiv paper listed below
exists only as v1 on arXiv except where stated (version history checked on the /abs/ pages).

## Answer first

1. **No earlier derivation found, now from full texts.** None of the papers I read in full derives sigma(k) = O(sqrt k),
   or sigma(k) <= k^{1/2+o(1)}, from Vojta, Bombieri-Lang, the uniformity conjecture, abc, or any Lang-type conjecture.
   None relates Rudin's conjecture to rational or integral points on moduli spaces of pointed curves, or to integral
   points on one fixed variety. The search was bounded (list in section 4).
2. **The record's comparison sentence is out of date, and it understates something.** "Bombieri-Lang gives
   sigma(k) << k^{5/7}" is correct as a citation of Cilleruelo-Granville, but it is not the best conditional bound in the
   literature. Shkredov-Solymosi (arXiv 2005.11559v1, Theorem 10, eq. (11)) prove |A ∩ P| << |P|^{2/3} for a set A of
   squares and any progression P, under the uniformity conjecture for curves y^k = f(x). They state it as a corollary of
   Bombieri-Lang via Caporaso-Harris-Mazur. Both conditional bounds, k^{5/7} and k^{2/3}, are **weaker than the
   unconditional** k^{3/5+o(1)} of Bombieri-Zannier. The record should say so: the Bombieri-Lang/uniformity route in
   the literature does not even reach the unconditional exponent.
3. **Inverse-sieve literature.** Nothing I read treats sets in an interval of length q^{1/2} that occupy at most two
   classes modulo each prime power exactly dividing one modulus q, or lattice points on short arcs. That literature
   concerns sets in [N] (or [N]^d) with few classes modulo **many primes** p (all p, or all p in a range), not modulo
   the prime factors of one q.

## 1. Squares in arithmetic progressions: what each source proves (item (a))

**READ. Cilleruelo-Granville, arXiv math/0608109v1 (only version).**
* Unconditional state, p.1: "Bombieri, Granville and Pintz [3] proved that σ(k) = O(k^{2/3+o(1)}), and recently Bombieri
  and Zannier [4] have proved that σ(k) = O(k^{3/5+o(1)})."
* Theorem 2 (under Bombieri-Lang): sum_n r_{E+E}(n)^2 << |E|^{11/4} for any set E of squares. The proof starts: "One
  consequence of [8] [CHM] is that there exists an integer B, such that if the Bombieri-Lang conjecture is true then for
  any polynomial f(x) ∈ Z[x] of degree five or six which does not have repeated roots, there are no more than B rational
  numbers m for which ..." [f(m) is a square].
* After Theorem 4: "theorems 2, 3 and 4 show that the Bombieri-Lang conjecture implies σ(k) ≪ k^{4/5}", then a direct
  count over the curves y^2 = x ∏_{j=1}^5 (x + n_j): "and we obtain σ(k) ≪ k^{5/7}. Anyway this upper bound was improved
  unconditionally in [3] and [4]."
* §9 (abc): under Baker's strong form (9.1) r(abc) >> c/exp((log c)^τ), τ = 1/2+o(1), Lemma 1 gives: if there are
  >> sqrt k squares among a+b, ..., a+kb then "a + b ≪ exp(k^{9/2+o(1)})". This bounds the size of the progression,
  not sigma(k).
* No occurrence of "Vojta" in the text (searched). Rudin follows from their congruence Conjectures 9 and 12 (the judge's
  reading of Theorems 11 and 12, v1 numbering, is consistent with the text I have; I did not re-verify the numbering).

**READ. Shkredov-Solymosi, "The Uniformity Conjecture in Additive Combinatorics", arXiv 2005.11559v1 (only version).**
* §3: "Throughout the paper we assume a corollary of the Bombieri–Lang Conjecture by Caporaso, Harris, and Mazur [7] –
  called the Uniformity Conjecture – in its particular form to hyper- and superelliptic curves." Conjecture 4: if
  y^k = f(x) has genus g >= 2, "then it has at most BL_g rational solutions where BL_g depends on genus g only."
* Theorem 10: for A a set of squares and a proper generalised progression H, a bound (10); "In particular, for any
  arithmetic progression P one has |A ∩ P| ≪ |P|^{2/3}. (11)" (Theorem 1 restates (2) = the same bound.)
* Their own comparison (after Remark 9): "Lemma 8, combining with the classical Burgess' method [6], implies that the
  intersection of a set of squares with any arithmetic progression P is, actually, O(|P|^{2/3}) and this is slightly
  stronger than unconditional result from [4] but weaker than the bound on such intersection from [5]." ([4] = BGP,
  [5] = Bombieri-Zannier.) This resolves the inconsistency the meaning referee saw in a snippet: they do **not** claim
  to beat 3/5 for ordinary progressions; their novelty is the generalised-progression bound (10).
* Remark 9: the uniformity conjecture via Lemma 8 "gives the estimate |A| ≪ k^{5/7} which coincides with the bound after
  Theorem 5 from [9]" [CG].
* Theorem 2 (higher powers, k > 3): at most N^{1/2}·exp(-O(log^{1/3} N)) k-th powers in a progression of length N
  (conditional, same hypothesis). This is for k-th powers with k > 3, not squares.
* No square-root bound for squares; no Vojta; no moduli spaces.

**READ. González-Jiménez and Xarles, "On a conjecture of Rudin on squares in arithmetic progressions", arXiv
1301.5122v1 (only version).**
* Abstract: "Rudin's conjecture asserts that Q(N) = O(√N), and in its stronger form that Q(N) = Q(N; 24, 1) if N ≥ 6.
  We prove the conjecture above for 6 ≤ N ≤ 52." The proofs are computational (covering collections, elliptic Chabauty);
  some steps are "under the Parity Conjecture" (§ on Parity Conjecture, line "some of them under the Parity Conjecture").
* Lemma 20(3): "If the Bombieri-Lang conjecture is true, there exists an r such that, for any set of r pairs (q_i, a_i)
  of coprime integers, ∩_{i=1}^r S(q_i, a_i^2) has at most 4 elements." Proof via genus-5 curves and CHM.
* Final remark: "the Bombieri-Lang conjecture implies that, for k ≥ 5, a constant c(k) should exist such that z_I ≤ c(k)
  for all I ⊂ N with #I = k."
* No asymptotic conditional bound on Q(N); no Vojta.

**READ. Xarles, "Squares in arithmetic progression over number fields", arXiv 0909.1642v1 (only version).**
* Abstract: "there exists an upper bound for the number of squares in arithmetic progression over a number field that
  depends only on the degree of the field. We show that this bound is 5 for quadratic fields". This concerns progressions
  consisting entirely of squares, unconditional.
* Corollary 18: o_{d,k}(N) k-th powers in any progression of length N over a field of degree d (via Szemerédi).
* Closing remark: this "can probably be improve[d] to some bound of the type O_{d,k}(N^{1-c_{d,k}}) ... One could even
  guess if c_{d,k} can always be taken equal to (k-1)/k, which will generalize Rudin's conjecture". Bombieri-Lang
  appears only for the separate 2-dimensional question (quadratic polynomials with f(i) square for i <= S).
* No conditional bound on sigma(k).

**READ (relevant parts). García-Fritz, "Quadratic sequences of powers and Mohanty's conjecture", arXiv 1610.05340v1**
(single author, Natalia Garcia-Fritz).
Abstract: "We prove under the Bombieri-Lang conjecture for surfaces that there is an absolute bound on the length of
sequences of integer squares with constant second differences, for sequences which are not formed by the squares of
integers in arithmetic progression." This is Büchi-type, about sequences made entirely of squares; not about sigma(k).

**READ (relevant parts). Dietmann-Elsholtz, "Hilbert cubes in progression-free sets and in the set of squares II",
arXiv 1207.1461v1.** Cites that "the Bombieri-Lang conjecture implies that d is absolutely bounded" (Hilbert cubes of
squares). Not about sigma(k).

**READ (relevant parts). Hegyvári, "Note on the sumset of squares", arXiv 2504.13230v2.** §3: "In 1960 Rudin conjectured,
that in any N−term arithmetic progression the number of square numbers is ≪ √N ... The proof of this conjecture still
seems hopeless today." No conditional derivation.

**SNIPPET / not opened.**
* Bombieri-Granville-Pintz (Duke 1992) and Bombieri-Zannier (Rend. Lincei 2002): not on arXiv; not opened. Their results
  (k^{2/3+o(1)}, k^{3/5+o(1)}) are as reported by CG, GJX, Xarles and SS (all READ, consistent).
* Helfgott-Venkatesh (Roth volume, CUP 2009, 224-234): not on arXiv; no open copy found. Known to me only through
  Green-Harper (READ: they cite "[12, pp 232–233]" for the speculation that quadratics are the only sharp examples for
  the large sieve) and search snippets.
* Pasten, "Powerful values of polynomials and a conjecture of Vojta" (JNT): host refused (403). Per snippets it is about
  Büchi's problem under a Vojta conjecture, not Rudin. SNIPPET.
* Caporaso-Harris-Mazur (JAMS 1997): not opened; used only as cited by CG, GJX, SS.
* Hajdu-Tengely; Alvarado-Koutsianas-Malmskog-Rasmussen-Vincent-West: no relevant paper found by search (NOT FOUND;
  I found no paper with the latter author list on squares in progressions).
* arXiv listing searches (arxiv.org/search, abstracts read for the hits that looked relevant): "squares in arithmetic
  progression" (100 hits, 2017-2026), "Rudin squares progression" (2 hits: GJX, CG), "Vojta squares arithmetic
  progression" (0 hits), "uniformity conjecture squares" (only SS relevant). Abstracts checked: 2606.01780, 2007.14324,
  1708.06951, 2602.03251, 2605.11104, 2610.37423 — none conditional on Vojta/Lang for sigma(k).

## 2. Vojta- or Lang-type derivations; moduli spaces (item (b))

* **Does any paper derive sigma(k) = O(sqrt k) or k^{1/2+o(1)} from a Vojta- or Lang-type conjecture? NOT FOUND.** In
  the full texts read (CG, SS, GJX, Xarles, 1610.05340), the strongest conditional statements are k^{4/5}, k^{5/7}
  (CG, Bombieri-Lang+CHM) and k^{2/3} (SS, uniformity). SS obtain N^{1/2-o(1)}-type bounds only for k-th powers, k > 3.
* **Relation to rational points on Mbar_{0,n}, or to integral points on a fixed variety: NOT FOUND.** The conditional
  arguments in the literature apply Bombieri-Lang/uniformity to **families of curves** (genus-5 curves C_I in GJX,
  y^2 = x∏(x+n_j) in CG, y^k = f(x) in SS), i.e. uniformity over a family, not Vojta on one fixed variety.

## 3. Inverse sieve literature (item (c))

* **READ. Green-Harper, arXiv 1311.6176v1.** Setting (§1): "sets A of integers with the property that the reduction
  A (mod p) occupies at most (p + 1)/2 residue classes modulo p for all sufficiently large primes p." Conjecture 1.5
  and Conjecture 1.7 ("almost equality in the large sieve implies quadratic structure") are about all primes in
  [X_0, X^{1/4}] or all large p. Theorem 2.3 is the larger sieve for a set of primes P.
* **READ. Walsh, "The inverse sieve problem in high dimensions", arXiv 1105.1551v2.** Abstract: "if a big set of integer
  points S ⊆ [0, N]^d, d > 1, occupies few residue classes mod p for many primes p, then it must essentially lie in the
  solution set of some polynomial equation of low degree."
* **READ (abstract and introduction). Walsh, "The algebraicity of ill-distributed sets", arXiv 1307.0259v2.** "every set
  S ⊆ [N]^d occupying ≪ p^κ residue classes for some real number 0 ≤ κ < d and every prime p, must essentially lie in the
  solution set of a polynomial equation of degree ≪ (log N)^C".
* **READ (abstracts). Paredes, arXiv 1901.00562v2; Menconi-Paredes-Sasyk, arXiv 1907.02049v3; Croot-Mao-Pohoata-Yip,
  "A sharp inverse theorem for the quadratic large sieve", arXiv 2607.15311v3.** All in the many-primes setting; the last
  assumes "|A_p| ≤ p/2 + O(1) for every prime p ≪ log N" and concludes A contains exp(c log N/loglog N) points of one
  quadratic m ± x^2.
* **READ (abstract). Aryan, "Distribution of squares modulo a composite number", arXiv 1502.05062v1:** gaps between
  squares modulo squarefree q (statistics of spacing), not bounded-count statements in sqrt(q) windows.
* **Conclusion: nothing found** about sets in an interval of length q^{1/2} occupying at most two classes modulo each
  prime power dividing q, or about lattice points on short arcs.
* *Remark (my own, not a literature claim).* The proof of Green-Harper Theorem 2.3 adapts verbatim to moduli
  p^e || q with |S_{p^e}| <= 2: summing e log p over p^e | x - y for x != y in an interval of length N gives
  |A|^2 (log q)/2 - |A| log q <= |A|^2 log N, so |A| <= log q / ((1/2) log q - log N). This is O(1) for
  N <= q^{1/2-δ} and gives nothing at N = q^{1/2}. So the two-residue statement at exponent 1/2+eps is exactly past the
  larger-sieve barrier; this is consistent with the record's view that TR(1/2-eps) is easy and TR(1/2+eps) is not.

## 4. What the record may say (item (d))

Keep, in this form:

> We found no earlier derivation of Rudin's conjecture, of σ(k) ≪ k^{1/2+o(1)}, or of a two-residue statement of this
> kind from Vojta's, Bombieri–Lang's or the abc conjecture. The conditional bounds we found are σ(k) ≪ k^{5/7}
> (Bombieri–Lang with Caporaso–Harris–Mazur; Cilleruelo–Granville, arXiv math/0608109v1, §3) and
> |A ∩ P| ≪ |P|^{2/3} (the uniformity conjecture for curves y^k = f(x), a consequence of Bombieri–Lang by
> Caporaso–Harris–Mazur; Shkredov–Solymosi, arXiv 2005.11559v1, Theorem 10), both weaker than the unconditional
> k^{3/5+o(1)} of Bombieri–Zannier. Searched in full text: Cilleruelo–Granville, Shkredov–Solymosi, González-Jiménez–Xarles
> (arXiv 1301.5122v1), Xarles (arXiv 0909.1642v1), García-Fritz on Büchi-type sequences (arXiv 1610.05340v1),
> Green–Harper (arXiv 1311.6176v1), Walsh (arXiv 1105.1551v2, 1307.0259v2); not opened: Bombieri–Granville–Pintz,
> Bombieri–Zannier, Helfgott–Venkatesh, Caporaso–Harris–Mazur.

Replace "Bombieri–Lang gives σ(k) ≪ k^{5/7}" (README_VIIIL §1.7 and the note text quoted in §5) by the two-bound sentence
above, or at least add the SS k^{2/3} bound and the remark that both are weaker than the unconditional 3/5.

Do not say "Vojta's conjecture is not known to imply Rudin" as a fact about the literature; say "we found no earlier
derivation" with the list above.

