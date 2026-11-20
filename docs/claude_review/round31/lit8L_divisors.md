# Literature status of LEN(theta) at and below 1/4 (block: divisors)

## Answer first

1. **One full text read: arXiv:2410.05030v1.** Its author is Jonathon Hales and its title is "Divisors in Residue Classes Revisited" (math.NT, submitted 7 Oct 2024; v1 is the only version). It proves nothing at or below exponent 1/4. It describes the regime near 1/4 only as "much less clear". It poses as open only questions at alpha >= 1/3 (is 6 optimal; can 7 divisors be reached with S^3 > N) and whether Lenstra's bound is sharp for alpha < 1/2. It gives no construction near or below 1/4, and it says nothing either way about whether the count is unbounded below 1/4. **READ.**
2. **None of the other sources could be opened.** That covers Lenstra 1984, Coppersmith-Howgrave-Graham-Nagaraj (CHN: Math. Comp. 2008, IACR ePrint 2004/339, IBM RC23466), Cohen 1982-83, the OEIS entry, and Hales's OSU thesis. Every host refused: eprint.iacr.org, ams.org, pub.math.leidenuniv.nl, scholarlypublications/openaccess.leidenuniv.nl, dominoweb.draco.res.ibm.com, cr.yp.to, library.slmath.org, etd.ohiolink.edu, api.semanticscholar.org and oeis.org all gave proxy 403 or "host not allowed". The arXiv export API was also blocked. Only arxiv.org itself, including its search page, could be reached. Statements about those papers below are **SNIPPET** or **NOT FOUND**.
3. **(b), the key question.** I found **no** construction, in what I read or in the search results, of unboundedly many divisors in an invertible class at any fixed exponent alpha > 0. The best families named in the literature I read give 6 divisors (Cohen; Hales's new examples), at alpha just above 1/3 or about 0.3584. The record's "16 divisors near 1/4" is the record's **own** elementary family (2^t divisors with t = 4), not a literature result. I found no source that proves c(alpha) bounded on [1/4 - kappa, 1/4], and none that proves c(alpha) -> infinity as alpha -> 1/4 from below. Status: **open as far as could be established, with one full text read (Hales) and the rest at snippet level.** This does not refute (L_8).

## Sources and what each establishes

### arXiv:2410.05030v1, J. Hales, "Divisors in Residue Classes Revisited". READ (full PDF, 14 pp. text)

Downloaded to `lit8L_divisors/dl/hales/p.pdf`. Author, title and version are from the PDF header and the abstract page: "arXiv:2410.05030v1 [math.NT] 7 Oct 2024", "JONATHON HALES". The paper is mainly about generalizing Lenstra's algorithm to Z[i]-type rings and Z[x]. Section 4, "A Note on Upper and Lower Bounds", is a survey section.

- **Upper bounds, §4.1.** Hales restates Lenstra's result: "the number of divisors of N congruent to r mod S is less than a constant c(alpha) depending only on alpha whenever S^alpha > N and alpha > 1/4." Note that he states it with gcd(N,S) = gcd(S,r) = 1.
- **Theorem 4.1 (Lenstra's combinatorial theorem), as stated by Hales.** "... alpha in R, alpha > 1/4 ... Then |D| <= c(alpha) ... In addition, c(1/3) = 12."
  - The value 12 is for the combinatorial set-system constant. The 11 in Lenstra's divisor theorem is a snippet only, so I cannot reconcile the two from the full text.
- **Lenstra's c(alpha), as reported by Hales.** "the upper bound c(alpha) is not given constructively, and a large portion of [4] is dedicated to proving that c(1/3) = 12 and that c(alpha) = O((alpha - 1/4)^{-2})."
- **Theorem 4.2 (CHN), as restated by Hales.** "Given alpha > 1/4 there is an integer n_0 > 0 such that for all integers n > n_0 and s > n^alpha and 0 < r < s < n with gcd(r,s) = gcd(s,n) = 1, the number of divisors d of n with d ≡ r mod s is bounded above by 2 + 2 pi alpha/(alpha-1/4)^{3/2} + 4 alpha/(alpha-1/4)."
  - The formula is reassembled from garbled PDF text extraction. It matches the earlier referee's snippet.
  - Hales also says the CHN bound holds "only ... [for] n sufficiently large (depending on alpha), showing that c(alpha) = O((alpha - 1/4)^{-3/2})."
- **Lower bounds, §4.2.** "The case is much less clear when alpha is closer to 1/4, the threshold of both the Lenstra and Coppersmith et al. methods with significant work only invested in the case alpha = 1/3."
  - Cohen [1] (as reported): for l >= 3, N(l) = (2l+1)(l^2+1)(l^2+l+1)(2l^2-l+1)(2l^2+l+1) and S(l) = 2l^3+l^2+2l give S(l)^3 > N(l), and "N(l) has exactly 6 divisors congruent to 1 mod S(l)."
  - "For alpha >= 1/3 It's currently an open problem as to whether 6 is the best possible."
  - A new example: N = 104254876089000, S = 105787, with 6 divisors ≡ 1 and alpha ≈ 0.3584.
  - A new polynomial family with "7 divisors congruent to 1 mod S(x) among both positive and negative divisors".
  - The case alpha = 1/2: "c(1/2) = 2 is achieved by infinitely many tuples".
- **Open problems, §5, question (3).** "Can new examples of integers with many divisors in a single residue class? Is Lenstra's upper bound sharp for values of alpha < 1/2? Specifically, can an integers N, S, r be found with gcd(N,S) = gcd(S,r) = 1 and S^3 > N so that N has at least 7 divisors congruent to r mod S?" (The phrasing is the author's.)
- **Not in the paper.** There is no statement about alpha <= 1/4 being bounded, unbounded or open. There is no construction near 1/4, no 2^t family and no "16". Nothing is said about non-invertible classes.

### Lenstra, "Divisors in residue classes", Math. Comp. 42 (1984) 331-340. Full text NOT READ (all hosts blocked)

- **SNIPPET.** Theorem for s > n^{1/3}, gcd(r,s) = 1: at most 11 positive divisors. Existence of c(alpha) for alpha > 1/4. The order O((alpha-1/4)^{-2}) is confirmed second-hand by Hales (READ, above).
- **SNIPPET (earlier referee) / NOT FOUND (here).** The non-invertible example: for odd n, the divisors of n^2 that are ≡ n mod 2n correspond to the divisors of n, so their number is unbounded. I could not find the exact sentence. Nor could I find any sentence by Lenstra saying the count is unbounded for alpha < 1/4, or that the alpha < 1/4 case is open. The task asks for such a sentence; it is **NOT FOUND**.

### Coppersmith, Howgrave-Graham, Nagaraj, "Divisors in residue classes, constructively", Math. Comp. 77 (2008) 531-545 (ePrint 2004/339). Full text NOT READ

- **SNIPPET.** The abstract: they "explicitly construct all such divisors, and incidentally improve this bound to O((alpha-1/4)^{-3/2})". Their theorem as restated by Hales is READ second-hand (above).
- **NOT FOUND.** I have no full-text check of what CHN say about alpha <= 1/4. That includes the "it is not presently known how to place a bound ..." sentence the earlier referee saw in a snippet.

### Cohen, "Diviseurs appartenant à une même classe résiduelle", Sém. Théorie des Nombres de Bordeaux 1982-83

- NOT READ. Content known only through Hales §4.2 (READ second-hand): the 6-divisor family at alpha -> 1/3 from above.

### Later work: arXiv searches (READ: search result lists only)

I searched arxiv.org/search for the phrases "divisors in residue classes", "divisors in a residue class", "divisors in residue" and "divisors in arithmetic progressions", and for the keywords Lenstra+divisors+residue+class and Coppersmith+Howgrave-Graham+Nagaraj.

- The only relevant hit was 2410.05030.
- "Divisors in arithmetic progressions" returned only 2110.00237, "Sums of divisors on arithmetic progressions". From its title it is about divisor sums, not counts in one class; I did not open it.
- Web search turned up the same items plus Hales's OSU thesis, which could not be fetched.
- **No later paper treating alpha <= 1/4 was found. This is a negative search result, not a proof of absence.**

## (b) Constructions: what the elementary family gives, and what is not known

**Trivial family (record's own; elementary, checked by hand).**

- Setup: pairwise coprime factors 1 + iq (i = 1..t) with t! | q, and n = ∏(1 + iq). All 2^t divisors are ≡ 1 mod q, and log q/log n -> 1/t from below.
- Every divisor d ≠ 1 in the class satisfies d >= q + 1. A product of j >= 1 factors each > q has size > q^j. So this family can only use t factors with roughly q^t < n <= q^{1/alpha}, that is t < 1/alpha.
- For fixed alpha it therefore gives about 2^{ceil(1/alpha) - 1}: bounded for each fixed alpha. It grows only as alpha -> 0.
- At alpha = 1/4 - kappa it gives t = 4, i.e. 16. This is consistent with a bounded c(alpha) on [1/4 - kappa, 1/4]. It does not show c(alpha) -> infinity as alpha -> 1/4 from below: the family is constant (16) on [1/5, 1/4).

**Other elementary ideas.** Take n made of many small primes and count divisors equidistributed over phi(q) classes. Then tau(n) = n^{o(1)}, while phi(q) is about n^alpha, so the expected count goes to 0 for fixed alpha. This gives nothing unbounded either. (My own heuristic remark, not a literature statement.)

**Literature.** No construction with c(alpha) -> infinity at any fixed alpha > 0 in an invertible class was found in the one full text read (Hales), nor in any snippet. No proof that c(alpha) is bounded on [1/4 - kappa, 1/4] was found either. The largest literature counts I could verify are both at alpha >= 1/3:

- 6 (Cohen; Hales's examples);
- 7, counting signed divisors (Hales).

## (c) What the record may say

1. "Lenstra (1984) proved LEN(alpha) for every alpha > 1/4, with c(alpha) = O((alpha - 1/4)^{-2}); Coppersmith, Howgrave-Graham and Nagaraj (2008) gave c(alpha) = O((alpha - 1/4)^{-3/2}) constructively, for n > n_0(alpha)." This is supported by Hales 2410.05030v1 §4.1 (READ). The primary sources are SNIPPET only.
2. "Hales (arXiv:2410.05030, §4.2) describes the case alpha near 1/4 as 'much less clear', with significant work only at alpha = 1/3; the open problems he poses concern alpha >= 1/3 and alpha < 1/2." (READ.)
3. "Whether LEN(1/4) or LEN(1/4 - kappa) holds is open as far as we could establish." This is acceptable only with the qualifier: it rests on Hales (full text) plus search snippets, and Lenstra 1984 and CHN 2008 were not read in full. The record should **not** attribute an explicit "open below 1/4" statement to any source: none was found.
4. "No construction of unboundedly many divisors in an invertible class at any fixed exponent is known." Downgrade this to "we found none in the literature we could check (Hales 2024 in full; Lenstra 1984 and CHN 2008 via snippets)". The stronger "is known" cannot be asserted without reading Lenstra 1984 and CHN 2008.
5. "The best known construction gives 16 divisors near 1/4." **Correct this.** It should read: "the elementary family ∏_{i<=4}(1 + iq) gives 16 divisors at exponents tending to 1/4 from below (our computation); the published examples we read (Cohen; Hales) give 6 divisors at alpha > 1/3." The phrase "best known" is not supported by any source read.
6. Small discrepancy to record: Hales's Theorem 4.1 gives "c(1/3) = 12" for Lenstra's combinatorial constant. The record (from snippets) says "at most 11 divisors for q > n^{1/3}". Both can be true, since the combinatorial constant can exceed the divisor bound, but the 11 is unverified here.

## Not resolved

- I did not see the full texts of Lenstra 1984, CHN 2008 / ePrint 2004/339 or Cohen 1982-83; the hosts were blocked.
- So the following remain SNIPPET or NOT FOUND:
  - Lenstra's exact remark on non-invertible classes, and any remark on alpha < 1/4;
  - CHN's remarks on alpha <= 1/4.
- To finish the check, allow eprint.iacr.org or ams.org (the CHN 2008 PDF is open access), or pub.math.leidenuniv.nl (Lenstra), in the session network settings.
