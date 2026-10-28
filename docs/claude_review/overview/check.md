# Adversarial accuracy check of overview.md

I checked the draft against `docs/claude_review/`: the top `README.md`, `two_thirds/proof.md`, `growth/*` and their
referee reports, `threshold/*`, the judges' READMEs of rounds 4–17, `lit18/README.md` and
`local_smoothing/README.md`. I also checked the branch note `inert_prime_cofactor_bound.md` and
`paper/jarnik_sublog.tex`.

Severity:
- **MAJOR**: a misstated theorem or status, a missing quantifier that changes the meaning, an unrefereed
  or evidence claim upgraded, or an overclaim.
- **MINOR**: an imprecise wording, pointer or attribution, a small omission or an inconsistency.

Locations are line numbers in `overview.md`. Paths are relative to `docs/claude_review/`.

**Overall.** Most statements match the referee-corrected versions closely. These all agree with their sources,
including the corrected hypotheses:
- the 2/3 theorem: Thm 7.5 with M≥2 and K_0<13.9092; Cor 8.1; the ℓ_0 values; §10;
- Vojta Thms 5.1–5.2, with R>c_0²C⁴², max(M_0, c_0C²²+1), M_α and ε_α;
- ptolemy Thm D and K·β;
- (L_5) and (L_6) false, and (L_7) open with expected dimension 3;
- the Walsh Thms B (M≥28), B'' and C;
- cotangent Thm A (√24, one prime suffices), Thm B (1.9433, 1.6245 calibration) and the corrected Thm C;
- stability Thm 3.1 and Cor 3.2;
- rank reduction Thm 1.3;
- determinant Thms 2.1 and 4.1, with n≥64, n−1≥300 ln 4r and W≥16 ln n+8log⁺C+14;
- the round-12 T=1 constants;
- H_light (q/3, 0.2349);
- the CIRCLE(q)/(W_q) numbers;
- the round-17 estimates.

The refuted items are handled correctly:
- QS(m≤3) is false, and QS is capped (but see 4);
- L4(K) is refuted;
- SL, SL* and ZE-braid are false;
- (L_6) is false;
- "exactly critical" and "√192" are corrected.

The problems are below.

## Findings

1. **MAJOR. Summary item 1 (lines 29–31).**
   - *Problem.* The draft says "No improvement in the growth rate log R/loglog R is proved, either
     unconditionally or from a standard conjecture." Item 2 and §4 then say that one Vojta instance gives
     the uniform bound. A uniform bound is a growth improvement from a standard conjecture, so as written
     the summary contradicts itself. The record's own sentence (`README.md` growth-rate status;
     `round7/README.md` §0) is always paired with the Vojta exception.
   - *Correct statement.* No growth improvement short of the uniform bound is known, unconditionally or
     from a standard conjecture. No intermediate rate o(log R/loglog R) is derived from any standard
     conjecture. The only standard-conjecture route is Vojta (D=0) on one fixed variety, and it gives the
     uniform bound directly (`README.md` §3; `round7/README.md` §0; `round8/README.md` §2.3(a)).

2. **MAJOR. Theorem 5.3 (lines 680–688).**
   - *Problem.* The draft drops the qualifier that the unconditional part produces fakes that are **not**
     at barrier scale. As written, a reader concludes that gap-admissible fakes at W=O(M log M) exist
     near every nondegenerate tangent-cone point.
   - *Correct statement* (`round12/README.md` §2.4 item 1). For P→∞ there are gap-admissible fakes within
     θe^{−W/4} of u*. "These fakes are **not** barrier-scale: W/(n log n) is unbounded." Only at the
     barrier scale W=O((a/ρ_1+1) r log n), and only under the Diophantine condition (2.2) for **all**
     relevant λ with ‖λ‖_1 ≤ r/ρ_1+1, does the same hold. Round 16 Thm 2.3 restricts (2.2) to light λ.
     The referee also notes that "neither easier nor harder than ZE-TC" is overstated.

3. **MAJOR. §1.4 line 159 ("Exhaustive scan, R²≤10¹² … exact") and §8.1 line 1180 ("scan to 10¹² … exact").**
   - *Problem.* The document defines "exact" as "a finite exact computation, reproduced independently".
     The referee did **not** reproduce the 10¹² scan: "I did not re-run the full 1e12 scan, which would
     take about a day" (`growth/referee_families.md` §5, around line 236). Its independent brute force
     reached only n≤10¹⁰. The referee also found that the "1e8" and "1e10" rows are really n<2²⁶ and
     n<2³³.
   - *Correct statement.* The 10¹² values (k=3..8: 0.2524, 0.5060, 1.9962, 4.0001, 5.2699, 5.8978; at most 3
     points on (1/2)√R arcs and 4 on √R arcs) are the author's computation, so they are **evidence**. They
     were reproduced independently up to 10¹⁰ (`referee_families.md` addendum). Only the explicit clusters
     and families (M(C)≥6, ≥8, ≥10) are proved.

4. **MAJOR. Status of lit18 Props A, B (route table line 774; §6.7 lines 989–991; §7.6 line 1135 and item 5,
   line 1146; §5.4 definition at line 713).**
   - *Problem.* `lit18/README.md` §2.2 heads these propositions "(unrefereed; not leads)": Prop A was recounted
     by a checker, and Prop B was re-derived by the judge. Row 13 of the draft labels them correctly, but
     other places upgrade them:
     - the route table lists QS as "capped" under a definition that requires "a refereed barrier theorem";
     - §6.7 says "Capped" with no label;
     - §7.6 opens with "the refereed barrier theorems say …" and then lists item 5 (Props A, B) among them.
   - *Correct statement.* QS is capped by Props A and B, which are proved in `lit18/check_tools.md` §2.2,
     re-derived by the lit18 judge and **unrefereed** (`lit18/README.md` §2.2). Item 5 of the round-17 list
     should carry that label.
   - The §5.4 definition ("Capped means: proved and refereed that the class contains fakes …") also fails
     for several rows:
     - row 6 is a Θ(log R) cap, not a fake cap;
     - rows 8, 9 and 15 are theorems on varieties, not fake witnesses;
     - row 13 is unrefereed.

     Either widen the definition or label these rows.

5. **MINOR. §5.4 row 12 (line 729).**
   - *Problem.* The draft gives the proved constant as (max(8,2/ε)+o(1)), but the proved statement uses
     4/ε. It also omits the quantifiers 0<ε≤1, 0<C_0≤1 and M=q+1≥M_2(ε,C_0).
   - *Correct statement* (`round8/README.md` §1.1 Thm C and §2.2 row 7). The proved bound is
     W≤(max(8,4/ε)+o(1))M log M + 8M log(1/C_0). The value 2/ε+o(1) comes from "the referee's
     refinement". The result is proved for the max form with Log i in the list only.

6. **MINOR. §5.3 line 690.**
   - *Problem.* "The route-5 hypothesis holds for actual clusters only if it holds on the fake pattern
     boxes" is too strong.
   - *Correct statement* (`round12/README.md` §0 item 5; `round16/README.md` §0). (R_H) and MIN^pair depend
     only on (J(A), V_H', G). So a box-level proof must also hold on the (Mar) fake pattern boxes. There a
     nondegenerate TC_H point that satisfies (2.2) yields gap-admissible fakes in G. §7.2 (lines 1080–1081)
     already states it this way.

7. **MINOR. §5.3 lines 706–709 (genuine labels).**
   - *Problem.* The draft says any lemma about labels alone "(L4, L4*) is satisfied by every fake". L4(K)
     itself is false (Thm 2), and the draft omits the regime qualifier.
   - *Correct statement* (`round10/README.md` §2 item 6). With genuine labels, the hypotheses **and**
     conclusions of L4(K), L4^forced and L4* hold for every fake **whose primes are in the regime**. So
     L4* uses no fake-violated input.

8. **MINOR. §5.4 caveat (lines 741–742).**
   - *Problem.* "Near-extremal actual clusters have all primes ≤ M^{O(1)} …, and no fake of any kind is known
     in that regime" is misleading. Fakes with all primes ≤ M^{O(1)} do exist: for `construct.md`,
     log p=(3+o(1))log M (`round8/README.md` §1.3 item 2).
   - *Correct statement* (`round9/README.md` §0 item 1; `round14/README.md` §0). Every construction of rounds
     5–7 has p_j > X. No fake is known in regime (1.4), where p_j ≤ M^{h_0}, X=3M^A and A is chosen after
     c_0.

9. **MINOR. §5.1 fact (ii) (lines 642–643).**
   - *Problem.* "All known fakes have W=c M log M" is false.
   - *Correct statement.* Some known fakes are larger:
     - the `outside.md` A.3 fakes have log R ≤ (15+o(1))M² log M;
     - the fake Thm 2.1 fakes have W/(n log n) unbounded;
     - the `p3G.md` Prop 6.1 fakes have W ≤ (66+o(1))MX.

     Say instead that the witnesses of the capped rows (1–4, 7, 10–12) are at W=O(M log M).

10. **MINOR. §6.1 lines 803–804 ("round 15 showed that no outcome of the rigidity question moves the growth rate").**
    - *Problem.* This is stated without a label. §6.3 (line 930) labels it correctly.
    - *Correct statement* (`round15/README.md`). The claim comes from `rigSLtheory.md`, which is **unrefereed**:
      its referee was cut off.

11. **MINOR. §6.1 lines 802–803 ("for actual circles its rigidity hypothesis is vacuous").**
    - *Problem.* This misquotes the record.
    - *Correct statement* (`round14/README.md` §1.2 item 2). "For actual circles alone the theorem is
      vacuous: it only routes P3-all* to statements about genuine labels." The corrected transfer theorem
      reads RigS up to (ε,t) + H^gen + L1 ⇒ W ≳ M log²M/loglog M for P3-all* members. Its size claim needs
      the extra hypothesis r ≥ M^γ, which regime (1.4) does not supply.

12. **MINOR. §6.6 lines 972–978 ("Every unconditional method in the record is excluded by a refereed barrier:
    … elimination …, single S-unit equations and abc …, Subspace / Ru–Vojta / Evertse (ptolemy.md §7)").**
    - *Problem.* `ptolemy.md` §7 is the author's informal "why every unconditional tool I know stalls". It is
      not a list of refereed barrier theorems.
    - *Correct statement* (`round17/README.md` §7.3). The judge's claim rests on items 3 and 5 of §7.1:
      - `threshold/L8.md` Thm A, which is refereed as a formal barrier in the precise form of
        `referee_L8.md` §2.3;
      - lit18 Props A and B, which are unrefereed.

      Separately, abc vacuity on squarefree relations is `round8/weakest.md` Props 4.1–4.2 (refereed).

13. **MINOR. §4.1 lines 539–554, heading "Unconditional by-products (proved, refereed)".**
    - *Problem.* The list mixes in two items that are not unconditional by-products:
      - the hyperbola Cor 5.3, which is conditional on V(7,1/4) for the split model and only sketched in
        the note;
      - "Other conjectures".

      Also, "literature status not checked" is inaccurate.
    - *Correct statement.* Move Cor 5.3 out of the list and label it conditional. `referee_conditional.md`
      supplied its details (§1) and ran a WebSearch-only literature check, which found no prior derivation
      from Vojta, Bombieri–Lang or abc (§5; direct fetches blocked). Per the same referee, minor item 4,
      "Bombieri–Lang gives no implication" is supported only as "nothing found" in `conditional.md`. The
      vacuity proof is `round8/weakest.md` Prop 6.2.

14. **MINOR. §5.2 item 2 and §5.4 row 2 (lines 663–664, 719).**
    - *Problem.* "(16+32A) at M^A" is given without its scope limits.
    - *Correct statement* (`round5/README.md`, `lower.md` referee corrections). `lower.md` Lemma 5.2 needs all
      primes of N to exceed X. The class should quantify over all four units. The witnesses meet twisted
      collisions at X=2M, but "the twisted extension at X=M^A, A>1, is open". Lemma S is not needed: item 393
      gives (32+o(1))M log M for b=16.

15. **MINOR. §1.3 item 5 (lines 125–126) and Corollary 3.8 (lines 403–405).**
    - *Problem.* The draft states as fact that Cilleruelo–Granville's open question concerns five points on
      O(R^{2/5}) arcs. The referee could not verify this against CG's papers.
    - *Correct statement* (`referee_cotangent.md` §4 "Literature"). The primary CG papers could not be
      fetched. The five-point reading comes from the project notes (`endpoint_continuation_status.md`) and is
      "consistent with" the four-point families at R^{1/3}. Label the attribution as unverified against the
      source. The equivalence itself (Cor 3.8) is refereed in both directions.

16. **MINOR. §1.3 table, row "(4+ε)… then (1+ε)…" (line 151).**
    - *Problem.* The row's status says "verified by the review" for both constants.
    - *Correct statement.* Only the (1+ε) bound of `inert_prime_cofactor_bound.md` (inequalities (6) and
      (10)) was verified (`README.md` §1). The constant 4 is "the previous prose leading constant", cited in
      that note, and was not audited.

17. **MINOR. §2.1 line 195 ("configurations with M=o(log R/loglog R) are not excluded").**
    - *Problem.* The phrase is copied from `proof.md` §10 but is logically empty.
    - *Correct statement.* The method does not exclude configurations with M ≍ log R/loglog R (up to
      (2/3+ε)log R/loglog R). Its total lower bound is (3/8)M² log M against the cut slack (M/4) log R.

18. **MINOR. §2.2 Step 5 (line 269): "only ℓ<M and p≤M contribute".**
    - *Problem.* This is the uncorrected wording. §8.2 lists the correction, but Step 5 does not use it.
    - *Correct statement* (`two_thirds/proof.md` review box item (v)). Only ℓ ≤ M−2 contribute.

19. **MINOR. §6.2 item 4 (lines 857–859) and §7.2 (line 1078).**
    - *Problem.* The implication MIN^pair ⇒ n=O(loglog R) is stated without its quantifiers.
    - *Correct statement* (`round12/README.md` §2.3 items 5–7 and the table at line 499). MIN^pair must hold
      for the cluster's own T=1 box at m'=⌈4w⌉, for **all regime profiles**, with n ≥ n_1, d ≥ 316 ln 4r and
      W > W_1. MIN* is false on the Sylvester box n=4096, so n_1 > 4096.

20. **MINOR. Theorem 5.4 (lines 694–698) and §6.2 item 6.**
    - *Problem.* The theorem omits the scale at which the test is posed.
    - *Correct statement* (`round17/README.md` §1). The test is posed at the barrier scale W=2qτ, τ ≥ τ_D,
      τ=O(log q), together with the thresholds of barrier Thm 2.3, and with H=H'=e^{m'W/40}. After round 17,
      τ_D is taken with c_0=1/6 (all q) or 0.2349 (regime).

21. **MINOR. §6.2 item 2 (lines 840–841), "Round 11 (refereed) … SL* … false".**
    - *Problem.* The refutation of SL* is attributed to the round-11 referees.
    - *Correct statement* (`round11/README.md` §0 item 1, §5.1). The refutation is the judge's **J1**,
      unrefereed in round 11. `round12/README.md` §2.1 item 10 confirms it ("J1, J2, J3 (R11) are confirmed").
      Cite it that way. §8.2 line 1280 already says "(judge J1)".

22. **MINOR. §5.5 (lines 744–756) omits unclassified items.**
    - *Problem.* `round10/README.md` §7.2 lists more unclassified inputs than the draft includes.
    - *Correct statement.* Add the following:
      - the sharp (X+K) large sieve, which is equivalent to the max-scale cross-gap and is still unclassified;
      - arithmetic fibre inputs (the P3-all class);
      - ZE and its algebraic-point form;
      - prime-weight concentration: whether near-extremal actual clusters must put their primes in one
        window [P, Pe^{1/(4r)}] is open (`round8/README.md` §2.5 item 6).

23. **MINOR. Omission: the `growth/auxiliary.md` barrier.**
    - *Problem.* The draft never mentions it, though the reference list cites Bombieri–Pila and Heath-Brown.
    - *Correct statement* (`growth/README.md`; there is no referee report, so it is **unrefereed**).
      Interpolation determinants on one circle factor as Vandermonde times a symmetric factor with no extra
      arithmetic. So Bombieri–Pila, Heath-Brown and Salberger reduce to pair data. Gap principles,
      larger-sieve variants and Stepanov give nothing beyond O(log R/loglog R). The refereed
      `round10/additive.md` §5 adds that sums, Sidon sets, the σ-cocycle, the F_q polynomial method,
      determinants and Bombieri–Pila all reduce to capped classes. Add a §5.4 row or drop the references.

24. **MINOR. Omission: further refereed lower bounds.**
    - *Correct statement* (`referee_families.md` addendum and §8). The referee found
      M(7.8804) ≥ 9 on the same n=1176852625. The golden family gives M_∞(C) ≥ 7 for C > 5.5365.
      The axis Pell family gives M_∞(C) ≥ 6 for C > 4, not only M(C) ≥ 6.

25. **MINOR. Omission: lw Theorem A(4).**
    - *Problem.* The draft does not state that the conditional LW route defeats the round-5/6 fakes.
    - *Correct statement* (`round8/README.md` §1.1). Under LW, `flipped.md` Thm L makes the genuine versions of
      the `construct`, `lower` and `sharp` (P3-char) fake profiles impossible at Paley orders M=q+1 ≥ M_0
      (ineffective). Only the skew double core escapes single-form LW. This belongs in §3.1, or next to
      rows 10–11.

26. **MINOR. §3.1 lines 337–339 and Theorem 3.3 (line 341).**
    - *Problem.* The comparison with earlier items follows the write-up, which the referee says understates
      them.
    - *Correct statement* (`referee_walsh-extraction.md` §1 and §6). Item 276 already gives an effective
      uniform bound for pure cores of **any** Hadamard type, Paley included, via Matveev, but only for
      C ≤ √2 (C ≤ 2 with a common unit). Item 278 covers pure Walsh for every C. What is new is every C for
      non-Walsh cores, with small explicit constants. For B'', item 277 already allowed about log₂M
      deletions in unit regimes. B'' is "automatic for M ≥ C^{6+o(1)}".

27. **MINOR. Pointers (line 593, and line 1185 of the §8.1 index).**
    - *Problem.* "`referee_ptolemy.md` §3.3" does not resolve relative to `docs/claude_review/`.
    - *Correct statement.* The file is `growth/referee_ptolemy.md` §3.3. Under the document's own path
      convention, the round-local pointers in §3.3 (`stability.md`, `weakest.md`) and §5.2
      (`lower.md`, `sharp.md`) also need their `round8/` and `round5/` prefixes where they are used alone.

## Counts

4 major, 23 minor.
