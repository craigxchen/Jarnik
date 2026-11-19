# Referee report, round 31, Part IV, block "lit": literature and context

Scope: the literature behind Part IV of `astra_free_IV/answer.md` (from line 1011): Langer's inequality and its predecessors, prior art for bounds on sections with prescribed contact, characteristic p, and the arithmetic analogue of (IV.2).

The line-by-line check of (IV.2) is in `ref4_bmy.md` (verdict: CORRECT) and the trace seed is checked in `ref4_trace.md`. I use their results and do not redo them.

My scripts are in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref4_lit/`:
* `weights.py` (sympy): Astra's ledger with a general section weight a.
* `fano.py`: the incidence structure of the trace seed.

**Access caveat, stated once.** arxiv.org and ar5iv.org do not resolve from this session (the host is not in the allowed network domains). I could not read any original paper. Every statement of a cited theorem below comes from WebSearch result snippets (secondary sources) or from memory. I mark each one accordingly. Nothing here was checked against an original text.

## Verdict table

| # | Claim / question | Verdict |
|---|---|---|
| a | Part IV uses Langer, PLMS 86 (2003), Thm 0.1, with weights 2/3 on sections, 1/2 on four fibres, 0 or 2/3 on exceptional curves, and K+B Q-effective. | **CORRECT** within the hypotheses as quoted by secondary sources. Not verified against the original text. |
| b | Is the bound already in the literature? | **Method known, statement apparently new.** The method is the Langer/Bojanowski/Pokora "Hirzebruch-type inequality", transplanted to Urzúa's arrangements of sections. I found no prior bound on the number of sections with prescribed mutual contact. No truncated function-field Vojta or SMT tool gives it. |
| b' | New computation: optimise the section weight a in Astra's ledger. | **COMPUTED.** a = 2/3 is optimal (bound 24). Reduced log BMY (a = 1, the classical Sakai/Miyaoka/Kobayashi form) gives only m <= 32. |
| c | Characteristic p. | The seed is the point–hyperplane design of PG(d-1,2), and for d = 3 it is the **Fano plane** (COMPUTED, d = 3..7). It is the additive-polynomial analogue of the all-F_q-lines arrangements that violate Hirzebruch's inequality. It violates (IV.2) from q = 32 on. Frobenius pull-back alone cannot produce a violation. |
| d | Arithmetic analogue. | It is an **arithmetic orbifold BMY** for a blown-up P^1_Z with weighted sections. This is the Parshin/Moret-Bailly type, open, and the naive version is false. Through K_{M_{0,m}}+∂ it becomes a statement about h_{K+∂} of the point of M_{0,m}. It is not a known special case of Vojta, and it is no weaker than the record's hypothesis. **No unconditional arithmetic case is known.** |

## (a) Langer's theorem and its predecessors

**Langer, Thm 0.1** (as quoted by Pokora, arXiv 1612.05141, and conic-line papers that apply it, e.g. arXiv 2111.12349). Let X be a normal projective surface and D = Σ a_i D_i a Q-divisor with 0 <= a_i <= 1. Assume (X,D) is log canonical and K_X + D is **Q-effective**. Then

  (K_X + D)^2 <= 3 e_orb(X, D).

The secondary sources add that if equality holds, K_X + D is nef. Here e_orb = e_top(X) - Σ a_i e_top(D_i^o) + Σ_x e_orb(x; X, D), with Langer's local orbifold Euler numbers.

In the smooth-SNC case the local values are:
* 1 at a general point;
* 1 - a_i along D_i;
* (1 - a_i)(1 - a_j) at a crossing.

Summing them gives exactly Astra's formula: e(Y) - Σ a_i e(D_i) + Σ_crossings a_i a_j. Ref4_bmy re-derived this independently.

**Predecessors** (from memory and secondary sources, not checked in the originals):
* **Miyaoka 1984**, Math. Ann. 268, "The maximal number of quotient singularities on surfaces with given numerical invariants". This is the orbifold MY for surfaces with quotient singularities, plus a log version with a boundary that Hirzebruch-type arguments use. Its standard form needs K + D nef (or a minimal-model situation).
* **Sakai 1980**, Math. Ann. 254. A reduced SNC boundary D with K + D big gives (K+D)^2 <= 3 e(X \ D).
* **Kobayashi–Nakamura–Sakai 1989**, Proc. Japan Acad. 65A. This covers standard coefficients 1 - 1/n_i, log terminal singularities and K + D nef (and big). Equality characterises ball quotients.
* **Megyesi 1999**, PLMS 78. This covers log canonical pairs with K + D nef.
* **Langer**, J. LMS 2001 and then PLMS 2003, which drops nefness. Two features matter for Part IV:
  * Q-effectivity replaces nefness. This is exactly the hypothesis Astra verifies; ref4_bmy (b) machine-checks the effective representative.
  * Arbitrary real coefficients in [0,1] are allowed, not just standard ones 1 - 1/n.

  **Part IV's weights are not standard orbifold weights. They need Langer's (or Megyesi's lc) generality, not Miyaoka/KNS.**

**Is Part IV's use within the hypotheses?** Yes, subject to the caveat on the original text.
* All coefficients lie in [0,1).
* The boundary is SNC. Ref4_bmy checked this in local coordinates and with its lattice engine.
* Hence (Y,B) is klt, so lc.
* K_Y + B is Q-effective as soon as the section weight a >= 2/m and all other weights are >= 0.

One remark that Astra does not make: the k >= 3 exceptional curves get weight 0, and **no log-crepant adjustment is needed**. Langer's theorem holds for any boundary on the blown-up surface. One does not need to pull back an lc pair from the unblown surface, which would fail at k >= 3 points where (2/3)k > 2.

## (b) Is the bound already known? What would give it more directly?

### (b1) The method is standard; the application is not

* **Hirzebruch-type inequalities via Langer.** Bojanowski (master's thesis, 2003) and Pokora (arXiv 1612.05141; survey arXiv 1808.09167) apply Langer's theorem to (P^2, αΣ lines). Bojanowski's inequality reads

  t_2 + (3/4) t_3 >= d + Σ_{r>=5} (r^2/4 - r) t_r   (d >= 6, t_r = 0 for r > 2d/3).

  The coefficient 3/4 on t_3 and the multiplicity-dependent penalties have the same structure as Astra's per-centre F(k,0,0) = 4 - 8k/9 and φ(k) = C(k,2) + 2k - 9.
* **Tangencies.** For conic and conic–line arrangements (Pokora–Szemberg, arXiv 2002.01760; arXiv 2111.12349; arXiv 2312.12950), tacnodes and A_{2k-1} points are handled with Langer's local orbifold Euler numbers. Astra's "blow up the full chain of e centres, weight a on the chain" is a hand-made SNC substitute for these. It is equivalent in spirit, and ref4_bmy checked it independently.
* **Arrangements of sections.** Urzúa, "Arrangements of rational sections over curves and the varieties they define", Rend. Lincei 22 (2011) 453–486, arXiv 0910.4928. This is the closest predecessor. According to its abstract (search snippets):
  * arrangements of d sections of P^1-bundles over curves correspond to irreducible curves in \bar M_{0,d+1}, via Kapranov;
  * over C the "extended" and "partially extended" arrangements satisfy log Miyaoka–Yau inequalities that are "no longer combinatorial";
  * in any characteristic there are simply connected surfaces of general type built from such arrangements that violate MY-type inequalities.

  As far as the snippets show, Urzúa uses log Chern numbers, that is, a reduced boundary (a = 1). He does not give a bound on the number of sections with prescribed contact. I could not read the paper; this is unverified.

**Searches found no prior statement** bounding the number of sections of a ruled surface, equivalently of points of P^1(C(t)), whose pairwise differences are prescribed-divisible. That is absence of evidence, not proof of novelty. My judgment: the transplant (Langer + Bojanowski-type weights + tangency chains on Urzúa's surfaces) is routine for a specialist. The inequality (IV.2) itself, the normalisation at infinity (IV.3), and the application to the record's cluster model are new as far as I can find.

### (b2) Which weight is best? (new computation, `weights.py`)

**Setup.** I redid Astra's ledger with a general section weight a. Centres with k >= 3 get weight 0 and chain exceptionals get weight a. Two inputs are used: the identity I = (m-1)S/2 and the bound (IV.3). The per-layer contribution is then

  H_a(s) = -(a^2/2) s^2 + (a^2 m/2) s + 4 - a m/2 + 3a^2 m/8 - a^2 m^2/8.

**Maximisation.** The maximum over s is attained at s = m/2 and equals

  4 - (m/4)(2a - 3a^2/2).

This is maximised in a at **a = 2/3**, where it gives 4 - m/6, i.e. **m <= 24**. For a = 2/3 the script reproduces Astra's per-layer identity 4 - m/6 - (2/9)(s - m/2)^2 exactly. Other weights give:

| a | 1/2 | 2/3 | 3/4 | 1 (reduced log BMY: Sakai / Miyaoka / KNS / Urzúa's log Chern numbers) |
|---|---|---|---|---|
| bound | m <= 25 | **m <= 24** | m <= 24 | m <= 32 |

So the classical reduced log-BMY inequality already gives a uniform bound, 32, in this ledger. Fractional weights, which are Langer's contribution, improve it to 24.

**Caveat.** This optimises only within Astra's scheme: weight 0 at k >= 3 centres, the same a on chains and four half-fibres. Other exceptional weights could change the constant. I did not explore them. The constant 4 per layer (the Euler term of a point) is what any such scheme must beat.

### (b3) Tools that do **not** give it directly

* **Truncated function-field Vojta / Nevanlinna with moving targets.** I checked three families:
  * Yamanoi (Acta Math. 2004, small functions; his function-field height inequality h_{K_X(D)} <= N^{(1)}(D) + d + ε h + O(1), quoted from a RIMS Kôkyûroku snippet; X. Chen, Canad. J. Math. 2011, reports that Yamanoi settled Vojta's 1+ε conjecture over function fields, and McQuillan gave an algebro-geometric proof);
  * Ru–Wang and J. T.-Y. Wang (moving-target Cartan/abc);
  * Noguchi.

  All of these have error terms O_D(1) or C(Σ h(a_j)) depending on the targets, or require slowly moving targets. Here the "targets" are the other sections, whose heights equal that of the point. So they give no bound uniform in D. This is the record's standing caution against template-dependent and target-dependent constants. The record had already found that these tools are vacuous (Mason on pairs) or critical (truncated SMT) for this model (`round29/kill_ff.md` items F2 and 189–193; `round30/README.md`).
* **Brownawell–Masser / Voloch / Mason S-unit bounds on M_{0,m}** (cross-ratios are S-units). These bound heights by C(m)(2g-2+|S|). They are not uniform cluster counts, because |S| grows with D and with the cofactor roots.
* **Canonical class inequality** K_{X/B}^2 <= (2g-2)(2b-2+s) (Vojta via Miyaoka–Yau; strict by Tan and K. Liu; higher-dimensional versions by Lu–Tan–Zuo, arXiv 1009.5246). This is truncated: it counts singular fibres s, not contact multiplicities. Its log genus-0 analogue counts distinct collision fibres. I did not derive a bound from it, and I expect none at fixed m, for the same reason the record found Mason vacuous. Part IV's gain comes from weighting each infinitely near contact (each layer) with a per-point Euler term. That is a BMY feature, not a canonical-class-inequality feature.
* **Miyaoka 2008**, Publ. RIMS 44, "The orbibundle Miyaoka–Yau–Sakai inequality and an effective Bogomolov–McQuillan theorem". This is the closest "BMY ⇒ height inequality" precedent: it bounds K·C <= (3/2)(2g-2) + o(g) for curves on surfaces with c_1^2 > c_2, using orbibundles with a parameter α ∈ [0,1]. Part IV is of this genre: BMY on the total space, used as a function-field height-type inequality.

**Conclusion (b).** No literature result I can find gives (IV.2), or any uniform bound, more directly. The most direct literature route is the classical reduced log BMY (Sakai/Kobayashi; Langer covers it). Run through the same accounting, it gives m <= 32 (b2).

## (c) Characteristic p

**Why BMY fails** (standard; secondary knowledge, not checked against originals). In characteristic p, Ω^1 can be Frobenius-destabilised, and Bogomolov's inequality needs strong semistability; Langer, Annals 2004, gives only a char-p version with correction terms. The known counterexamples are:
* Raynaud's surfaces (1978), counterexamples to Kodaira vanishing;
* W. Lang's surfaces with vector fields (1983);
* Szpiro's Frobenius-twisted fibrations, where K^2 grows like p while c_2 does not (Szpiro 1979, on Parshin–Arakelov rigidity);
* Moret-Bailly's non-isotrivial families over P^1 (Astérisque 86, 1981);
* Easton, Proc. AMS 2008;
* density of Chern slopes in (3, ∞) (arXiv 1509.05260).

Combinatorially, the arrangement of all q^2+q+1 F_q-rational lines in P^2(F_q) has t_{q+1} = q^2+q+1 and t_2 = 0. It violates Hirzebruch's inequality for q >= 3. The Fano plane (7 lines, 7 triple points, t_2 = 0) violates Melchior's inequality t_2 >= 3 + Σ(r-3)t_r. Urzúa (above) builds MY-violating surfaces from arrangements of sections in any characteristic.

**The trace seed is an instance of this phenomenon** (`fano.py`, COMPUTED for d = 3..7). With marks a ∈ F_q^* and primes x ∈ F_q^*, the incidence "Tr(ax) = 0" is the point–hyperplane design of PG(d-1,2):
* blocks have size q/2 - 1;
* any two blocks meet in exactly q/4 - 1 points;
* each block together with 0 is an F_2-subspace.

For d = 3 it is the symmetric 2-(7,3,1) design, the **Fano plane**. The sections are built from the additive polynomial Tr(u) = Σ u^{2^j}, which satisfies Tr(u)^2 + Tr(u) = u^q + u (Artin–Schreier). The contact sets are zero sets of additive polynomials. In characteristic 0 these sets would be linear, and the incidence could not occur (cf. Sylvester–Gallai).

**Violation of (IV.2).** With one layer per prime of size q/2 - 1, L = q - 1 and δ = E = 0 (ref4_trace confirms constant cofactors), the left side of (IV.2) is:

| q | 8 | 16 | 32 | 64 | 128 |
|---|---|---|---|---|---|
| LHS of (IV.2) | -87.5 | -97.5 | +170.5 | +1858.5 | +9842.5 |

So **(IV.2) fails in characteristic 2 from q = 32 (m = 31) on**. This is consistent with ref4_bmy's remark that the uniform bound fails over F_{q^2}.

**Frobenius pull-back alone cannot do this.** Base change t → t^p multiplies every root multiplicity, hence every layer count, and L, δ and E by p, so (IV.2) is homogeneous under it. The violation therefore comes from the F_2-linear (additive-polynomial) incidence, i.e. from the "strange arrangement" mechanism, and not from inseparable base change.

I know of no named family in the literature of exactly this form (sections f_a = (t^{q/2} - a^{-q/2} Tr(at))/t). Its incidence geometry is classical (PG(d-1,2), Fano). The two-sided split-norm reading in ref4_trace §4 (hyperbola UV = N) is consistent with this: it is a "split" configuration and not a circle-type one.

## (d) The arithmetic analogue

### (d1) The dictionary

* m primitive columns P_i ∈ Z^2 are m rational points of P^1(Q), i.e. sections C_i of P^1_Z. If the points are distinct, they define one point Φ ∈ M_{0,m}(Q).
* v_p(D_ij) is the local intersection (C_i·C_j)_p. This is exact, as in Part IV, because the columns are primitive.
* A layer h = (p, depth r) with marked group T contributes one unit of intersection, weighted by log p, of the curve Spec Z → \bar M_{0,m} with the boundary divisor D_T. The size bound log|D_ij| <= (log n_i + log n_j)/2 - (log N)/4 + δ is the archimedean place, replacing t = ∞.

### (d2) The quadratic term is K + ∂ of \bar M_{0,m} (new observation; formula from memory)

* By Keel–McKernan / Pandharipande, K_{\bar M_{0,m}} = Σ_T ( |T|(m-|T|)/(m-1) - 2 ) D_T. Hence K + ∂ = Σ_T c_T D_T with c_T = |T|(m-|T|)/(m-1) (= κ_1).
* Since (s - m/2)^2 = m^2/4 - s(m-s), (IV.2) at δ = E = 0 is **equivalent** to

  Σ_h s_h(m - s_h) >= (m^2/4 + 3m/4 - 18) L,

  i.e. (1/(m-1)) Σ_h s_h(m-s_h) >= ((m+4)/4 - 17/(m-1)) L.
* When each layer is a single contact group, the left side is the finite part of deg φ^*(K+∂). The trivial bound s(m-s) <= m^2/4 caps it at ((m+1)/4 + 1/(4(m-1))) L. The two bounds are incompatible exactly when m > 24.
* In words: **BMY plus the archimedean budget force more boundary degree of K + ∂ at the primes of N than the m points can possibly carry.**

### (d3) Number-field form

The arithmetic analogue would read: for primitive P_1..P_m ∈ Z^2 and n_i | N, with the minor size condition at the archimedean place,

  (3/4)(m-24) log N + Σ_{p|N} Σ_{r <= v_p(N)} (s_{p,r} - m/2)^2 log p <= m(m+3) δ + 4 Σ_{i<j} E_ij + O_m(1).

Here δ and E_ij are now in log units, and the O_m(1) must be independent of N and of the configuration. Gaussian-gcd and denominator costs (factors of 2, units) enter δ and E as O(1). At fixed m they are harmless, but they must be kept.

### (d4) What it is

* **It is an arithmetic orbifold BMY**, (K+B)^2 <= 3 e_orb in Arakelov theory, on the arithmetic surface obtained from P^1_Z by blowing up the contact tree, with boundary (2/3)ΣC_i. The arithmetic "e_orb" carries one Euler term log p per blown-up layer; this is the source of the constant 4 per layer.
* This is the **Parshin / Moret-Bailly** type of conjecture: Parshin, Sém. Théor. Nombres Paris 1986–87, Progr. Math. 75 (1988) 299–312; Moret-Bailly, Astérisque 183 (1990).
  * Per search snippets, versions of it are related or equivalent to effective Mordell, Szpiro and abc.
  * The **most naive arithmetic BMY is false**: genus-2 counterexamples are reported in Compositio, "A discriminant and an upper bound for ω^2 for hyperelliptic arithmetic surfaces", 1999.
  * Miyaoka's 1988 announcement of a proof of an arithmetic Bogomolov–Miyaoka–Yau inequality, with Fermat as an application, was withdrawn after a gap was found. This is common knowledge; the searches only surfaced his 1988 talk titles. No arithmetic BMY of any kind is proved.
* **It is not a known special case of Vojta's conjecture.**
  * The closest Vojta statement is the Main Conjecture for (\bar M_{0,m}, ∂), where K + ∂ is ample for m >= 5 (Keel–McKernan), or its truncated "general abc" form h_{K+∂}(Φ) <= N^{(1)}_∂(Φ) + ε h + O(1) off Z.
  * Vojta's statements have an exceptional set Z_ε and constants depending on ε. With S fixed, the constants may also depend on S, which is the record's caution against assuming uniform exceptional sets.
  * The analogue of (IV.2) asks for **no exceptional set and an explicit constant**. Part IV's function-field proof delivers exactly that, but only because BMY is a global inequality on a single surface.
  * In the function field, Parshin-type arguments derive Vojta-type height inequalities from BMY (Miyaoka 2008; the canonical class inequality). The arithmetic converse is unknown.

### (d5) Comparison with the record's hypothesis (expert note §2: Vojta with D = 0 on \bar M_{0,8} ⇒ (L_8) ⇒ uniform bound)

* **That hypothesis** is one fixed variety and has no counting function (D = 0, so S is irrelevant). It needs the exceptional set Z and the balanced-selection steps.
* **The analogue of (IV.2)** is uniform over varying arithmetic surfaces, needs no selection, and gives m <= 24 for integral fakes directly (hence for genuine clusters). It is formally incomparable with Vojta-D=0-on-M_{0,8}. In strength it is at least of abc/Szpiro type: arithmetic BMY is the parent conjecture of those.
* **So it is not a weaker or more accessible hypothesis than the one the record already uses.** It is a different, and arguably stronger, conditional route.

### (d6) Unconditional arithmetic cases

**None known.**
* The ingredients that do transfer are proved:
  * the identity 2C_i·C_j = C_i^2 + C_j^2 and the elementary-transformation normalisation (IV.3), whose arithmetic analogue is the arithmetic Hodge index theorem (Faltings–Hriljac) together with lattice reduction or choice of integral model (the record's star matrix and Minkowski);
  * the arithmetic Noether formula;
  * arithmetic Bogomolov–Gieseker for hermitian bundles (Moriwaki, Amer. J. Math. 1995; Soulé).
* The missing ingredient is the c_2/Euler-number side applied to a cotangent object. Spec Z has no Ω^1, which is the record's "no derivation" point (`round29/moon.md`). Moriwaki's Bogomolov–Gieseker therefore has nothing to act on.
* Unconditional abc-type bounds (Stewart–Yu) are exponentially too weak.
* For m <= 24 the inequality has no content.

## Novelty and significance (blunt)

* **Novelty.**
  * The method is standard: Langer's theorem with Bojanowski-type fractional weights, tangency chains as in conic-arrangement papers, on Urzúa's surfaces of sections.
  * The specific inequality (IV.2), the infinity normalisation and its use on the record's model are new to the record. The record's only BMY computation is `round4/L7.md`, a seven-curve "no obstruction" check. I found them nowhere in the literature.
  * The trace seed's incidence is classical (PG(d-1,2), Fano), and its BMY failure is of the known "F_q-linear arrangement" type. Its realisation by polynomial sections with constant cofactors appears new.
* **Significance.**
  * For the function-field model it is decisive: ref4_bmy derives M(D,0) <= 24 and sup_D M(D,c) < ∞, both listed OPEN in rounds 29–30. The record predicted that this would need "truncated function-field Vojta". It does not: BMY suffices, and truncated tools would not have done it.
  * For the integer problem it gives nothing. Its arithmetic analogue is an arithmetic orbifold BMY, a conjecture at least of abc/Szpiro strength, whose naive form is false. No unconditional case is known. It is not easier than the record's existing Vojta-on-\bar M_{0,8} hypothesis.
  * The F1/F2 picture is unchanged. The one useful new reading is (d2): the obstruction says that the finite part of h_{K+∂} on \bar M_{0,m} would have to exceed its trivial maximum. This connects Part IV to the same divisor K + ∂ that drives the record's Thm D, and an arithmetic attack would have to bound exactly that quantity.

## Sources (via search snippets; originals not read)
* [Pokora, arXiv 1612.05141](https://arxiv.org/pdf/1612.05141); [Pokora survey, arXiv 1808.09167](https://arxiv.org/pdf/1808.09167); [Langer, arXiv math/0012180](https://arxiv.org/abs/math/0012180)
* [Conic-line arrangements, arXiv 2002.01760](https://arxiv.org/pdf/2002.01760); [arXiv 2111.12349](https://arxiv.org/pdf/2111.12349); [arXiv 2312.12950](https://arxiv.org/html/2312.12950v2)
* [Urzúa, Rend. Lincei 2011](https://ems.press/journals/rlm/articles/5209), [arXiv 0910.4928](https://arxiv.org/pdf/0910.4928); [Urzúa, arXiv 0711.0765](https://arxiv.org/pdf/0711.0765); [Chern slopes in char p, arXiv 1509.05260](https://ar5iv.arxiv.org/html/1509.05260)
* [X. Chen, On Vojta's 1+ε conjecture, arXiv 0705.1727](https://arxiv.org/pdf/0705.1727); [Yamanoi, Nagoya Math. J.](https://www.projecteuclid.org/euclid.nmj/1114631982); [RIMS Kôkyûroku 1319-3](https://www.kurims.kyoto-u.ac.jp/~kyodo/kokyuroku/contents/pdf/1319-3.pdf)
* [Lu–Tan–Zuo, arXiv 1009.5246](https://arxiv.org/pdf/1009.5246); [Miyaoka, Publ. RIMS 2008](https://ems.press/content/serial-article-files/41034?nt=1)
* [Kollár, arXiv math/0602562](https://arxiv.org/pdf/math/0602562); [BMY, Wikipedia](https://www.wikipedia.org/wiki/Miyaoka%E2%80%93Yau_inequality); [Moret-Bailly, Astérisque 183](https://www.numdam.org/book-part/AST_1990__183__37_0/); [hyperelliptic ω^2, Compositio](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/09FD9ED09D419B1562A733AE76E809B2/S0010437X99000032a.pdf/a-discriminant-and-an-upper-bound-for-w2-for-hyperelliptic-arithmetic-surfaces.pdf)
