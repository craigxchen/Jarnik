# Referee report, block "meaning": novelty, literature and consequences of Part VIII

Material refereed: Part VIII of `round31/astra_lemma_VIII/answer.md` (lines 1742-1929), `final8.txt`, and the
coordinator's question `reply8.txt`. I read Astra's `viii_check.py` but did not import or run it. My scripts and
their outputs are in `round31/ref8L_meaning/`: `m1.py` → `out_m1.txt`, `m2.py` → `out_m2.txt`. Record paths are
relative to `docs/claude_review/`.

## 0. Results

| # | Statement | Label |
|---|---|---|
| R1 | "(L_8) ⇒ Rudin" is **new to the record**. The record never derived Rudin, CG Conjecture 9 or a two-residue statement from any hypothesis. It did say the opposite in places (R2). | Proved (by reading the record) |
| R2 | The record says in three places that Cilleruelo–Granville (CG) prove the direction "arcs ⇒ Rudin-type statements". That is **inaccurate** for the text I read (arXiv math/0608109v1). CG's lattice-point route to Rudin is Theorem 12, from **Conjecture 12**, a congruence condition on one coordinate of lattice points on one circle. It is not an arc statement. CG relate their arc conjectures 13–15 only to one another and to the L^4 Conjecture 19 (Theorem 16). | Proved (by reading) |
| R3 | Literature. I found no source deriving Rudin's conjecture σ(k)=O(√k), CG Conjecture 9 or the two-residue statement from Vojta's conjecture. The known conditional results are weaker: Bombieri–Lang with Caporaso–Harris–Mazur gives σ(k)≪k^{4/5} through sumsets and σ(k)≪k^{5/7} directly (CG, read). A uniformity-conjecture bound ≪|P|^{2/3} (Shkredov–Solymosi) is known to me only from a search snippet. So "V_8 ⇒ Rudin" appears to be new, subject to the caveats in §3. | Literature (bounded search) |
| R4 | Astra's citation should be corrected. In arXiv v1, "Conjecture 9 ⇒ Conjecture 1 (Rudin)" is **Theorem 11**. Theorem 12 is "Conjecture 12 ⇒ Conjecture 1". I could not open the "final text" Astra cites. | Correction |
| R5 | Astra describes the two-class statement as going "a fixed power beyond the critical version of their two-residue problem". That description is accurate. It is the d=2 case of CG Conjecture 10 at the exponent 1/2+ε for one fixed ε>0. It lies strictly between CG's "critical" form (exponent 1/d = 1/2), which it implies, and their full Conjecture 10 (exponent 1−ε for every ε), which implies it. The step "two-class ⇒ Conjecture 9" is CG's own remark ("Conjecture 9 would follow easily from Theorem 9 if we could get the exponent 1/2"), completed by the nonunit reduction of their first proof of Theorem 10. Astra's genuinely new step is "(L_8) on square-block data ⇒ two-class at q^{1/2+ε}". | Proved (re-derived) |
| R6 | The two-residue problem and the arc problem have **identical Liouville exponents**. CG's Theorem 10 forces k+1 roots to span more than b^{1/2−1/(2ℓ)}, with ℓ the largest odd integer ≤ k+1. Cilleruelo–Córdoba forces k+1 lattice points to span an arc of length at least about R^{1/2−1/(4⌊k/2⌋+2)}. The exponents agree for every k. The CG proof is the same cut and Vandermonde count as the record's Theorem 1.1. | Computed (k ≤ 59, exact) and Proved |
| R7 | **Logical map** (§4). Proved arrows: V_8 ⇒ (L_8) ⇒ U; (L_8) ⇒ (L_8)^sq ⇒ TR(1/2+ε) ⇒ SQ ⇒ RUD; (L_8) ⇒ LEN(1/4−κ). Every other arrow asked about is open: U ⇒ RUD, RUD ⇒ U, TR ⇒ U, SQ ⇒ U, QB ⇒ RUD, (L_8^pc) ⇒ RUD. | Proved / Open as marked |
| R8 | The uniform lattice bound is **not known** to follow from the two-residue statement. CG do not claim it, and I found no reduction. | Open |
| R9 | Astra is right that Rudin does not obviously follow from the Gaussian short-arc bound alone. The natural route is CG's sum-of-two-squares device. That device puts both coordinates in **one** interval, and the record proved this one-interval case unconditional (round 4, Prop. 2.3(2)). So a critical-arc bound adds nothing to that route. This is not a proof that U fails to imply Rudin. | Proved (for the route) / Open (in general) |
| R10 | Part VIII is no evidence against (L_8). Every consequence is believed true and none is known to fail. My data: at most 3 roots of x²≡a mod q lie in any interval of length √q, for q ≤ 3000 and all a. It is strong evidence that Question A is **at least as hard as Rudin's conjecture** and as a divisor bound below the Coppersmith/Lenstra 1/4 barrier. | Assessment (Heuristic) + Computed |
| R11 | The standing sentence "no intermediate growth bound is known from a weaker standard conjecture" **stands**. Part VIII refutes one candidate (VIII.4). The hypotheses of VIII.5–6, which would give growth, are not standard and are not compared with Vojta. A new, separate sentence is needed on the hardness side (§7). | Assessment |

I also re-derived VIII.1 (the probabilistic extraction) and the identity (VIII.3) by hand, because the meaning of the
implication depends on them. My own exact check of the square-block identity passed 79,632 pair checks with no
failure, on random two-class CRT eight-tuples with coalescing residues (`out_m1.txt`). The blocks `rudin` and
`extract` referee these steps in detail. I treat them as correct; nothing in my reading contradicts them.

## 1. Definitions (nodes of the map)

* **V_8.** Vojta's Main Conjecture on X = M̄_{0,8} over Q with D = 0, A = K+B, and ε = 1/43 (`growth/ptolemy.md` §6).
  Off a proper closed Z_ε it says h_K ≤ ε h_A + C_ε.
* **(L_8).** As in `expert_note.md` §2: there are c, η_0 and Z such that every anchored datum with the chain condition,
  η-balanced on the 119 visible blocks and with Φ ∉ Z, has H_Σ ≥ cw − O(1).
* **(L_8)^sq.** (L_8) asserted only for data whose assigned blocks n_T are pairwise coprime perfect squares, with
  arbitrary nonzero integer cofactors t_ij. The assigned block need not equal the contact integer d_T(Φ).
* **(L_8^pc).** The lattice form of `expert_note.md` §2: pairwise coprime e_T, pair cofactors of size e^{w(1±η)}.
* **TR(θ).** There is B such that, for every q, if one prescribes at most two residues mod p^e for each p^e ∥ q,
  then every real interval of length q^θ contains at most B integers lying in the prescribed classes at every p.
  TR(1/2+ε) for one fixed ε > 0 is Astra's VIII.2. CG Conjecture 10 with d = 2 is TR(1−ε) for every ε > 0. CG's
  "critical" belief for d = 2 is TR(1/2).
* **SQ.** CG Conjecture 9, for all a, units and nonunits alike: an interval of length √q contains O(1) roots of
  x² ≡ a (mod q), uniformly in a and q.
* **RUD.** Rudin's conjecture, CG Conjecture 1: a nonconstant k-term integer progression contains O(√k) squares.
* **CG12.** CG Conjecture 12: there are δ > 0 and m such that if a_i² + b_i² = n and a_i² ≡ a_1² (mod q) for
  i = 1, …, m, then q = O(n^{1−δ}).
* **LEN(θ).** There is K such that an invertible residue class mod q ≥ n^θ contains at most K divisors of n.
  Lenstra proved LEN(α) for every α > 1/4, with a constant depending on α.
* **U.** M(C) < ∞ for every C (the critical uniform Cilleruelo–Granville conjecture). **CG14** is their
  R^{1−ε} arc conjecture.
* **QB.** Question B, (L_5°) on the norm-one torus (`expert_note.md` §3). QB ⇒ U is in the record.

## 2. Novelty against the record

**What the record said** (I read every file named in the task and grepped for Rudin, two residue, arithmetic
progression, Lenstra and divisors).

1. `round4/hardness.md` §3, Theorem 2 and Prop. 2.3. Every consequence of U that passes through additive energy or
   sumsets is unconditional at exponent 1/2. Corollary 2.2 gives at most √(6K) squares in {a+qn : n < K} when
   a ≥ (qK/2)^{4/3}. U is equivalent to a two-interval L^∞ statement whose one-interval case is trivial. The author
   wrote "This closes the Cilleruelo–Granville link to Rudin's conjecture as a hardness route at exponent 1/2."
   `round4/referee_hardness.md` M3 struck this: the author had not read CG, and Theorem 2 covers only
   second-moment consequences. The corrected statement is "every consequence of U for squares from square-root
   intervals that factors through the second moment is unconditional."
2. `round4/hardness.md` Prop. 6.1 and Theorem 4 (with the referee's quantifier corrections). Residue classes at or
   above the Gaussian Lenstra exponent 1/4 hold at most one cluster point. The divisor problems of Erdős–Rosenfeld
   and Ruzsa are the **split twin** of the circle problem, on the split form of the same toric variety. No
   Q-homomorphism of the tori exists in either direction. `round29/clothes.md` and `round29/kill_clothes.md`
   describe the residue twin as "Lenstra at exactly α = 1/4". No implication between the twin and the arc problem
   is known.
3. `overview/overview.md` line 180, `lit18/direct.md` §2.10 ("CG prove only the direction arcs ⇒ squares-in-AP and
   Rudin-type statements. An improvement of BZ would not transfer.") and `local_smoothing/chains.md` row F ("CG:
   arc / short-interval boundedness conjectures **imply** Rudin-type conjectures. The direction is arcs ⇒
   squares."). `lit18/gh/cg_index.md` lists CG's Theorem 12 correctly ("Conjecture 12 … implies Rudin's
   Conjecture 1") but does not mention Conjecture 9, Conjecture 10 or Theorem 11.
4. `round31/ref5_lit.md` (b): CG connect squares in progressions and Rudin to Bombieri–Lang (from the abstract
   only). `growth/referee_conditional.md` and `growth/referee_ptolemy.md` §4: CG use Bombieri–Lang and abc only for
   sumsets and squares in progressions.
5. `growth/ptolemy.md` §6.2 Remark 3 and §7: the proof of Theorem D "never uses the Gaussian norm structure beyond
   the extraction", and (L_8) "involves no Gaussian integers and no circle."

**Verdict on novelty.** No file in the record derives Rudin's conjecture, CG Conjecture 9, a two-residue statement
or a Lenstra-type statement below 1/4 from (L_8), V_8 or any other hypothesis. The record only recorded barriers
(items 1–2) and a direction of implication (item 3). So **"(L_8) ⇒ Rudin" is new to the record**, and so is
"(L_8) ⇒ LEN(1/4−κ)". Both concretely realise item 5: the hypothesis is blind to the norm form, so it also governs
the split (rational) configurations. Rudin and Lenstra live on that side (item 2).

**Correction to the record (R2).** For the text I read, item 3 misdescribes CG. The facts in arXiv math/0608109v1
are these.
* §4: Theorem 11 is "Conjecture 9 ⇒ Conjecture 1 (Rudin)". Conjecture 9 is a quadratic-congruence statement, not
  a lattice-point statement.
* §5: Theorem 12 is "Conjecture 12 ⇒ Conjecture 1". Conjecture 12 is about m representations n = a_i² + b_i² with
  all a_i² congruent mod q. Its content is a congruence on one coordinate, at moduli q ≫ n^{1−δ} ≥ (√n)^{2−2δ}. It
  is again a two-residue-type statement at exponent 1/2+δ′. No arc appears in it.
* The arc conjectures 13, 14 and 15 are shown equivalent to one another (p. 11). Theorem 16 shows Conjecture 13 ⇒
  Conjecture 19, an L^4 bound for squares in [N², (N+N^α)²], which CG call "trivial for α ≤ 1/2".
* Neither flowchart (pp. 10 and 15) has an arrow from an arc conjecture to Rudin.

So "CG prove arcs ⇒ Rudin" should read: "CG prove that two congruence conjectures imply Rudin: Conjecture 9 (square
roots mod q in intervals of length √q; Theorem 11) and Conjecture 12 (lattice points on one circle with congruent
first coordinates; Theorem 12). Their arc conjectures 13–15 are linked to an L^4 conjecture for squares
(Theorem 16), not to Rudin." The sentence "an improvement of BZ would not transfer" remains correct. I found no
implication from Rudin to arcs.

## 3. Literature

**What I read in full text.** Only the extraction `ref6L_pell/dl/cg07.txt` of Cilleruelo–Granville, arXiv
math/0608109v1 (3 Aug 2006). I read pp. 1–12 in place (§§1–5, including the proofs of Theorems 2, 8, 9, 10, 11, 12
and 13 and the flowchart), §6 up to Theorem 16 and the second flowchart, and §9 (abc) up to Lemma 1. Direct fetches
of arxiv.org failed: the host does not resolve in this session, so arxiv.org would have to be allowed in the
session's network settings. I did not open the published CRM version or the "Conjavi.pdf" final text.

**From CG (read).**
* Bombieri–Lang with Caporaso–Harris–Mazur gives a uniform bound B on rational points m with f(m) a square, for
  deg f ∈ {5, 6}. This yields Σ r²_{E+E} ≪ |E|^{11/4} (Theorem 2), hence σ(k) ≪ k^{4/5} through Theorems 3 and 4.
  A direct counting argument yields σ(k) ≪ k^{5/7}. CG remark that both are superseded unconditionally by
  Bombieri–Granville–Pintz (k^{2/3+o(1)}) and Bombieri–Zannier (k^{3/5+o(1)}).
* §9: the abc conjecture in Baker's strong form (9.1) bounds the size of progressions with many squares; it does
  not give Rudin.
* §4: Theorem 9 is unconditional, the TR exponent 1/d − ε_d(k). Theorem 10 and Konyagin–Steger give O(log q)
  roots in an interval of length √q. CG "strongly believe" Conjecture 10 at exponent 1/d. Poonen asked whether
  Conjecture 10 holds at exponent 1/2 for d = 4, because that would imply uniform boundedness of rational
  preperiodic points of quadratic polynomials.
* CG nowhere derive Rudin or Conjecture 9 from Vojta, Bombieri–Lang or abc.

**From search snippets only (not read).**
* Shkredov–Solymosi, "The uniformity conjecture in additive combinatorics" (arXiv 2005.11559; SIAM J. Discrete
  Math. 2021). Under the uniformity conjecture, a set of squares meets a progression P in ≪ |P|^{2/3} points.
  The snippet's claim that this improves Bombieri–Zannier's 3/5 is inconsistent, since 2/3 > 3/5; I could not
  resolve this without the paper.
* González-Jiménez–Xarles (arXiv 1301.5122; LMS J. Comput. Math.): Rudin's strong form is verified for
  6 ≤ N ≤ 52, with extremal progression 24n+1.
* Xarles (arXiv 0909.1642): an unconditional bound, depending on the degree, for squares in arithmetic progression
  over number fields; six for quadratic fields.
* Hales, "Divisors in residue classes revisited" (arXiv 2410.05030), §4: Lenstra's c(α) for α > 1/4;
  Coppersmith–Howgrave-Graham–Nagaraj O((α−1/4)^{−3/2}).
* Caporaso–Harris–Mazur (JAMS 10, 1997): Bombieri–Lang implies uniform bounds B(d, g) on rational points of
  curves of genus g ≥ 2.

**Not read at all.** Bombieri–Granville–Pintz (Duke 1992) and Bombieri–Zannier (2002), apart from CG's account and
the snippets.

**Conclusion (with caution).** My searches covered Vojta with Rudin, Vojta with squares in arithmetic progression,
Bombieri–Lang and CHM with Rudin, and square roots mod q in √q intervals. None returned a derivation of RUD, SQ or
TR(θ ≥ 1/2) from Vojta's conjecture, from Bombieri–Lang or from abc. The strongest conditional bounds I found are
k^{5/7} (Bombieri–Lang, CG) and, per a snippet, k^{2/3} (uniformity, Shkredov–Solymosi). Within the limits of a
bounded web search with arxiv.org unreachable, **"Vojta's Main Conjecture for the single variety M̄_{0,8} (D = 0,
ε = 1/43) implies Rudin's conjecture" appears to be new**. It should be called "new to us; we found no earlier
derivation", not "new".

There is also a qualitative reason it is unlikely to be classical. The Bombieri–Lang and CHM route bounds rational
points on the curves y² = ∏_{j≤5}(x + n_j) for each fixed 5-subset of the progression. Each curve then contributes
O(1) points, and counting 5-subsets cannot produce √k. Part VIII's route instead uses a single 5-fold and a
GCD-type (D = 0) inequality whose heights add contacts over all primes of q. That is closer to the
Silverman / McKinnon–Roth blow-up mechanism than to uniformity on curves.

## 4. Logical map

```
                      V_8  (Vojta, D=0, on Mbar_{0,8}, eps=1/43)
                       |  Thm D (record, conditional, refereed)
                       v
   QB = (L_5°) ===> U <=== (L_8) ===> LEN(1/4-kappa)            [VIII.3]
   (record)        ^  (record)  |
                   |            | restriction (trivial)
               ?? (open)        v
                   |        (L_8)^sq
                   |            |  VIII.1 + VIII.2 (Astra; refereed)
                   |            v
   CG Conj 10 (d=2, all eps) ==> TR(1/2+eps) ==> TR(1/2) ("critical")
                                 |
                                 |  CG remark + nonunit reduction (CG Thm 10, 1st proof; Astra)
                                 v
                    CG12 -----> SQ = CG Conj 9
                      \          |  CG Thm 11 (v1)
          CG Thm 12    \         v
                        `-----> RUD (Rudin) ---> sigma(k) = O(k^{1/2+eps})
   Bombieri-Lang (+CHM) ---> sigma(k) << k^{5/7}        (CG, read)
```

| Arrow | Status | By whom / where |
|---|---|---|
| V_8 ⇒ (L_8) | proved, conditional | record, `growth/ptolemy.md` Thm D; c = 6/69, η_0 = 1/381 |
| (L_8) ⇒ U | proved | record, `ptolemy.md` §7 with Steps 0–5 of §6.2 |
| QB ⇒ U | proved | record, `round20/gauss.md` Thm 3.4′ (expert note §3) |
| (L_8) ⇒ (L_8)^sq | trivial | restriction |
| (L_8)^sq ⇒ TR(1/2+ε) for some absolute ε > 0 | proved (conditional theorem) | Astra VIII.1–2; blocks `extract` and `rudin`; my re-derivation and check (79,632 pairs) |
| TR(1/2+ε) ⇒ TR(1/2) | trivial | |
| CG Conj 10 (d = 2) ⇒ TR(1/2+ε) | trivial (special case) | |
| TR(1/2+ε) ⇒ SQ (all a) | proved | Astra VIII.2 (nonunit reduction); in substance CG's §4 remark and the first proof of their Thm 10 |
| SQ ⇒ RUD | proved | CG Theorem 11 (arXiv v1) |
| CG12 ⇒ RUD | proved | CG Theorem 12 (arXiv v1) |
| TR(1/2+ε) ⇒ CG12 | not checked | plausible via the coordinates a_i ≤ √n with q ≥ n^{1−δ}. Nonunit a_1 needs care because the reduction only reaches length √q; Gap |
| (L_8) ⇒ LEN(1/4−κ) | proved (conditional theorem) | Astra VIII.3; block `lenstra`; κ = 10^{−17} under V_8 |
| V_8 ⇒ RUD, V_8 ⇒ LEN(1/4−10^{−17}) | proved (composition) | new to the record; not found in the literature |
| Bombieri–Lang ⇒ σ(k) ≪ k^{5/7} | proved | CG p. 4 |
| U ⇒ RUD | **open** | no route known (§5) |
| RUD ⇒ U, SQ ⇒ U, TR(θ) ⇒ U for any θ < 1 | **open** | CG claim none; I found no reduction (R8) |
| QB ⇒ RUD | **open** | QB lives on the anisotropic torus; RUD lives on the split side |
| (L_8^pc) lattice form ⇒ RUD | **not shown** | VIII.2 needs nonzero cofactors up to e^{(a/4+λ)W}, not necessarily pairwise coprime; the lattice form does not obviously cover them |
| RUD ⇒ (L_8)^sq, SQ ⇒ (L_8)^sq | **open** | (L_8)^sq demands a gain c·w over Liouville on balanced 8-tuples; SQ gives only O(1) at √q |

**Does U follow from the two-residue statement? (R8).** CG do not say so; neither flowchart links §4 to §5's arc
conjectures. In the obvious dictionary a primitive lattice point (x, y) of norm N corresponds to a root t of
t² ≡ −1 (mod N), a two-residue set, through x ≡ ty. Nearby angles do not correspond to nearby t; the map passes
through lattice reduction. So no interval statement about Ω(N) transfers to arcs. The two problems share their
Liouville mechanism exactly (R6), but structurally they are twins on different tori: the split and anisotropic
forms of `round4/hardness.md` Theorem 4. I know of no transfer in either direction. **Open.**

## 5. Rudin and the Gaussian short-arc bound (R9)

Astra writes that she "could not deduce Rudin from the Gaussian short-arc bound alone". That is the right
statement, and the record already contains the reason the obvious route fails.

* **CG's lattice route (Theorem 12) uses a congruence, not an arc.** Squares x_1² < … < x_ℓ² from a + b, …, a + kb
  give points v = (x_i, x_j) on circles n = x_i² + x_j². Both coordinates lie in one interval
  I = [x_1, x_ℓ] of length L ≤ (k−1)b/(x_1 + x_ℓ).
* **Both coordinates in one interval is the trivial case.** Points with both coordinates in a single interval of
  length c·n^{1/4} number at most 2(c²/√2 + 2), unconditionally (`round4/hardness.md` Prop. 2.3(2), refereed).
  That is the one-interval case of the box form of U.
* **So U adds nothing to this route.** Its only content beyond what is proved is the two-separated-interval case,
  which this route never produces.
* **The known proof needs a congruence.** CG's proof of Theorem 12 needs the congruence information r_j, s_j ≡ ±w
  mod p^e ∥ b at the primes of b. That information is encoded in Conjecture 12, which is a two-residue statement,
  not an arc bound.
* **Scale check.** The arc of the points v has length about kb/√n. It is a critical arc (≤ C n^{1/4}) only when
  b ≲ C n^{3/4}/k, i.e. in the far regime a ≳ (bk)^{4/3}. There Rudin is already unconditional by energy
  (`round4/hardness.md` Cor. 2.2).

**Label.** Proved: U gives nothing for Rudin through the sum-of-two-squares route. Open: whether U implies Rudin by
any other route. The hardness that Part VIII exhibits belongs to the route through (L_8), not to U.

## 6. What this means for (L_8) and for Question A (R10)

**Credibility.** Part VIII does not make (L_8) less credible. The consequences are Rudin's conjecture, CG
Conjecture 9 for all a, TR(1/2+ε) and LEN(1/4−κ). All four are believed true.
* CG "strongly believe" TR(1/2).
* Rudin's strong form is verified for N ≤ 52 (snippet).
* My exhaustive count finds at most 3 roots of x² ≡ a (mod q) in any half-open interval of length √q, for all
  2 ≤ q ≤ 3000 and all a (`out_m2.txt`, Computed).
* The `lenstra` block found nothing beyond the 2^k pattern below 1/4.
* No construction in the literature I saw beats TR at any exponent below 1. CG's Theorem 8(3) puts 2^k elements in
  an interval of length εq, which is consistent.

A refutation of any of these would refute (L_8) and V_8. That is a new and **cheap test** of (L_8): search for
eight-element two-class configurations of length q^{1/2+10^{−16}}. Under V_8 such a test is hopeless in practice,
since the margins are 10^{−16}.

**Feasibility. This is the real content.**
* **Rudin.** A proof of Question A, even of its restriction (L_8)^sq, would prove Rudin's conjecture. Rudin's
  conjecture has been open since 1960. The best exponent, 3/5 + o(1) (Bombieri–Zannier), comes from Faltings'
  theorem and Diophantine approximation on genus-5 curves. Even Bombieri–Lang with CHM gives only k^{5/7}.
* **Lenstra.** Full (L_8) would also give a uniform divisor-in-residue-class bound below Lenstra's exponent 1/4.
  1/4 is exactly where the Coppersmith lattice method, and every combinatorial argument of Lenstra type, stops.
* **The new difficulty is not specific to the circle.** It comes from the fact, noted in `ptolemy.md` Remark 3,
  that (L_8) is blind to the norm form. It therefore governs the split (rational) twins of the circle problem:
  two residues per prime power, divisors in classes. The record had already identified those twins as recognised
  hard problems at the same critical exponent (§2 items 2–3). Part VIII turns that analogy into implications.

**Blunt assessment for an expert reader.** A number theorist would conclude that Question A is not a reasonable
target on its own. It is a reformulation that contains Rudin's conjecture and a below-1/4 Lenstra bound. Any proof
of it would be a landmark well beyond the lattice-point problem. Question A remains valuable as a precise,
Vojta-predicted statement and as a test object (refuting it refutes V_8). The route U ⇐ (L_8) is now known to be
strictly harder than anything the uniform bound is known to imply. Anyone attacking the uniform bound should prefer
inputs that use the Gaussian norm structure, as Question B does. Question B is not known to imply Rudin, but that
absence is no evidence of easiness: QB is also a Vojta instance with no unconditional gain over Liouville
(`expert_note.md` §3, §5 item 12).

## 7. Proposed wording for the expert note

**§2, new paragraph after "Theorems around (L_8)" (proposed).**

> **Strength of $(L_8)$** (round 31, Part VIII; refereed). $(L_8)$ implies Rudin's conjecture: every nonconstant
> $k$-term integer progression contains $O(\sqrt k)$ squares. The proof uses $(L_8)$ only on data whose assigned
> blocks are pairwise coprime perfect squares, with arbitrary nonzero cofactors. It runs through three statements.
> (i) If one prescribes at most two residues modulo each exact prime power $p^e\parallel q$, then any interval of
> length $q^{1/2+\epsilon}$ contains $O(1)$ integers in the prescribed classes, for one absolute $\epsilon>0$. This
> is the case $d=2$ of Cilleruelo–Granville's Conjecture 10, a fixed power beyond its critical exponent $1/2$.
> (ii) Hence any interval of length $\sqrt q$ contains $O(1)$ roots of $x^2\equiv a\pmod q$, for every $a$
> (their Conjecture 9). (iii) Hence Rudin (their Theorem 11 in arXiv v1). Full $(L_8)$ also implies that an
> invertible residue class modulo $q\ge n^{1/4-\kappa}$ contains $O(1)$ divisors of $n$, for one absolute
> $\kappa>0$; this is below the exponent $1/4$ of Lenstra and Coppersmith–Howgrave-Graham–Nagaraj. Under $V_8$ one
> may take $\epsilon=10^{-16}$ and $\kappa=10^{-17}$, so Vojta's conjecture for the single variety
> $\overline M_{0,8}$ implies Rudin's conjecture. We found no earlier derivation: Bombieri–Lang gives
> $\sigma(k)\ll k^{5/7}$ (Cilleruelo–Granville). None of these consequences is known to fail, so they measure the
> strength of $(L_8)$, not its plausibility. They use the fact that $(L_8)$ does not see the norm form. The
> uniform bound itself is not known to imply Rudin's conjecture.

Also amend the Question A box by adding at the end: "Any proof of $(L_8)$, even restricted to pairwise coprime
square blocks, proves Rudin's conjecture."

**§5, item 12 (round 31), new bullet (proposed).**

> * *Hardness of Question A* (Part VIII; refereed). $(L_8)$ restricted to pairwise coprime square blocks implies
>   Cilleruelo–Granville's two-residue statement at length $q^{1/2+\epsilon}$, their Conjecture 9, and Rudin's
>   conjecture. Full $(L_8)$ implies $O(1)$ divisors of $n$ in an invertible class modulo $q\ge n^{1/4-\kappa}$.
>   This is new to the record, and we found it nowhere in the literature. A polynomial-in-$|S|$ count of
>   exceptional subspaces above an **absolute** height cutoff is false in every fixed dimension, and in dimension
>   five even with equal norms. Subquadratic counts above a cutoff **relative** to $\prod_{p\in P}p$ would give
>   $M=o(\log R/\log\log R)$; these are not standard conjectures, and none is derived from abc or Bombieri–Lang.

**Correction elsewhere.** In `overview/overview.md` §1 item 5 (line 180), `lit18/direct.md` §2.10 and
`local_smoothing/chains.md` row F, replace "arcs ⇒ Rudin-type statements" with "two congruence conjectures (their
Conjectures 9 and 12) ⇒ Rudin; the arc conjectures 13–15 are linked only to an $L^4$ conjecture for squares
(Theorem 16)".

**Sentences to avoid.**
* "The uniform bound (or the lattice-point conjecture) implies Rudin's conjecture." This is not known.
* "Cilleruelo–Granville show that arcs imply Rudin." This misdescribes Theorem 12.
* "Vojta's conjecture is known to imply Rudin." Say "we found no earlier derivation".
* "Question A is equivalent to Rudin" or "to a two-residue problem." Only one direction is proved.
* "Part VIII makes (L_8) implausible" or "is evidence against Vojta."
* "(L_8) crosses Lenstra's threshold" without "conditionally" and "a cardinality bound, no algorithm".
* "The lattice form (L_8^pc) implies Rudin." This is not shown.
* "Theorem 12" as the citation for Conjecture 9 ⇒ Rudin, without naming the edition.
* "Part VIII gives growth from a weaker hypothesis."

## 8. The standing sentence (R11)

The note says: "no intermediate growth bound is known from a weaker standard conjecture." Part VIII does not
change this.
* VIII.4 is a negative result. The polynomial-|S| exceptional-subspace hypothesis with an absolute cutoff, one of
  round 8's candidate weaker hypotheses, is false. (Round 8 had shown this only for m ≤ 3.)
* VIII.5–6 are conditional on relative-height counting hypotheses: f(s) = o(s²) subspaces, or F_δ(P) polynomial
  in |P| with an exponent independent of δ. These are not standard conjectures. Astra states that she cannot
  derive them from abc or Bombieri–Lang, and she does not compare them with Vojta.
* So the sentence stands as written.

I recommend adding a second sentence on the hardness side: "In the other direction, the sharpest Vojta-type form
of the problem, $(L_8)$, implies Rudin's conjecture on squares in arithmetic progressions (§2); the uniform bound
itself is not known to imply any recognised open problem." The second clause is the corrected round-4 verdict and
still holds.

## 9. Method and computations

* `m1.py`, check (1), Computed. Builds random q with 1–5 prime-power factors, two residues per prime power
  (coalescing at a random level, or coinciding), and eight distinct CRT lifts. It then forms Astra's reduced grid,
  columns P_i and square blocks independently of her script, and verifies the following on 79,632 pairs with no
  failure:
  * primitivity and P_0 = (1, 0);
  * the two-class unit-difference property after reduction;
  * det(P_i, P_j) = t_ij ∏_{T∋i,j} n_T;
  * the exact identity t_ij² = e^{D_ij}(y_j − y_i)²/Q².
* `m1.py`, check (2), Computed. Exact equality of the Cilleruelo–Córdoba and CG Theorem 10 exponents for
  k = 1, …, 59 (R6).
* `m2.py 3000`, Computed. The maximum over 2 ≤ q ≤ 3000 and all a of the number of roots of x² ≡ a (mod q) in a
  half-open real interval of length √q is 3, first attained at (q, a) = (105, 46). It is the same for unit a.
* By hand I re-derived VIII.1:
  * the Gram bound;
  * E m_A² ≤ E G² for even |A|;
  * the 127-term Markov step;
  * D_ij ≤ (1+a)W/2;
  * the 28-pair sum against cW/128 − O(1);
  * the Schwartz–Zippel escape, with a polynomial depending only on Z.

  All choices precede the modulus, so no exceptional set depends on q or S: the standing caution is respected. I
  also re-derived the diameter bound and the cofactor inequality of VIII.2.
