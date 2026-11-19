# Referee report, round 31, Part VII, block "meaning" (analyst): what Part VII changes, and what to do next

Sources read: Part VII of `astra_free_VII/answer.md` (lines 1461-1665), `reply7.txt` and `final7.txt`; the record
(`expert_note/expert_note.md` §5 item 12, `round31/README_VI.md` §§1.5, 3(f), 3(g), 5, `README_IV.md` §0,
`round10/README.md` and `round10/determinant.md` §3, `round8/weakest.md` §4, `growth/auxiliary.md` Lemma B,
`two_thirds/proof.md` Lemma 4.2 item 4); and this round's other Part VII referee reports
(`ref7_cap2.md`, `ref7_ramified.md`, `ref7_obstruct.md`, `ref7_census.md`, `ref7_hadamard.md`), which I use for
correctness. My own scripts and logs are in
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref7_meaning/`
(`sep_check.py`, `run_m3_L469.log`, `run_m3_L469_pad.log`). All the arithmetic is mine: F_49 = F_7(i) polynomials
are pairs of numpy integer arrays, with my own Euclid. I did not import Astra's scripts.

## 0. Verdicts in brief

| item | verdict |
|---|---|
| (a) Updated filter (F3+++ below) | CORRECT WITH CORRECTIONS. Cap two, finite ramification and "conjugation is not a coefficient Frobenius" are no longer escapes. This rests on VII.1 being correct, as cap2 and ramified find. The quantitative threshold stays at the ternary scale: an argument valid there cannot prove f(W) = o(log W/loglog W). It can still prove anything above that, including far more than the integer growth target. One escape remains provisional: ramified infinity is computed (m = 3) but not proved in Part VII. |
| (b) Integer meaning of "squarefree in odd characteristic is open" | No argument in the record uses squarefreeness and odd characteristic together in an essential way. Every squarefree-specific step in the record is combinatorial, and each holds verbatim over F_2[t] with squarefree N, where clusters are unbounded. The one route that turns squarefreeness into a strong conclusion (route 5, Thm 5.7) gets its power from ZE, which is archimedean and transcendence-type. Its target O(log W) is not even contradicted by the odd-characteristic families. **The last algebraic escape has no client among the known strategies.** |
| (c) Is the VII.1 device portable? | No, and it does not need to be. Separation (squarefree, disjoint T - v) is generic in characteristic 0. The device only repairs an inseparability defect that exists in characteristic p. The mechanism is still the VI.4 identity P_V = A_j^q - c_j^{q-1}A_j for every kernel, i.e. one parent with many same-degree right composition factors (Lemma 1 of README_VI §1.5). Confirmed by hand and by my own computation. |
| (d) Next question | **Recommend closing the function-field thread with a write-up.** Estimates are in §4. If one more exchange is spent on this thread, spend it on the write-up (first step in §4), not on squarefree-ternary or the "fixed exponent-one fraction" variant. |
| (e) Literature | Nothing found on worst-case short-arc clusters, or on unbounded same-norm clusters of bounded multiplicity, over F_q[t]. The F_q[t] sums-of-two-squares literature is statistical. Its standard model A^2 + TB^2 has **ramified infinity**, which is exactly the case Part VII covers only provisionally. |
| (f) Expert-note wording | §5 below. |

## 1. (a) The filter after Part VII

### 1.1 What VII.1 removes, and checks against the standing cautions

VII.1 is CORRECT according to `ref7_cap2.md` (constant model, every p ≡ 3 mod 4, certificate expanded) and
`ref7_ramified.md` (ramified model, census reproduced). I read the proof myself (answer.md 1476-1516) against the
project's standing cautions:

* **gcd costs.** The marks are Z_j = F_j/G(T), the degrees used are those of the reduced marks
  (D_0 = (n-g)L), and the width estimate is made "after precisely the formal division by G". Padding multiplies by
  a real e, which raises D and the parent exponent to exactly 2 on new primes. OK.
* **Distinct, equal degree, one norm in the base ring.** Distinctness comes from the disjoint zero sets of the T - v
  together with two nonproportional forms. Equal degree is (n-g)L. The norm is in F_p[t] by descent (conjugation
  negates labels; the constant automorphisms permute coefficient vectors). OK.
* **Conjugates and units.** The marks are monic and distinct. A pair z_i, unit·z̄_j would still count as two
  distinct lattice points, as in VI.4, which the judge accepted.
* **Infinity inert.** Yes in both models of the statement. Ramified infinity is an extension; see 1.3.
* **Thresholds depending on m.** D ≥ C_{p,Δ} 2^m n^2 depends on m. The statement is a liminf in D for each m, not a
  uniform claim. OK.

My own small instance (`sep_check.py`; p = q = 7, k = 2, n = 49, h = 7, C = t^7 - t, T = 1 + C^7 t^420,
L = 469, forms x, y, x+y, all signs +):
* every T - v is squarefree, and the 49 factors are pairwise coprime;
* R = ∏(T - v) has degree 22981 and lies in F_7[t];
* the three marks have degree D = 16884, one norm in F_7[t], and max deg(z_i - z_j) = 8211 ≤ D/2 (slack 462);
* the parent exponents on the T - v are 0 (7 labels), 1 (12 labels) and 2 (30 labels), so the cap is 2.

At L = 109 the same script gives max deg = 2811, exactly Astra's bound ((n+h-2g)L + b(n-h))/2 with g = 13, b = 56.
This bound is not critical there, so the stated width bound is sharp at this size too.

### 1.2 The new filter (proposed wording, "F3+++")

> **F3+++.** Call an argument *function-field valid* if every step holds verbatim, with deg in place of log and the
> place at infinity in place of the archimedean place, in **either** of the following settings:
> * (i) F_2[t] with the norm of F_4(t)/F_2(t), for all N or for squarefree N only (Parts V, VI.5);
> * (iii) F_p[t] for any fixed p ≡ 3 mod 4, with **either** the constant extension (A^2 + B^2) **or** a ramified
>   order F_p[t, η], η^2 = Δ, with deg Δ even and positive, lc(Δ) a nonsquare, at least one finite ramified prime,
>   conjugation η ↦ -η fixing the constants, and infinity inert. Here N ranges over parents **all of whose prime
>   exponents are at most 2** (and a fortiori over all N). The unpadded parents are prime to the discriminant.
>
> Then no function-field-valid argument proves the uniform bound M ≤ M(C). Quantitatively:
> * an argument valid in (i) cannot prove M ≤ f(W) with f(W) = o(sqrt(W/log W)), even for squarefree N;
> * an argument valid in (iii) at a given p cannot prove M ≤ a·log W/loglog W for any a < (1/2)log_2 p, at
>   **every** sufficiently large W (all large D are covered, not a sparse sequence), even with exponents capped at 2,
>   finite ramification, and a conjugation that is not a coefficient Frobenius;
> * an argument valid in (iii) for all p ≡ 3 mod 4 with p-independent constants therefore cannot prove
>   M = O(log W/loglog W). *(Heuristic sharpening, mine, unrefereed: letting p ≈ 2^m with k = 2 in VII.1's
>   threshold suggests even M = o(log W) is out of reach for such p-uniform arguments, with roughly M ≥ (1/6)log_2 W.
>   This is from Astra's stated threshold and not checked in detail. Do not record it as proved.)*
>
> An argument valid in (iii) can still prove, for example, M ≪ W^ε or M ≪ log W. **Growth-rate arguments at the
> integer scale (W/log W) are not filtered at all.**

Setting (ii) of F3++ (constant extension, all N) is now contained in (iii). An argument that needs "exponents at
most 2", "finitely ramified", or "conjugation not a Frobenius on coefficients" is filtered along with all the earlier
items.

**Admissible, hence insufficient (additions to F3++'s list):**
* any fixed multiplicity cap e_0 ≥ 2;
* parents prime to the discriminant;
* ramification at a finite prime, with a ramified prime that need not be principal (`ref7_ramified.md`: the class
  number is even and the order is never a PID with inert infinity);
* a conjugation fixing constants.

### 1.3 What remains available to a proof over Z

None of these is known to suffice, and the list is in my order of relevance:

1. **Characteristic 0 used essentially.** This covers Newton linearisation, BMY (Part IV: ≤ 24 over C(t)), and
   Lemma 1 (uniqueness of same-degree right composition factors). Since VI.4 every function-field counterexample
   rests on q-linearised polynomials. This is the only item with a working proof behind it, and only in the
   function-field model.
2. **Archimedean order, curvature, real-place balance; Baker/transcendence input** (route 5's ZE is of this type).
3. **Exactly squarefree N in odd characteristic.** Open. Data:
   * M_3^sf(D,0) = 5 at D = 12, 14, 16, 18 and 4 at D = 13, 15, 17, 19 (`ref7_census.md` extended the table to
     D = 19);
   * it recurs at every D = 12L;
   * every specialisation of the original VI.4 marks fails (VII.2, CORRECT per `ref7_obstruct.md`; numerically no
     squarefree case for any m ≥ 2 with T nonconstant);
   * the unrestricted census showed M_3(D,0) ≤ 4 for D ≤ 11 while sup M_3 = ∞, so **small-D census data are no
     guide**.
4. **A fixed positive fraction of deg N on exponent-one primes** (Astra's proposed question). This is open. In VII.1
   the exponent-two primes are exactly the mixed labels, and the exponent-one share is O(m/q) = O(m·2^{-m}) → 0.
   My q = 7, m = 3 instance has share 1/6, but such a fixed small instance says nothing about m → ∞.
5. **Ramified infinity (deg Δ odd, e.g. A^2 - tB^2, O = F_3[u] with u ↦ -u).** This is *provisionally* no longer an
   escape. `ref7_ramified.md` computed q = 27, m = 3 (critical at D = 210210) and states that the proof transfers
   line by line: anisotropy holds because the two leading orders differ by a half-integer. It is not proved in
   Part VII. Record it as "computed (m = 3); proof expected to transfer; not refereed as a theorem".
6. **p ≡ 1 mod 4** (χ(-1) = +1 breaks F_j σF_j = P_V^2; the half-system obstruction). This is a calibration only:
   -1 is a square in the residue field there, which has no clean integer meaning.
7. **Characteristic-2 ramified (Artin–Schreier) models.** Untouched.

Item 3 is the only genuinely Z[i]-shaped algebraic escape left. Item 1 is the only one with a theorem behind it.

## 2. (b) What "squarefree in odd characteristic remains open" means for Z

**How much of the integer problem is squarefree?** For the growth rate, almost all of it:
* every recorded fake is squarefree (`round8/README.md` §2.4, `weakest.md` §4);
* the hard growth-regime profiles are squarefree sign profiles;
* the 2/3 proof's chain lemma loses only the factor c/λ = 1 - O(λ) against **the squarefree model**
  (`two_thirds/proof.md` Lemma 4.2 item 4), so squarefree is the *extremal* case of the proved bound, not an easy
  case;
* abc over Q(i) is **vacuous on squarefree profiles** (`weakest.md` Prop 4.2).

For uniform boundedness, a counterexample over Z need not be squarefree. The record's natural low-rank adversaries
(`determinant.md` §3) need e_p ≥ 2. So "n squarefree" is a real sub-case, and for growth it is the core of the
problem.

**Does any integer argument in the record use squarefreeness essentially, so that the open F_3 case could support
it?** I grepped every occurrence of "squarefree" outside round 31. The squarefree-specific steps are:

| step | what it uses | valid verbatim over F_2[t], squarefree N? |
|---|---|---|
| Lemma 3.1, `round10/determinant.md` (rank ≥ n/4 - 2) | Gram matrix of ±1 sign columns plus the pair gap (I4) | yes: real linear algebra and degree bookkeeping |
| Prop 3.4, `round8/stability.md` (ω(N') > M/4) | counting | yes |
| Lemma B, `growth/auxiliary.md` (height isolation) | triangle inequality on exponent distances | yes |
| Lemma A5, `growth/auxiliary.md` (built-in divisibility is zero) | a *loss*, not a gain | n/a |
| `round7` dichotomy (squarefree members) | residues and pinning (P3-G) | local and characteristic-free |
| abc / n-term abc, `weakest.md` Prop 4.2 | shows abc is *vacuous* there | a negative statement |
| Plotkin M ≤ 2D (`round29`, over R(t)) | characteristic 0 | not relevant to odd p |
| route 5, Thm 5.7 (M = O(loglog R) on squarefree norms under ZE) | Lemma 3.1 plus Thm 4.1 plus **ZE** on the angle torus | ZE is archimedean and transcendence-type |

Consequences, stated bluntly:
* Every squarefree-specific step that is proved is valid in setting (i). There, squarefree clusters are unbounded
  with M ≍ sqrt(W/log W). So none of these steps can bound M even in combination. None of them uses "2 ≠ 0"
  together with squarefreeness in a way that would separate F_3 from F_2.
* Route 5 is the only strategy whose conclusion depends on squarefreeness. Its non-filtered input is ZE, which is
  not an algebraic statement about F_p[t]-type arithmetic. Its target M = O(log W) is **not contradicted by any
  odd-characteristic family**: VI.4 and VII.1 have M ≍ log W/loglog W = o(log W). So even an unbounded squarefree
  ternary family at the VI.4 rate would not test route 5. Only the binary families (M ≍ sqrt(W/log W) > log W) do,
  and route 5 must already use something invalid in characteristic 2.
* Hence **the last algebraic escape is of no use to any known strategy.** If sup M_3^sf were proved finite, the
  proof would have to use odd characteristic or reflection conjugation essentially (since F_2 squarefree is
  unbounded), together with squarefreeness. It would be the first uniform theorem in a positive-characteristic
  model and would be interesting in itself. But nothing in the record is waiting for it. If the supremum is
  infinite, the filter extends to squarefree N, and the integer problem is again unaffected, because no integer
  argument was relying on that escape.

## 3. (c) Portability of the VII.1 device; the mechanism is still Lemma 1

**The device.** The labels v ∈ η^p C^p F_q[t^p] are p-th powers, so dv = 0. For each v, d(T - v) = T' dt, and T' ≠ 0
off the collision set {C = 0}, where every v vanishes and T = 1. Hence the T - v have simple, pairwise disjoint
zeros, and R = ∏(T - v) is squarefree. My check confirms this for p = q = 7 (§1.1).

**Portability.**
* In characteristic 0 nothing needs separating. For a polynomial family with simple-rooted P_V, squarefree and
  pairwise coprime specialisations T - v are the *generic* behaviour (Hilbert irreducibility, squarefree values).
  Over Z they are the expected behaviour by sieve/CRT, or under abc for high degree.
* In characteristic p, by contrast, inseparability creates forced common zeros. Examples are the collisions of
  v - v' at the zeros of t^{q^j} - t, and the multiple zeros where d(T - v) = T' - v' vanishes. These are the
  multiplicities VII.2 measures.
* So the device imports a characteristic-0 feature (generic separability) into characteristic p. There is nothing
  to export: the integer side already has separation and lacks the formal family.

**Mechanism.** VII.1 reuses VI.4's formal identities unchanged:
* P_V = A_j^q - c_j^{q-1}A_j for **every** kernel polynomial A_j = P_{ker ℓ_j} (all of the same X-degree h), and
* F_j = A_j(A_j^s - c_j^s)^2 = P_V + 2c_j^{q-1}A_j - 2c_j^s A_j^{(q+1)/2}, which gives
  deg_X(F_i - F_j) ≤ (n+h)/2.

I verified P_V(T) = A_j(T)^7 - c_j^6 A_j(T) exactly for all three forms in my instance (deg A_j(T) = 3283,
deg P_V(T) = 22981). So the parent P_V has m different right composition factors A_j of the same degree, with the
outer polynomial Y^q - c^{q-1}Y. By Lemma 1 (README_VI §1.5; Astra repeats the argument at answer.md 1520) this is
impossible when char ∤ deg f, i.e. always in characteristic 0. Confirmed: **Part VII changes the label space and the
specialisation, not the mechanism.** Every unbounded family on record (V.1, VI.4, VI.5, VII.1) rests on additive
(q-linearised) polynomials.

## 4. (d) Next question: estimates and recommendation

The probabilities below are my own guesses. They are calibrated against this thread's history: seven exchanges, each
settling the question posed, with constructions arriving much faster than obstructions.

| candidate for the next exchange | P(decided in one exchange) | P(the answer affects the integer problem) | comment |
|---|---|---|---|
| Astra's: cap two with exponent-one primes carrying a fixed fraction δ > 0 of deg N (normalised clusters) | ~15% (almost all of it "unbounded by a new device") | < 0.1% | A pure model refinement. It **must** be posed for normalised clusters (gcd 1). Otherwise a common nonreal factor of degree up to the slack adds exponent-one primes for free (`sep_check.py` padding run: still critical after a common factor of degree 120). In VII.1 the slack is only ≈ 2^{1-m}·D, so this alone cannot give a fixed δ. |
| squarefree ternary: is sup_D M_3^sf(D,0) finite? | ~10% (≈ 8% construction, ≈ 2% proof) | ≈ 0.3% (only through a bounded answer whose proof then ports; §2 shows no known strategy is waiting) | VII.2 rules out the natural route; the census is no guide. |
| binary growth scale (W(M) between M log M and M^2 log M over F_2) | ~5% | ≈ 0.3% | The only function-field question at the integer *growth* scale, but it is a filter only. VII.4 and the two-core computations close classes without trending. |
| direct integer question using the lesson (e.g. an arithmetic substitute for "uniqueness of same-degree right composition factors") | < 5% that it is even well posed in one exchange | unclear | Over Z there is no composition structure to be unique. Part IV already says the integer counterpart class (polynomial templates) is ≤ 24, and the record identifies the integer substitute for BMY as Szpiro/abc strength, which is vacuous on squarefree profiles. This is the barrier the record has documented since round 8; another exchange will not move it. |
| abc/Vojta reformulations | low | low | Already in the record (barrier 6; LW support test equivalence in Part VI). |

**Recommendation: close the function-field thread.** The thread has done its job. The filter F3+++ now states
precisely that a uniform proof over Z must use characteristic 0 (as Part IV does), archimedean or transcendence
input, or exact squarefreeness in odd characteristic. No known strategy uses the last of these (§2). Every further
model refinement (cap with δ, squarefree, p ≡ 1 mod 4, Artin–Schreier) has integer relevance near zero.

**If the coordinator wants one more exchange with Astra, make it the write-up, not a new question.** Ask for a
self-contained note, "Unbounded short-arc clusters for quadratic norm forms over F_q[t]", with these first steps:
1. one theorem statement covering V.1 (F_2), VI.4/VI.5 and VII.1 (constant, finitely ramified, and, with a
   complete proof, **ramified infinity, deg Δ odd**, which is the A^2 + TB^2 model of the statistical
   literature), with exact hypotheses: p ≡ 3 mod 4, nonsplit infinity, cap 2 or squarefree (F_2);
2. the mechanism isolated as a lemma: one parent with m same-degree right composition factors gives an m-point
   cluster after specialisation, plus Lemma 1 for its characteristic-0 impossibility;
3. the explicit open list: exactly squarefree odd characteristic (with the census to D = 19), fixed exponent-one
   fraction, p ≡ 1 mod 4, Artin–Schreier.

Writing the ramified-infinity proof out is the one item that would turn a provisional filter entry into a theorem.

If the coordinator insists on a research question instead, squarefree ternary is the only one with any integer
relevance (≈ 0.3%). Its first step is a census of squarefree parents carrying **5** marks at D = 12-19 (now known),
inspecting whether they share a composition structure (P = f(A) = g(B) with deg A = deg B) or are sporadic. Do not
pursue the δ-fraction variant.

## 5. (e) Literature

Searches (no personal data used): "lattice points on short arcs function field F_q[t] Cilleruelo Granville";
"sums of two squares in short intervals polynomial ring finite field Bank Bary-Soroker Fehm"; "polynomials with the
same norm close together, linearized subspace polynomial, quadratic character, norm equation".

* No function-field analogue of the Cilleruelo–Granville or Cilleruelo–Córdoba worst-case short-arc problem was
  found. No construction of many same-norm elements of F_q[t, η] at small mutual distance was found either,
  bounded multiplicity or not.
* The F_q[t] literature on sums of two squares is **statistical**:
  * [Bank–Bary-Soroker–Fehm](https://arxiv.org/pdf/1509.02013): density in short intervals, large q;
  * [Bary-Soroker–Fehm](https://pdf.arxiv.org/pdf/1701.04092): autocorrelations;
  * Gorodetsky–Rodgers: variance;
  * [Yiasemides](https://arxiv.org/pdf/2309.01290): UE^2 + VF^2, thin annuli.

  These count representations on average. They do not bound the worst case.
* Their model is A^2 + TB^2 = norm from F_q[sqrt(-T)], which has **ramified infinity**. That is exactly the case
  VII.1 does not state and `ref7_ramified.md` covers only by computation at m = 3.
* Background on linearised and subspace polynomials (Ore; Ben-Sasson–Kopparty–Radhakrishnan; Gow–McGuire on Galois
  groups of linearised polynomials over F_q(t)) is classical. The combination is new as far as searches show (search
  snippets only; no full texts read).
* Integer side: [Cilleruelo–Granville, Close lattice points on circles](https://www.cambridge.org/core/journals/canadian-journal-of-mathematics/article/close-lattice-points-on-circles/E14DBCBF07B4B98FECF3282E2BA672B5)
  and [Lattice points on small arcs](https://arxiv.org/pdf/2107.09991) concern Z only.

**Novelty against the record (grep):**
* cap two / bounded multiplicity was recorded only as **open** (expert note lines 913, 958; README_VI caveat 1);
* the ramified model was "untested" (README_VI §5 item 2), and A^2 - tB^2 was "outside the model"
  (`ref6_ternary.md` 220).

So VII.1, the VII.3 tables and VII.4 are new. The separation device is the new idea; the rest is VI.4.

**Significance for the integer goal: nil.**
* The growth rate is unchanged at (2/3+ε) log R/loglog R, and the uniform bound M(C) is unchanged.
* Questions A and B and route 5 are unchanged.
* Part VII removes three escape clauses (cap two, finite ramification, non-Frobenius conjugation) that no integer
  argument was using.

## 6. (f) Wording for the expert note

**Part VII bullet (proposed, for §5 item 12):**

> ***Part VII: cap two, ramified models*** (`round31/README_VII.md`). (1) **Proved** (VII.1; referees cap2,
> ramified): for every p ≡ 3 mod 4, in the constant extension (A^2+B^2) and in every ramified order F_p[t,η],
> η^2 = Δ (deg Δ even > 0, lc Δ nonsquare, a finite ramified prime, conjugation fixing constants, infinity inert),
> critical clusters at c = 0 on parents with **every prime exponent ≤ 2** are unbounded:
> M^{(≤2)}(D) ≥ ((1/2)log_2 p - o(1)) log D/loglog D for every large D. The new device is to take labels in
> η^pC^pF_q[t^p] (p-th powers, so dv = 0), which makes the T - v squarefree and pairwise coprime. The mechanism is
> unchanged from VI.4 (one parent P_V with a same-degree right composition factor for every kernel; impossible in
> characteristic 0 by Lemma 1). The exponent-two primes carry a share 1 - O(m·2^{-m}) of deg N. Ramified infinity
> (deg Δ odd, e.g. A^2 - tB^2) is computed at m = 3 and the proof is expected to transfer (not refereed as a
> theorem). (2) **Proved** (VII.2): no specialisation of the VI.4 marks (m ≥ 16, q ≥ 27) has squarefree norm after
> any common division; sign editing and kernel or translate replacements fail for q ≥ 11. (3) **Computed**
> (exhaustive, reproduced): M_3^sf(D,0) = 2,4,3,4,3,4,4,4,4,4,4,5,4,5,4,5,4,5,4 for D = 1..19 (five first at D = 12,
> recurring at D = 12L); in the ramified order η^2 = 2t^2+1, all parents 2,3,2,4,3,4,3,4,4,4 and squarefree
> 1,2,2,3,3,4,3,4,3,4 for D = 1..10. (4) **Proved** (VII.4): over F_2[t], a Hadamard matrix of order m ≥ 8 with all
> rows and a missing nonconstant column cannot give a critical cluster (block form, disjoint supports; false at
> m = 4); the eight-row x, x^3 two-core profile is impossible at e = 10..18. **Open:** squarefree N in odd
> characteristic; the binary growth scale. **No change for the integer problem.**

**Filter:** replace the F3++ paragraph with F3+++ (§1.2), including both quantitative thresholds and the sentence
"growth-rate arguments at the integer scale are not filtered". Replace the "still available" list with §1.3, items
1-3, plus "a fixed positive share of exponent-one primes (open)", and mark ramified infinity as provisional.

Add one sentence to the record's squarefree discussion: "No squarefree-specific step in the record (Lemma 3.1 rank,
ω(N') > M/4, height isolation, the 2/3 chain lemma, abc vacuity) uses odd characteristic; all hold over F_2[t] with
squarefree N, where clusters are unbounded. Route 5's squarefree conclusion rests on ZE, not on squarefreeness."

**Sentences to avoid:**
* "bounded multiplicity does not help over function fields" without "p ≡ 3 mod 4, cap ≥ 2, nonsplit infinity";
* "squarefree N is unbounded in odd characteristic" (open);
* "false for sums of two squares A^2 + TB^2 over F_q[t]" (ramified infinity: computed at m = 3 only);
* "the census suggests squarefree boundedness" (the unrestricted census said ≤ 4 to D = 11 and was wrong about
  sup);
* "squarefreeness is the last hope of a proof over Z" (no integer strategy uses it; the remaining real inputs are
  characteristic 0, archimedean, transcendence);
* "evidence against the conjecture over Z";
* "the exponent-two primes are an artefact" (they carry a 1 - O(m·2^{-m}) share of deg N, and VII.2 shows they are
  forced for this family);
* "Part VII is a new mechanism" (same mechanism as VI.4; a new separation device).

## 7. Reproduction

```
cd /tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref7_meaning
prlimit --as=3500000000 python3 sep_check.py 420     # L = 469: critical, cap 2, identity checks (about 1 min)
prlimit --as=3500000000 python3 sep_check.py 60      # L = 109: identity checks, width bound attained, not critical
```
