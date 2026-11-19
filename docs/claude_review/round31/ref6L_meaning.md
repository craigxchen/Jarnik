# Round 31, Part VI (Question B thread): referee block "meaning" (analyst)

Source: Part VI of `astra_lemma_VI/answer.md` (l.1319-1563), `final6.txt`, `reply6.txt`. Record:
`expert_note/expert_note.md` §§1, 3, 5; `round31/README_VL.md` §§3(b), 5; `round20/gauss.md` §3;
`growth/families.md` §§1, 2, 5; `round18/theory.md` Thm 2.5; `lit18/tools.md`.
Scripts: `round31/ref6L_meaning/checks.py`. It checks the cut algebra and LP optimum behind VI.4/VI.6, Astra's
v_n relation (exact, no proper vanishing subsum at n = 10, 100, 1000), the two-term unit form of the golden block,
and Fibonacci gcd overlap. Downloaded text: `ref6L_meaning/dl/cz.pdf` (Corvaja-Zannier, arXiv math/0311030).

This block judges meaning, scope and novelty. I did check the VI.1 proof structure at the level needed for the
scope questions. That covers the Subspace inequality in the zero-sum step, the size reduction, the degree count
in (VI.5) and the integrality/product-formula contradiction, and I found no error. The line-by-line refereeing of
VI.1, VI.3 and VI.5 belongs to the other blocks.

## Verdicts

| claim | verdict |
|---|---|
| VI.6 forward: existence in (*) refutes A_5(delta) for every delta<1/10, not (L_5°), not M(C) | **CORRECT** |
| VI.6 converse: no Zariski-dense (*) sequence <=> A_5(delta) for SOME delta in (0,1/10) | **CORRECT** (with the precise quantifier wording below) |
| Reconciliation with README_VL 3(b) | consistent; VI.6 upgrades a one-way implication to an equivalence; "strictly weaker" must read "logically weaker (no separation known)" |
| VI.1 legitimacy w.r.t. the standing cautions | **CORRECT**, honest: K, S, b are fixed independently of the parent; Z, rho, H_0 ineffective |
| VI.1 novelty | **new to the record** as a theorem on genuine tuples. It is not a literal instance of a published theorem I could find. Its engine (VI.3) is Evertse 1984 / Corvaja-Zannier 2005 Prop. 1 with almost-S-unit coefficients |
| "VI.1 makes fixed-order unit templates unable to produce unbounded clusters" | **CORRECT WITH CORRECTIONS**: true only in a weak/trivial form; the correct unconditional content is a five-point no-go for constructions (Cor. U below) |
| (c) degree-uniform Subspace / field-uniform VI.3 | **FALSE** as a statement (v_n relation, checked exactly). The version with a discriminant term (Vojta, algebraic points) is conjectural and costs as much as it gains. Growing S gives nothing |
| VI.5 scope | correct as a ten-pair-profile statement only; it does not cover balanced profiles (unlike VI.1) |

## (a) VI.6 and README_VL §3(b)

Definitions (expert note §3): A_5(delta) says there are a proper closed Z and an R_0 such that every primitive
five-tuple on a circle of radius R''>=R_0 with all distances <= R''^{1/2-delta} has its class in Z. Any unit
classes are allowed. With H = R''^2 and Delta = max|q_i - q_j|, the condition is Delta <= H^{-1/4-delta/2}, which
at delta = 1/10 is H^{-3/10}.

**Forward.** Take a (*) sequence that escapes every proper closed set: eventually outside any fixed Z, with
Delta = H^{-3/10+o(1)}. For delta < 1/10, -3/10+o(1) < -1/4-delta/2 eventually, so the distances are at most
R''^{2/5+o(1)} <= R''^{1/2-delta}. That contradicts A_5(delta). (L_5°) needs |log n_A - w| <= eta_0 w on all 15
cuts, and the five singleton cuts of (*) are o(W), so no balanced tuple arises. Five points say nothing about M(C).
**CORRECT.**

**Converse.** Suppose A_5(delta) fails for every delta in (0,1/10). For each delta and Z there are then infinitely
many primitive tuples outside Z: finitely many could be added to Z, and bounded R'' gives finitely many tuples.
Enumerate the countably many hypersurfaces over Q, take delta_n increasing to 1/10, and pick the n-th tuple
outside the first n hypersurfaces with R''_n increasing to infinity. Then every proper Q-closed set contains
finitely many terms. Q-closed suffices: a geometric closed set containing infinitely many Q-points (Q(i)-points of
the twisted torus) contains them in a Q-closed proper set, by linear algebra on a fixed degree-d space of
equations. Two-sidedness and the profile are forced (my checks.py reproduces both):

* |q_i - q_j| >= e^{-d_ij/2}, by xi C - C-bar being a nonzero Gaussian integer. This holds across unit
  classes.
* sum d_ij = 6W - 2W_sing (symbolic check).
* max_{w>=0, W=1} min_ij d_ij = 3/5 (LP). At the optimum the solution is unique: singletons 0, pairs W/10.

So Delta = H^{-3/10+o(1)} and the ten-pair profile follow, and the sequence is a (*) sequence. Primitivity
excludes the primes 2 and q ≡ 3 (4), so the cut dictionary applies.

Exact statement to record: *"There is no sequence of primitive five-tuples in the ten-pair-cut profile with
Delta = H^{-3/10+o(1)} that eventually leaves every proper closed subset of T_5" is equivalent to "A_5(delta)
holds for some delta in (0,1/10)", hence for all delta in (delta_0,1/10).* (A_5 is monotone in delta.)

**Reconciliation.** README_VL §3(b) said "non-existence off Z is implied by A_5(delta) for any delta<1/10, much
weaker than Question B". VI.6 adds the converse, so non-existence in (*) *is* A_5 just below the
Cilleruelo-Córdoba exponent. It is logically implied by Question B's A_5(delta<1/30), because A_5(delta) implies
A_5(delta') for delta'>delta. "Strictly weaker" cannot be asserted: no separation is known. It is independent of
(L_5°) as far as is known, since the profiles are disjoint. Calibration that matters: A_5(delta) for *some*
delta<1/10 is an unconditional improvement of Cilleruelo-Córdoba at five points. gauss.md calls that "the first
improvement ... since 1992", and the record found no unconditional theorem beating Liouville (tau_e = 10/3) on a
Zariski-dense set at a rational point of any variety (expert note §3, round 22 J3). **So (*) is not an easier
sub-question of Question B. It is the "beat Liouville at all" problem, with the worst (unbalanced) profile.**

## (b) What VI.1 and VI.5 really give for genuine tuples

**First genuine-tuple statements with exponent >= 4 off a closed set: yes.** The record's genuine-tuple results
are of two kinds:
* Liouville/Seshadri-type bounds (exponent <= 10/3);
* family exclusions for bounded-degree curves: gauss Thm 4.3(ii), census Prop 7.5, VI.3 here.

The record's Ridout theorem (round18 theory Thm 2.5, expert note §2 item 2) assumes **all blocks are S_0-units,
S_0 = the parent primes**, so its constant C(S_0) depends on the parent and it "is useless here". lit18/tools.md
classifies all fixed-S gcd/Vojta theorems (CZ, Levin, Xiao, ...) the same way, "S containing the primes of N".
VI.1 inverts the roles: S consists of the fixed places of the **additive ingredients**, and the parent primes
(divisors of G_e = sum a_j u_j) are completely free and outside S. Hence the exceptional set depends only on
(K,S,b,eps), not on the parent. Against the standing cautions:
* No uniform exceptional set across varying S is used: S is fixed by the class.
* No gcd is discarded: the hypothesis prod|G_e|^2 <= H is required. Note: in general the theorem reads
  Delta <= (prod_e|G_e|^2)^{-1/4-eps}, so block overlap costs exactly as much as it should.
* No drifting phase is set to zero: eta_i of height rho log H stay in the coefficient field k.
* The threshold H_0 depends on (K,S,b,eps) only.

All of Z, rho, H_0 are ineffective (compactness/diagonal argument plus Subspace), and Astra says so. This is
honest.

**What it covers.** VI.1 holds for **all fifteen cuts**, so on its class it gives A_5(delta) for every delta>0,
hence (L_5°) by gauss Prop 3.1, which is pointwise. That is exponent 4, the Vojta prediction, above the needed
15/4. It covers:
* golden blocks. J_r = f_{r+1}+i f_r = a phi^r + b phibar^r with a = (phi+i)/sqrt5, b = -(phibar+i)/sqrt5
  (checks.py, r<12 exact): b = 2, K = Q(i,sqrt5), S = archimedean places;
* silver/Pell blocks, all `units_search.py` recurrences and Pell orbits, in finitely many fixed orders;
* exponents chosen independently per block, and seeds/linear maps of subpower height.

**How large is the class compared with what the conjecture needs.** The number of Gaussian integers of norm <= X
that are sums of <= b S-units with coefficients of height <= rho log X is X^{O(b rho)} (log X)^{O(b)}. This is a
subpower-density set, and every Gaussian integer is a one-term sum whose coefficient has full height (1/20) log H.
In the reduction A_5 => uniform bound (gauss Thm 3.4/3.4′, expert note §3), the five-tuple is an arbitrary
gcd-rich subtuple of an arbitrary cluster on an arbitrary circle. Nothing gives it this structure. **For the
uniform conjecture VI.1 is conditional on structure. It is a no-go theorem for constructions, not a step towards a
proof.**

**The unconditional statement, formulated precisely (Cor. U).** Fix K ⊇ Q(i), finite S, b, and a finite
coefficient set F ⊂ K^×. Let a five-tuple be built from a "template": z_k = eps_k gamma prod_{s=1}^m
J_s^{w_{ks}} J̄_s^{1-w_{ks}}, with each J_s a sum of <= b elements of F·O_{K,S}^× (exponents free), and the
nonconstant J_s pairwise coprime up to H^{o(1)}. Assume moreover m <= m_0 fixed. Then for every eps>0 there is a
proper closed Z (ineffective) such that all such primitive five-tuples with distances <= R''^{1/2-eps} lie in Z
once R'' is large. Proof: group the J_s by induced cut, so each G_e is a product of <= m_0 bounded sums, which is a
bounded sum. The coprimality makes prod|G_e|^2 = H^{1+o(1)}. Apply VI.1.

Corrections to the "unbounded clusters" reading:
1. With m <= m_0, a template has at most 2^{m_0}·4 points, so "cannot produce unbounded clusters" is trivial.
   The content of Cor. U is that **no Zariski-dense five-point family at arc R^{1/2-eps}, and hence no
   counterexample to Question B, (L_5°) or A_5(delta) for any delta>0, comes from bounded fixed-order unit
   templates with honest blocks.** This is the five-point analogue of families.md Thm A (fixed-shift quadratic
   templates: at most 4 rows with C -> 0). Thm A is stronger in arc (o(sqrt R)) and weaker in class (one unit,
   fixed shifts, one-parameter). Cor. U is new to the record.
2. Two hypotheses are necessary.
   * m bounded: products of unboundedly many bounded sums are unbounded sums, and the dimension d in the proof
     is then unbounded.
   * Subpower overlap: with free exponents, gcd(f_{2r+1}, f_{2s+1}) = f_{gcd}, which can be all of the smaller
     block (checks.py: 45 | 135).
   The record's resonant cyclotomic layers (families.md §1, `cyclotomic_five_row_subendpoint_family.md`) exploit
   exactly the overlap. Their honest cut blocks are cyclotomic factors Phi_d(lambda,lambda') with unbounded term
   count. They lie **outside** VI.1, so there is no conflict. This is inference: that note is not in the
   checked-out record.
3. Cor. U is ineffective, and Z depends on (K,S,F,m_0,eps). It does not give M(C) bounds.

**VI.5.** It applies only to the **ten-pair profile** (its proof uses the ten pair cubics and the II.5 two-moment
computation), not to balanced profiles. It is the five-point counterpart of the record's A_4 counterexample
mechanism (gauss Thm 4.1: Gaussian linear blocks H_T(s) = q_T(s-r_T), i.e. affine in one parameter). Read
together they say: **affine one-parameter blocks refute A_4 but cannot reach exponent 4 at five points in the
endpoint profile.** That is useful calibration and new. It has no consequence for Question B, whose profile is
balanced.

## (c) Weakening the fixed-field hypothesis

* **Bounded degree.** a+ib = (v_a+v_a^{-1}) + i(v_b+v_b^{-1}) with v_n = (n+sqrt(n^2-4))/2. So with b = 4 and
  degree <= 8, *every* Gaussian integer is covered. A VI.1 uniform over fields of degree <= D would therefore
  prove A_5(delta) for all delta>0 for all tuples, and with it the uniform conjecture (gauss Thm 3.4). That would
  be far too strong to be cheap. It is false at the VI.3 level: v_{n+1}+v_{n+1}^{-1}-v_n-v_n^{-1}-1 = 0 has no
  proper vanishing subsum (exact, n = 10, 100, 1000) and pairwise ratio heights ≍ log n.
  The same points also show a "degree-uniform Subspace Theorem with fixed exceptional subspaces" is false:
  * The forms are fixed (coordinates, with the largest one deleted at each place). With S the archimedean
    places, the points (v_{n+1}, v_{n+1}^{-1}, -v_n, -v_n^{-1}) satisfy the Subspace inequality with margin about
    H^{-5}.
  * Yet no fixed linear relation among them holds for infinitely many n. The leading terms force all
    coefficients to vanish.
  * The absolute Subspace Theorem (Evertse-Schlickewei, Crelle 548 (2002)) does not apply. Its forms at places
    of Q̄ must be Galois-conjugates of forms attached to places of the fixed field of definition. The
    point-dependent choice "delete the largest coordinate at each embedding of K_n" is not allowed. This is my
    reading of its hypotheses, not re-derived from the paper.
* **Counting uniformity.** Evertse-Schlickewei-Schmidt (Ann. Math. 155 (2002)) bounds the number of nondegenerate
  unit-equation solutions in terms of n and the rank only. That is uniform in field and coefficients, but each
  v_n field contributes about one solution. A count per group says nothing about a sequence that changes group, so
  no structure follows.
* **Conjectural uniform version.** Vojta's general conjecture for algebraic points of bounded degree carries the
  discriminant term d(P) = log|disc|/[K:Q]. For the v_n representation, log disc Q(i, sqrt(a^2-4), sqrt(b^2-4)) is
  ≍ log|G_e| ≍ (1/20) log H per block, the same order as the height being saved. The Vojta-with-truncation /
  "uniform Schmidt" forms are conjectural and give nothing here. **Known: no. Conjectured: yes, with a
  discriminant term. Useful: no.**
* **S growing with H.** In the quantitative Subspace Theorem (Evertse-Schlickewei; Evertse-Ferretti, Ann. Math.
  177 (2013)), the number of exceptional subspaces is bounded in terms of n, delta, |S| and the degree. The
  subspaces themselves depend on S and are not controlled in height. A sequence whose S changes has no common
  subspace, so a single Z is impossible. This is exactly the standing caution. lit18 Prop A also records
  quantitative lower costs >= exp(s^{1/2-o(1)}).
* **Honest intermediate class.** Fixed K, S and b, with coefficient height up to rho(eps) log H. rho is
  ineffective and probably small, so the class stays subpower. There is no visible way to reach coefficient
  height (1/20) log H, where every block qualifies.

## (d) Literature

* **(VI.3)** is Evertse's estimate (Compositio 53 (1984) 225-244). In the two-variable form it is Corvaja-Zannier,
  Monatsh. Math. 144 (2005), Prop. 1 (2.2): sum_{mu in S} log^-|f(u,v)|_mu < -eps max h(u,v) only on finitely
  many translates of subtori. They call it "essentially a corollary of Evertse's paper" (read in
  `ref6L_meaning/dl/cz.pdf`, l.285-300).
* Astra's extension to coefficients of height o(L) is the "almost S-unit" setting of Xiao (arXiv 2110.01751),
  already in lit18 #7. Her proof via the projective Subspace Theorem H(P)^{-N-delta} (the CZ statement,
  l.241-250) is standard. I checked the zero-sum step: the log of the double product is -s h + o(L), against the
  required -(s-1+delta) h.
* **VI.1 itself** is not, to my search, an instance of a published theorem. CZ 2005 and Levin (Invent. 215 (2019))
  bound heights and gcds of f(u),g(u) at S-unit points of tori. VI.1 instead bounds the archimedean distance of
  Phi(u) ∈ T_5 to the rational point e, where Phi is a rational map from a torus of S-units to the
  anisotropic torus and the height is that of the image. The image height comes from primes outside S. The
  logarithmic-Jacobian step, with a Laurent polynomial F of negative total S-degree, is a Wronskian-type
  argument in the CZ tradition (cf. CZ, Acta Math. 193 (2004) on ||q alpha^n||). It is best described as **"an
  application of Evertse/CZ to a new target"**: new in this problem, built from standard parts. I did not search
  Autissier or the McKinnon-Roth follow-ups beyond the record. The record (expert note §3; papers_check) already
  states that no unconditional alpha^ess > tau_x at a rational point is known **for Zariski-dense sets of
  arbitrary points**. VI.1 does not contradict this, because its point set is thin (subpower).
* Cilleruelo-Córdoba and Cilleruelo-Granville are consistent: VI.6 is exactly the CC exponent (delta = 1/10,
  alpha = 10/3). Bourgain-Rudnick gives the log R bound, which is not relevant to five points.

## (e) Recommendation

**The integer-model and (*) sub-thread is exhausted for Question B.** (*) ⇔ A_5(some delta<1/10) is the
"beat Liouville at all" problem. Every unconditional lever found in Parts III-VI needs imposed structure: fixed
ingredients (VI.1), one affine parameter (VI.5), or fixed curves (VI.3). The remaining step Astra names, forcing
one chord scalar Im C_ij >= H^c off Z for unrestricted tuples, is a restatement of A_5 near the endpoint.

**One next question with a real, limited chance.** Extend VI.1 to *unbounded* products in one fixed quadratic
order (the resonant regime):

> (Q-res) Fix a quadratic unit lambda and a finite F. Five-tuples whose cut blocks are products of arbitrarily
> many factors a lambda^r + b lambda'^r (a,b ∈ F, r free), including cyclotomic quotients, satisfy
> Delta >= H^{-1/4-eps} off a proper closed Z.

* First step: rewrite each block as an integer polynomial in t = lambda^{g} (g = gcd of exponents) of degree
  ≍ log H with bounded coefficients, i.e. Laurent in one unit.
* Then try to replace VI.3 by a one-variable estimate. This is where it could fail: one variable, many terms,
  so Evertse fails.
* Test numerically against the record's cyclotomic five-row family first. If that family reaches arc
  R^{1/2-delta} Zariski-densely, (Q-res) is false and the family is a candidate counterexample to A_5. Either
  outcome is informative about whether resonance can refute Question B.

It does not touch the uniform bound for general clusters.

**Exact wording for the expert note (§3, after "Covering families are excluded"):**

> *Structured genuine families (round 31, Part VI; proved, ineffective).*
> (i) Let K ⊇ Q(i) be a number field, S a finite set of its places and b fixed. Consider five-tuples whose class
> is eta·prod_e (G_e/Ḡ_e)^{±} over cuts e, with each G_e a sum of at most b S-units of K, all coefficients and
> eta of height o(log H), and prod_e|G_e|^2 <= H^{1+o(1)}. For every eps>0 such tuples with
> Delta <= H^{-1/4-eps} lie in a proper closed Z = Z(K,S,b,eps). The exponent is 4, on all profiles. The parent
> primes are unrestricted; only the additive ingredients are fixed. This is the opposite of the fixed-prime-set
> (Ridout) results of §2. Consequently no Zariski-dense counterexample to A_5(delta) (any delta>0) or to
> (L_5°) comes from boundedly many unit-template blocks in finitely many fixed orders with subpower overlap.
> Resonant (overlapping, cyclotomic) templates and moving fields are not covered. The hypothesis cannot be
> relaxed to bounded degree: v_{n+1}+v_{n+1}^{-1}-v_n-v_n^{-1}-1 = 0 with v_n = (n+sqrt(n^2-4))/2.
> (ii) Ten-pair profile, G_e = a_e B + c_e with one Gaussian B and subpower coefficients: exponent >= 4 off a
> closed set of dimension <= 2. This is the five-point counterpart of the affine A_4 counterexample (Thm 4.1).
> (iii) A Zariski-dense genuine sequence on the ten-pair profile at Delta = H^{-3/10+o(1)} exists iff A_5(delta)
> fails for every delta<1/10. Its non-existence is therefore A_5 just beyond Cilleruelo-Córdoba. Question B
> implies it; it is independent of (L_5°) as far as is known.
> None of (i)-(iii) changes the growth rate (2/3+eps) log R/log log R, Question B, or M(C).

## Significance (blunt)

Growth rate, M(C) and Question B are all unchanged. VI.1 is the first theorem in the record that beats Liouville
(and reaches Vojta's exponent 4) on genuine five-tuples whose parent primes are unrestricted. The class is thin
(subpower) and defined by additive structure that no cluster is forced to have. Its value is negative-constructive:
it closes the natural route to refuting Question B by unit templates. VI.6 is a clean and correct bookkeeping
equivalence. It shows the (*) question is the hard core ("beat Liouville"), not a shortcut. Score for progress on
the goal: 3/10. For correctness and honesty of scope: 8/10.
