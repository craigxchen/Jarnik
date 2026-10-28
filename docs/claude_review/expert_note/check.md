# Adversarial check of `expert_note.md`

Checked against: `overview/overview.md` (§§1–2, 4, 5.4–5.5, 6.6, 7.1, 7.1A, 7.6), `growth/ptolemy.md` (§1, §6, §7),
`round18/README.md`, `round18/theory.md`, `round18/referee_theory.md`, `round19/README.md`, `round19/theory.md`,
`round20/README.md`, `round20/gauss.md`, `round20/kill_gauss.md`, `round21/README.md`, `two_thirds/proof.md`
(Thm 7.5, Cor 8.1), `round8/weakest.md` §4, `threshold/construct.md`. All paths are relative to `docs/claude_review/`.
Line numbers refer to `expert_note.md` as it stands now. The note was not edited.

Severity: **Major** means a false or overstated mathematical claim, a status overclaim, a wrong pointer for a
load-bearing claim, or a missing definition without which an outsider cannot work on the stated problem.
**Minor** means wording, missing caveats, missing quantifiers or constants, or a gap in self-containedness that a
reader can repair.

## Findings

1. **Major. §2 Data, line 81: "singletons and $[7]$ occur in no minor".** This is false for $[7]$. Since
   $[7]\ni i,j$ for all $1\le i<j\le7$, $n_{[7]}$ divides all 21 minors $\det(P_i,P_j)$ with $i,j\ge1$; only the
   singletons occur in no minor. The record's reason for ignoring $[7]$ is different: it can be inflated by a fresh
   prime ($z\mapsto nz$, i.e. $\mathrm{diag}(n,1)$) without moving $\Phi$, and it enters no Plücker relation.
   **Fix:** "Singletons occur in no minor. The top block $[7]$ divides every minor with $i,j\ge1$, but it can be
   inflated by a fresh prime without moving $\Phi$, so balance is imposed only on the 119 visible blocks, which is
   equivalent." Sources: `round18/README.md` §1 ("Balance" bullet); `round18/theory.md` §1.1; `overview/overview.md`
   §7.1 ("enter no Plücker relation").

2. **Major. §2, lines 102–103: "$(L_8^f)$ … So a slower-growing variant is not easier in substance."** This repeats
   the overstatement that referee item A withdrew ("'No easier than $(L_8)$' and 'the only freedom is in the
   constants' are therefore overstated"). What is proved: any proof of $(L_8^f)$ must prove $(L_8^c)^{\rm cat}$, which
   is the whole $H_\Sigma=0$ content of $(L_8)$; $(L_8^f)$ is not shown to be equivalent to $(L_8)$.
   **Fix:** "Every proof of $(L_8^f)$ must prove $(L_8^c)^{\rm cat}$ (the $H_\Sigma=0$ content of $(L_8)$), so a
   slowly growing $f$ is no shortcut; $(L_8^f)$ is not known to be equivalent to $(L_8)$." Sources:
   `round18/referee_theory.md` §2.3 and §2.9 item A; `round18/README.md` §2.3 item 3.

3. **Major. §5 item 9, lines 373–374: "Bounded-degree covering curve families with infinitely many points per member
   are excluded at $k=8$ and, for circles, at $k\ge5$."** As written this is false, and it is exactly the misreading
   that correction M1 removed. The $\mathbb Q$-lines of Kapranov's $\mathbb P^5$ form a bounded-degree covering
   family, and every member carries infinitely many rational points. What is excluded is families whose members
   *each carry infinitely many balanced small-$H$ points* (at $k=8$, genus $\le1$ by census Prop 7.5; genus $\ge2$ is
   Faltings), and, for circles, *infinitely many A-points* (gauss Thm 4.3(ii), $k\ge5$, $\delta>0$). §2 of the note
   (lines 156–158) states this correctly. §5 drops the qualifier.
   **Fix:** insert "balanced small-$H$" (for $k=8$) and "A-" (for circles) before "points per member". Sources:
   `round18/referee_census.md` §7 (M1, verbatim in `overview/overview.md` §7.1 item 4); `round18/README.md` §2.1
   item 4; `round20/gauss.md` Thm 4.3(ii); `round20/kill_gauss.md` §5 item 7.

4. **Major. §3 and §5, lines 239 and 358: the pointer `round5/construct.md` Cor C3 is wrong.** `round5/construct.md`
   has no Corollary C3. "Corollary C3 (the cube profile defeats every auxiliary polynomial)" and the mean-width
   Theorem B are in `threshold/construct.md`; `round20/gauss.md` §4.4 cites `threshold/construct.md` Theorems B, C2,
   C3. (Both `overview.md` and `round20/README.md` write only "`construct.md` Cor C3", which is ambiguous, and that is
   probably the source of the slip.) **Fix:** replace both occurrences with `threshold/construct.md` Cor C3 (with
   Thm B). The pointer row on line 395 (`round5/construct.md` for fakes) is correct and should stay.

5. **Major. §3, Question B, lines 256–257: the suggested refutation "by a rational 2-parameter family of balanced
   5-tuples at exponent $\ge7/15$, as at $k=4$".** Refuting $(L_5^\circ)$ (and hence $V_E(5)$) needs balanced
   bounded-residue points that are Zariski-dense in the 4-dimensional $T_5$, off every proper closed $Z$. At $k=4$ the
   2-parameter family plus the curve parameter $s$ gives a 3-dimensional dominant map onto the 3-dimensional $T_4$
   (Thm 4.1 needs this dominance, via the exact Jacobian rank 3). At $k=5$ a 2-parameter family plus $s$ has image of
   dimension $\le3<4$, so a proper closed $Z$ can absorb it, and it refutes nothing. The record also predicts that any
   such search fails: by the killer's dimension count, one-variable polynomial families cannot be dense at $k=5$ for
   any profile, and Thm 4.3(ii) excludes bounded-degree covering families. The phrase comes from `round20/README.md` §7,
   which is loose, but the note must not present it as a refutation route.
   **Fix:** "Or refute it by a Zariski-dense supply of balanced bounded-residue 5-tuples at exponent $\ge7/15$ (for a
   rational family, one dominating $T_5$). By the killer's count, one-variable polynomial-block families cannot do this
   at $k=5$ (**heuristic**), and bounded-degree covering families are excluded (Thm 4.3(ii))." Sources:
   `round20/gauss.md` Thm 4.1 proof (dominance); `round20/kill_gauss.md` §4 (dimension count) and §6.3;
   `round20/README.md` §7.

6. **Major. §3, the hypothesis "within one unit class" (lines 192, 199, 203) is never defined.** It is essential:
   the integrality of $t_{ab}$, and hence $H_\Sigma\ge0$ termwise, holds only inside one unit class (killer minor 1).
   An outside expert cannot work on $(L_5^\circ)$ without the definition. **Fix:** define it once in §1 or §3, for
   example: "for a fixed choice of the $\pi_p$, write $\zeta=\varepsilon\prod_p\pi_p^{a(p)}\bar\pi_p^{e_p-a(p)}$ with a
   unit $\varepsilon\in\{\pm1,\pm i\}$; two points lie in one Gaussian unit class when their units agree,
   equivalently $\zeta_a/\zeta_b=P/\bar P$ with $P\in\mathbb Z[i]$." Also give the alternative the record allows:
   $t_{ab}:=|c_{ab}|/2$ real, with conclusion $H_\Sigma\ge cw-O(1)$. Sources: `growth/ptolemy.md` §1
   ($z_i/z_0=P_i/\bar P_i$ within the retained class); `round20/kill_gauss.md` §1.2 item 2 and §7 item 1;
   `round20/README.md` §1.1 item 1; `overview/overview.md` §4.3. (`docs/claude_review/` has no formal definition; it
   lives in item 532 on the branch, so say so.)

7. **Minor. §3 opening, line 176: "Everything in this section is proved (killer-checked) or conditional."** The
   section also contains open statements ($A_k(\delta)$ for $k\ge5$, $(L_5^\circ)$, the arc form, Question B),
   heuristic statements (density count $e^{-w}$, "$\delta=0$ is counting-critical") and an un-literature-checked
   remark. **Fix:** "Every theorem in this section is proved (killer-checked) or conditional; the statements
   $A_5$, $(L_5^\circ)$ are open." Source: `round20/README.md` §1.1 Outcome.

8. **Minor (label consistency). §3, line 233, Lemma K1 "proved, killer-checked"; line 218, the $k=3$ corollary.**
   By the note's own legend (line 8–11), "killer-checked" means checked by a killer other than the author, and
   "unrefereed" means checked by its author only. Lemma K1 and the $k=3$ projection corollary are the killer's *own*
   results ("Nothing here is refereed", `kill_gauss.md` header), so no second party checked them. `overview.md` uses
   the same "killer-checked" label, so this is a legend problem, not a record conflict. **Fix:** label them "proved by
   the killer (`kill_gauss.md` §6.1; §3 `k3_projection.py`), not independently checked", or widen the legend.
   Source: `round20/kill_gauss.md` header and §7; `round20/README.md` §6 ("Lemma K1 … proved new content").

9. **Minor. §3, lines 202–205: the arc form is presented as what $(L_5^\circ)$ "reads as" at $k=5$.** The record keeps
   two forms. The residue form ($H_\Sigma\ge cw$) follows from the arc form by the argument of Prop 3.1, and the
   converse is not recorded as proved. **Fix:** "A sufficient arc form (it implies $(L_5^\circ)$ by the argument of
   `gauss.md` Prop 3.1; the converse is not recorded) reads: …". Also note that "$|t_{ab}|\ge1$" (line 205) needs the
   one-unit-class hypothesis. Source: `overview/overview.md` §7.1A (parenthesis after the Question);
   `round20/gauss.md` Prop 3.1 and §4.4.

10. **Minor. §4, lines 317–320: misattributed round-18 remark.** Round 18 called the gcd "a Vojta-type gcd problem in
    its own right" for pencils whose base locus **meets the interior**. These are outside $P^\circ$. For base locus in
    codimension-2 boundary strata (that is, $P^\circ$), round 18 asserted without proof that the excess gcd vanishes on
    $\mathrm{Bal}^{pc}$. Round 19 §4.1 calls that extension plausible but unproved and gives an unrefereed sketch of
    the obstruction. **Fix:** "Round 18 asserted, without proof, that for base locus in codimension-2 boundary strata
    the excess gcd vanishes on $\mathrm{Bal}^{pc}$ (`round18/theory.md` §2.4, after Prop 2.6); round 19 §4.1 records
    this as unproved. (For base locus meeting the interior, round 18 calls the gcd a Vojta-type gcd problem in its own
    right.)" Sources: `round18/theory.md` §2.4; `round19/theory.md` §4.1; `round19/README.md` §1.1 item 5.

11. **Minor. §4 Question C(c), lines 339–340.** The record already says Theorem 2.4 "holds verbatim for any movable
    class $\gamma$ with 31 replaced by $\max_k\Lambda_k\cdot\gamma$", and referee m7 notes that the maximum must be
    taken over orderings and that the remark is unused. **Fix:** cite it, and ask instead whether the inequality is
    known or sharp in that generality. Sources: `round19/theory.md` §7 item 4; `round19/README.md` §1.1 item 6 (m7).

12. **Minor. §5 item 6, lines 365–366: "This is vacuous, because the blocks are squarefree and pairwise coprime."**
    Actual blocks need not be squarefree. The proved statement is vacuity on squarefree profiles. Every recorded fake
    is squarefree, so abc through pair, triangle, Ptolemy or $n$-term relations, together with the fake-satisfied local
    inputs, cannot prove growth. **Fix:** "abc through relations among products of chords is vacuous on squarefree
    profiles (Prop 4.2), which contain every recorded fake (Cor 4.3), so it cannot prove growth." Source:
    `round8/weakest.md` Props 4.1–4.2, Cor 4.3 and "Precise scope"; `overview/overview.md` §5.4 row 15.

13. **Minor. §5 item 7, line 367: gcd normalisation is "invisible to every method above (Lemma K1)".** K1 covers
    forced-divisibility, auxiliary-form, Liouville, pair, cut, character and single-form inequalities. It does not
    cover items 2, 5 or 8, and the record says Vojta's $D=E$ term and the heuristic count do see the reduction.
    **Fix:** list the covered classes explicitly. Source: `round20/kill_gauss.md` §6.1; `overview/overview.md` §7.6
    item 10.

14. **Minor. §2, lines 94–96: "Roughly, $(L_8)$ says that a rational point … cannot have all 119 boundary contacts of
    the same size."** That is $(L_8^c)^{\rm all}$, which implies $(L_8)$ but is not known to follow from it. $(L_8)$
    gives it only for caterpillar points. The sandwich follows, but "roughly" may still mislead. **Fix:** "The
    strongest natural form, which implies $(L_8)$, says …; $(L_8)$ itself gives this only at caterpillar points."
    Also add the constants and the one-sided converse that the record states: $(L_8^c)^{\rm all}(\eta_1)\Rightarrow
    (L_8)$ with $\eta_0=\eta_1/2$, $c=\eta_1/4$; and $(L_8)\Rightarrow(L_8^c)$ only for $\Phi$ admissible from anchor 0
    (any anchor, with symmetrised $Z$). Sources: `round18/README.md` §2.1 item 1; `overview/overview.md` §7.1 item 2.

15. **Minor. §2 Question A, lines 171–172: the refutation clause lacks its quantifier.** Correction B reads
    "$\mathrm{Bal}^{\rm cat}(\eta)$ dense **for all $\eta$** $\Rightarrow$ not $(L_8)$ $\Rightarrow$ … $\Rightarrow$
    not $V_8$". **Fix:** "for every $\eta>0$ and $w_1$, a Zariski-dense set of balanced points with $w\ge w_1$ that are
    caterpillar up to $o(w)$ of residue". Source: `round18/referee_theory.md` §2.7 item B; `round18/theory.md` §3.1.

16. **Minor (self-containedness of $(L_8)$, $(L_8^{pc})$). §2.**
    * *Anchoring* is used (line 166: "the anchoring … vacuous") but never named. Say that column 0 is the anchor,
      normalised to $P_0=(1,0)^T$, and sent to $\infty$ in the Kapranov chart $\mathbb P^5=\mathbb P(\mathbb Q^7/\mathbb
      Q\mathbf 1)$. Also say that $Z$ depends on the labelling, so a statement with an arbitrary anchor must symmetrise
      $Z$ (census m1). Source: `round18/README.md` §1 "Anchor"; `round18/theory.md` §2.2.
    * *Caterpillar* is used (lines 99–100, 172) without definition. Add: "a compatible split system is a caterpillar if
      one side of each split can be chosen so that the chosen sides form a chain; equivalently, all internal edges lie
      on one path" (`round18/theory.md` §1.2). Add caveat H: the block system equals the $p$-adic tree of $\Phi$ only at
      primes dividing no residue (`round18/referee_theory.md` §2.1).
    * $\ell_T(p)$ in $d_T=\prod p^{\ell_T(p)}$ (line 93) is undefined. One clause suffices: "$\ell_T(p)$ = the number
      of levels $t$ at which $T$ is a residue class of $x\bmod p^t$" (`round18/theory.md` Lemma 2.2 proof;
      `census.md` Lemma 1.1).
    * Thm D (line 138): define $B=\sum_TD_T$ (the boundary) and say that $A=K+B$ is ample by Keel–McKernan
      (`growth/ptolemy.md` §6).
    * Line 111: the index product is over $3\le|T|\le6$, not "$|T|\ge3$" (which would include $[7]$); write the range.

17. **Minor (self-containedness of $(L_5^\circ)$, $V_E(5)$). §3.**
    * $E_A$ (line 187) is the toric ray divisor for the ray $\mathbf 1_A$ of the permutohedral fan in
      $\mathbb Z^k/\mathbb Z\mathbf 1$. Say so (`round20/gauss.md` §1.1).
    * $c_A$ and $A$ in $V_E$ (lines 224–227) are undefined: $A$ ample, $A=\pi^*L-\varepsilon'E$ with $L=\sum a_\rho
      E_\rho$, $a_\rho>0$, and $c_A=\max a_\rho$ (`round20/gauss.md` Thm 4.2 proof).
    * $\hat\gamma$ and $D_{PQ}$ (line 228) are undefined. $D_{PQ}=E$ (marks $1..k$ colliding away from $\pm i$ in the
      twisted $\overline M_{0,k+2}$), and $\hat\gamma$ is the balanced circle class of `threshold/L8.md` §6.1. Either
      define them or drop the sentence and keep "non-vacuous (`gauss.md` Thm 4.2 Remark 3)".
    * Thm 4.3(ii) (lines 229–231) omits its hypotheses $k\ge5$, $\delta>0$ and the definition $B=\pi^*\sum_\rho E_\rho$.
      Note also that Thm 4.3(i) was narrowed by the killer (`kill_gauss.md` §1.2 item 4).
    * The similarity torus is defined correctly. Add that $[\zeta]$ is invariant under $\zeta\mapsto\gamma\zeta$,
      $\gamma\in\mathbb Q(i)^*$, which is why $Z$ must live in $T_k$ (`gauss.md` §1.1 and Thm 3.4′ check (b)).

18. **Minor. §3 approximation reading, lines 252–254.** "The only predictions beyond Liouville are Vojta and
    McKinnon's conjecture" is also the killer's unchecked statement ("beyond Liouville/Seshadri … conjectural"). Only
    the last sentence carries the "not literature-checked" label. **Fix:** put the label on both sentences, and keep
    "Liouville/Seshadri". Source: `round20/kill_gauss.md` §6.2.

19. **Minor. §4, lines 303–306 (Thm 3.3).** $n_1$ is undefined: $n_1=\#\{T:\nu_T(B)>\mu_T\}$, i.e. the boundary
    divisors in the fibre over 1. Also add the sharper "$\le-16$ if $k(u)=4$". In the necessary conditions
    (line 326), $n(f)=n_0+n_1+n_\infty$ (boundary divisors in the three special fibres) is undefined and differs from
    $n'(f)$. Source: `round19/theory.md` §1 and Thm 3.3; `round19/README.md` §4 J1.

20. **Minor. §4, lines 331–333 (J2).** Add the judge's own caveat: "the radical also counts $p$ once, which might
    compensate; that is not analysed". J2 shows a gcd excess. It does not show that the nested-chain extension fails.
    Source: `round18/README.md` §5 R18-J2.

21. **Minor. §5 item 8, lines 368–369.** The genus-0 Arakelov reduction is to rows 2, 8 **and Lim** (Mordell–Weil
    packing $e^{O(\mathrm{rank})}$; bounded rank is equivalent to the uniform theorem). In the unlikely-intersection row,
    the Rémond-cone case of the $\omega(N')$ ceiling is **heuristic**. **Fix:** add both. Source:
    `overview/overview.md` §5.4 rows 18–19; `round20/README.md` §3 rows 18–19.

22. **Minor. §1, lines 59–60.** The residue-design caps carry the referee's scope limits: all primes of $N$ must exceed
    $X$, and the twisted extension at $X=M^A$, $A>1$, is open. The note gives these only by reference in §5 item 1.
    Add them here or point to §5. Source: `overview/overview.md` §5.2 item 2; `round5/README.md`.

23. **Minor. §1, lines 34 and 63–64: undefined symbols.** The cut identity uses $S_{p\tau}=\#\{i:a_{ip}<\tau\}$, the
    product over $p\mid N$, and the normalised cluster, none of which is defined. The window uses $P_+(N')$ (largest
    prime factor of the normalised norm $N'$). Define them in one line each. Source: `overview/overview.md` §2.2
    Step 2, §5.5.

## Checked and found accurate (no action)

* Theorem 1.1, with its constant and quantifiers (`two_thirds/proof.md` Thm 7.5, Cor 8.1). The lower bounds
  $M(C)\ge10$ ($C\ge8.2394$) and $M_\infty(C)\ge8$ ($C>5.89771$).
* The three Vojta routes and their varieties. "No $o(\log R/\log\log R)$ bound from any standard conjecture."
* $(L_8)$ data, chain condition (= Lemma 1.1(a)), 119 visible blocks, $H_\Sigma$ with $7+21$ terms, balance on
  visible blocks, quantifier order of $(L_8)$.
* The sandwich chain and "no arrow known to reverse". $K=-6H+\sum(|T|-2)E_T$, $H\cdot\beta=31$, $K\cdot\beta=10$,
  index $\prod d_T^{|T|-1}$, scales $e^{31w}$, $e^{49w}$, $e^{30w}$, the gap $e^{18w}$, heuristic $98-108=-10$.
* The lattice form of $(L_8^{pc})$, matching `round18/theory.md` §5(2) and `overview.md` §7.1 item 3, with item G.
* Thm D: $Z_{1/43}$, $\eta_0=1/381$, $h_K\ge10w-282\eta w-60H_\Sigma-O(1)$, $K\cdot\beta=2^{k-3}(k-8)+k+2$. Step 1
  constants $2\sqrt2\cdot255/\sqrt{M-7}$, $3584w/(M-1)$, $w=W/128$. Lemma 6.1 $c(Z)M^7$.
* Thm E numbers (24 primes, $59.2$, $10^{49.4}$). Prop 2.3 (finite part $\kappa$). M1 in §2. Lift Thm 4.2 scope.
  Ridout Thm 2.5. Correction B. Lemma 7.7 direction.
* $T_k$, Lemma 1.1 heights, $d_{ab}$, $t_{ab}$ and its cross-class correction, $A_k(\delta)$, $(L_k^\circ)$.
  Thm 3.4′ range and the convexity regime $\delta\le2^{-k-1}$. Converse not proved. Thm 4.1 and the $k=3$ corollary.
  $K_{Y_k}$, $\kappa=(k-1)(1+2\delta)/4-1$, $\kappa=2\delta$ at $k=5$, $\varepsilon<(1/32)/c_A$ via $A_5(1/64)$.
  $(K+D_{PQ})\cdot\hat\gamma=1/16$; $D=0$ needs $k\ge6$. The $7/15$ arithmetic. Liouville thresholds $1/30$ and
  $1/10$. The comparison with $(L_8)$.
* Gain definitions and identities. Prop 2.6 with $\lambda\in\mathbb Q^*$ (item F). Thm 2.4 with
  $\Lambda_k\cdot\beta=16,24,28,30,31$ and its equality cases. Cor 2.5. Lemma 3.2 ($N(u)\ge32$). Thm 3.3 and Cor 3.4.
  Exact $-32$ for $D\le6$. Bruno–Mella labelled. The four necessary conditions and the definition of $P^\circ$.
* Every cited file exists on branch `claude/upbeat-einstein-p0peju`. The exception is the wrong theorem pointer in
  item 4.
