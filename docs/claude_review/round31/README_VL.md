# Round 31, Part V (Question B thread) -- judge report: local isotropy at five marks

Source: Part V of `astra_lemma_V/answer.md` (heading "# Part V", §§V.1-V.4 and the twelve-line summary; scripts
`v_five_test.py`, `v_isotropic_star.py`, `v_star4_certificate.py`; question `reply5.txt`, Astra's summary
`final5.txt`). This is the Question B ("lemma") thread. Three referees each took one block:

| key | report | scope |
|---|---|---|
| full | `ref5L_full.md` | V.1 (factorisation, the circle formula (V.1), common multipliers, two-sided width), V.2 (size-free full isotropy) |
| partial | `ref5L_partial.md` | V.3 (IV.2 pencil fails isotropy, specialisation lemma, four-block star, flips, five-block curve) |
| exact | `ref5L_exact.md` | V.4 (rigidity at defect $P^{2-\varepsilon}$, the transfer (V.2), the nine-block remark), plus assessment |

Judge script: `judge_VL/check.py` (sympy). It independently rebuilds the four-star at the parametrised certificate
shape $(v,t,k)=(3/2,2,2)$ and the five-block curve, and checks the minor identities, the isotropic incidences, $d$
and the defect degrees. **Every claim labelled [J] is the judge's own and is unrefereed.**

Labels as in Parts III-IV. "Proved" means a complete proof that a referee checked line by line, or an exact
computation that a referee reproduced with independent code. "Computed" means exact data reproduced independently.
"Evidence" means not certified. **No referee found a false statement or a gap in any claim Astra made.** The one
GAP verdict (exact (b')) concerns a question the referee posed ("what is known for 5 to 9 isotropic blocks"). It
records an open state that Astra herself declares open. I rule on it in §1.3. The two CORRECT WITH CORRECTIONS
verdicts on "the two-sided question is Question B itself" concern a wording in the coordinator's summary. Astra's
text is more careful (§3(b)).

## 0. Verdict in five lines

1. **Nothing changes for the integer problem.** The growth rate stays $(2/3+\varepsilon)\log R/\log\log R$. There is
   no uniform bound $M(C)$. Question B ($A_5(\delta)$ for some $\delta<1/30$, i.e. $\alpha^{\rm ess}_e>15/4$ off a
   proper closed set) is exactly where it was.
2. **Everything Astra claimed is correct.** The serious corrections are of scope. The V.3 isotropic families reach
   width $N^{-3/10}$ only with defects $\delta_i\asymp n_i^{7/8}$ (relative error $N^{-1/20}$; [J] confirmed, exact
   degree 7 against $\deg n_i=8$). Incidence counts depend on the orientation, so the intrinsic data are 4 isotropic
   blocks (a star, dense off any $Z$) and 5 blocks (one curve, no escape).
3. **Full isotropy is genuineness in disguise.** Without a size bound it is the genuine problem with the height
   decoupled from $N$, so it gives no width bound (V.2). With subpower multipliers it is the genuine problem at
   height $N^{1+o(1)}$, on the ten-pair-cut profile, at the Liouville/Cilleruelo-Córdoba endpoint $\alpha=10/3$.
   That slice is **not** Question B. Excluding it off $Z$ is strictly weaker than $A_5(\delta<1/30)$.
4. **The one new theorem with content is V.4 rigidity.** At defect $\le P^{2-\varepsilon}$, three isotropic
   endpoint blocks not forming a star force exact norms. The referee extended it to a peeling hierarchy
   $P^1/P^2/P^3$ for 1, $\ge5$ and $\ge8$ blocks. The nine-block remark is correct with exponent $(20-k)/3$. All of
   these characterise genuine tuples. None gives an inequality for them.
5. **Score for Part V on the goal: 1/10.** It is a correct, honest calibration of the isotropy-versus-size trade-off,
   with no progress on Question B. The full-isotropy line should be closed.

## 1. Referee-corrected statements per block

### 1.1 full: V.1 and V.2 (referee: no error; one overstatement; partly not new)

**Proved**
* **Propagation.** At $p\mid J_e$ the marked columns span one line mod $p$ (primitivity plus $p\mid g_{ij}\mid\det$).
  One isotropic member makes the whole line isotropic, and the marked dot products vanish. This holds prime by prime
  for squarefree composite blocks, in both orientations. Elementary; implicit in the record (free thread Part III).
* **Factorisation and local multiplier rule.** Under full isotropy, $V_i=G_iU_i$ with $|U_i|^2=c_i$. Here $c_i$ is
  $2^b$ ($b\le1$) times split primes, each in a single orientation, which must be $\pi_p$ on marked rows. The
  converse holds. **Not new**: `astra_free_IV/answer.md` l.866-868 ("full-isotropy endpoint", "exact star norms force
  isotropy"), refereed in `README_III.md`.
* **(V.1), corrected statement.** *$H=\prod_{p\ \rm split}p^{\max_is_{ip}-\min_is_{ip}}$ is the least norm of a
  primitive Gaussian 5-tuple with the given doubled ratios, and $N/C^2\le H\le NC^2$ for $C=\mathrm{lcm}\,c_i$; the
  points $B(G_i/\bar G_i)(C/c_i)U_i^2$, $B=\prod_{p\mid N}\bar\pi_p$, are integral of norm $NC^2$.* This was checked
  against Gaussian-gcd computations in 300 random models. The range formula is already in the record
  (`astra_free_III/iii_graph_notes.md` l.75). The sandwich and the explicit circle are new and elementary.
* **Common multiplier.** $H=N\prod_{p\in D}p^{2v_p(c)}$. A single-orientation common factor cancels, and $c=2$ is
  equivalent to $c=1$ on odd parents. Referee's addition, accepted: mixed $c_i\in\{1,2\}$ is the record's
  unit-class phenomenon. New, elementary.
* **$c^4\mid\prod_{i<j}|t_{ij}|$.** Five rows in two orientation classes have at least 4 same-orientation pairs. New,
  elementary.
* **Two-sided width forces $c_i=N^{o(1)}$** (equal-weight blocks, $1\le|t|\le N^{o(1)}$, full isotropy). The extremal
  pair has $|W_i||W_j|=N^{o(1)}$, and Cramer coefficients bound the other $W_k$. [J] I checked this. The upper half
  $\Delta\le N^{-3/10+o(1)}$ is automatic under full isotropy because $|W_i|\ge1$, so the hypothesis that does the
  work is the **lower** bound on $\Delta$, which is an archimedean size bound in disguise (exact (c)).
* **Converse.** If $c_i=N^{o(1)}$ and $\Delta\le N^{-3/10+o(1)}$, then all ten pair-cut weights are $W/10$ (from
  $S_{ij}\ge3W/5$, $\sum S_{ij}=6W$). Standard averaging.
* **Dictionary, corrected statement** (both full and exact): *a two-sided fully isotropic model with
  $c_i=N^{o(1)}$ is a genuine primitive 5-tuple $\zeta_i=\lambda V_i/\bar V_i$ at height $H=N^{1+o(1)}$. It has the
  same anchored ratio point, weight $W/10+o(W)$ on each of the ten pair cuts and $o(W)$ on the other cuts, and width
  $H^{-3/10+o(1)}$, i.e. distances $R''^{2/5+o(1)}$. These are A-points for every $\delta\le1/10-o(1)$, so
  $A_5(\delta)$ for any $\delta<1/10$ puts them in $Z$.* Excluding them is strictly weaker than Question B, and
  $(L_5^\circ)$ does not apply (the profile is unbalanced). The index-one remark is true and irrelevant to $A_5$.
* **V.2 (size-free full isotropy), as stated.** *For every proper closed $Z$ and every $\epsilon(N)>0$ there are
  arbitrarily large squarefree split parents with equal-weight blocks, fixed cofactors, index one, all 20 endpoint
  (or all 25 balanced) incidences isotropic, $y\notin Z$ and $0<\Delta<\epsilon(N)$.* The ingredients are a CRT lift
  to $\mathrm{SL}_2(\mathbb Z)$, Zariski density of $\Gamma(N)A$, a shear $(x,y)\mapsto(x+mNy,y)$ with finitely many
  bad $m$, and a cycle flip that preserves $t_{ij}$, $|V|^2/n$ and the doubled ratios. Reproduced on the referee's
  own 193-digit parent and on Astra's. The escape from $Z$ is inherited from V.3's star. New, but nearly trivial.
* **Why V.2 contradicts nothing** (referee's quantification, accepted). On the referee's instance, $H\approx10^{2104}$
  against $N\approx10^{193}$. Along the shear $H\sim m^{10}$ and $\Delta\sim m^{-2}$, so $\Delta\approx H^{-1/5}$, above
  $H^{-1/4}$: the tuple is not even an A-point at $\delta=0$.

**Computed.** The 191-digit parent (27 split primes): all identities, both CRT shears, the cycle flip, and all 20 and
all 25 incidence congruences. Pre-asymptotic: at $P=500D$ the "fixed" cofactors ($\sim10^{43}$) exceed the blocks
($\sim10^{19}$), and the base width is $10^{-13}$ against $N^{-3/10}=10^{-57}$ ($\log\Delta/\log N\approx-0.07$). The
instance also has 13 dividing both $t_{25}$ and $J_{25}$, contrary to the specialisation sieve. That is harmless
here: only $p\nmid t_{ak}$ for $a\in e$, $k\notin e$ is needed, and it holds. The instance certifies identities, not
exponents.

### 1.2 partial: V.3 (referee: no error; four scope corrections)

**Proved**
* **IV.2 pencil fails isotropy.** Mod $J_{ij}$, $|V_i|^2\equiv q^4\beta$ with
  $\beta=(-b_{ij}(h_i-b_{ij})+b_0x_i)^2+l^2(h_i-b_{ij})^2>0$, and for $q=P^{o(1)}$ this lies in $(0,J_{ij})$. Proof,
  not just the earlier numerical 0/20. This covers the IV.2 pencil only (not IV.3).
* **Specialisation lemma.** From five monic degree-4 Gaussian-rational columns with $[F_i,F_j]=d_{ij}J_{ij}$ and ten
  distinct monic quadratic blocks, two dilations plus an elementary squarefree sieve give infinitely many integral
  primitive models with odd squarefree split pairwise coprime blocks coprime to the fixed cofactors. Polynomial
  isotropy persists, and polynomial non-isotropy fails at the whole-block level for large $T$. No prime-tuple input
  is used. Cosmetic correction: add the primes of $\mathrm{Res}(J_e,\ |F_i|^2\bmod J_e)$ to the bad set to remove
  stray prime-level isotropy on non-isotropic blocks.
* **Four-star, corrected statement.** *For every proper closed $Z$ there is a rational shape $w$ with
  $s_3=s_1s_2$, $s_1s_4\ne0$ (chosen after $Z$) whose specialisations are endpoint-orientation models with ten
  fixed nonzero cofactors, squarefree split parents, index one, $\Delta\asymp N^{-3/10}$, ratio point outside $Z$,
  and exactly the four isotropic blocks $12,13,14,15$ (8 incidences). Here $\delta_1=0$ identically and
  $\delta_j\asymp n_j^{7/8}$ for $j\ge2$, so $|V_j|^2/n_j=1+O(N^{-1/20})$. The width constant is of the size of the
  fixed cofactor.* [J] Rebuilt at the parametrised shape $w=(1,-21/16,45/16,-153/32,-39/32)$:
  $s_3-s_1s_2=0$, $[G,A]=s_4$, all ten minor identities, exactly 8 isotropic incidences on blocks $12,\dots,15$,
  $\delta_1=0$, $\deg\delta_j=7$ for $j=2,\dots,5$ (`judge_VL/check.py`). The referee also checked $s_3=s_1s_2$
  identically on the $(v,t,k)$ chart and the Jacobian at three points. The quantifier order (shape and cofactors
  after $Z$; a diagonal sequence with growing subpower cofactors) is right.
* **Flips to a mixed orientation, corrected statement.** *Flipping two central blocks gives 10 isotropic incidences
  among the marked ones in a mixed orientation, with the same ratios, cofactors, chords and $|V|^2/n$. The
  intrinsic content is still 4 isotropic blocks, and the transfer bound stays $H\le N^{8/5}$.* So ten incidences
  carry exactly the information of eight. The $n_i$ then range over $N^{1/5}$ to $N^{3/5}$.
* **Incidence counts, corrected.** In the endpoint orientation an isotropic incidence forces its partner, so exact
  counts are even. "One" and "five" incidences mean "at least". Exact counts 2, 4 and 6 off closed sets were not
  constructed as exact counts.
* **Five-block curve.** $[G,B]=(x-1)^6$. $d=-418284581/1009043070$ is real. $L_2(i)=L_3(i)=0$. The change of
  parameter $x=1+1/T$ gives monic columns and rational sum-of-squares blocks. The isotropic blocks are exactly
  $12,13,14,15,23$ (10 incidences), with fixed cofactors, $|V|^2/n\to1$ (again $\delta_1=0$ and
  $\delta_j\asymp n^{7/8}$) and $\Delta\asymp N^{-3/10}$. [J] All of these identities rechecked independently,
  including that all ten minors are nonzero constants times $J_{jk}$ in the $T$-chart and that the defect degrees
  are $7$. **One shape, one rational curve, no escape from closed sets** (Astra says so). $A_5$ may contain that
  curve in $Z$.

**Computed.** Astra's three certificates rerun and pass. Two new lemma-faithful instances exist (a shape
$(1,-9,-5,7,8)$ with 823/883-digit parents, and the five-block curve with parents of over 3000 digits): exact
incidences, gcd 1 at the non-isotropic incidences, and $\Delta N^{3/10}$ constant, of the size of the fixed cofactor
($\sim10^{232}$, $\sim10^{959}$). All instances are pre-asymptotic in the width constant.

### 1.3 exact: V.4 (referee: no error; extension; one GAP on a referee question)

**Proved**
* **V.4 rigidity, statement as verified (quoted from `ref5L_exact.md` §1.1).**
  > *Let $J_e$ be pairwise coprime with $J_e=P^{1+o(1)}$, in the endpoint orientation. Let $V_1,\dots,V_5\in\mathbb Z^2$
  > be primitive, with $J_{ij}\mid[V_i,V_j]=t_{ij}J_{ij}$, $0<|t_{ij}|\le P^{o(1)}$ and $\delta_i=|V_i|^2-n_i$. Suppose
  > $\max|\delta_i|\le P^{2-\varepsilon}$, and that the set $S$ of isotropic edges contains three edges not all through
  > one mark. Then for $P$ large, $\delta_1=\dots=\delta_5=0$. In particular this holds whenever $|S|\ge5$.*

  The proof uses IV.4 item 1 only: the linear identity (IV.5), the existence of $Q_*$, and the integrality of $DE$,
  valid up to $P^{3-\varepsilon}$. It does not use the rank-one dichotomy, the width, index one, squarefreeness or
  primality. [J] I read the proof (answer.md l.1292-1298) and the referee's line-by-line check. The step
  $\gcd(J_{jk},[V_a,V_j])\mid t_{aj}$ uses that $J_{jk}$ is coprime to $J_{aj}$, which holds because the blocks are
  private and pairwise coprime. Correct.
* **Isotropy kills the Part IV rank-one branch** (referee, accepted). If $E=\lambda\ell^2$ and two isotropic edges
  meet at $a$, then $\ell\propto[V_a,\cdot]$ and $|\delta_j|\ge P^{2-o(1)}$. So below $P^{2-\varepsilon}$ the
  $\mathbb Z[\sqrt k]$-cluster branch of `README_IVL.md` §3(c) is empty once two isotropic edges meet, and below
  $P^{1-\varepsilon}$ once there is one incidence.
* **Peeling extension (referee's own theorem).** *At defect $\le P^{s-\varepsilon}$, $s\le3$: exact norms follow if
  there are marks $a,j,k$ with $\deg_S(a)\ge s$, $\deg_{S-a}(j)\ge s-1$, $\deg_{S-\{a,j\}}(k)\ge s-2$.* At $s=3$ the
  condition holds for 125/210 six-edge sets (including $K_4$), 110/120 seven-edge sets (the failures have degrees
  $(2,2,2,4,4)$), and every set of $\ge8$ edges. Calibration only.
* **Transfer (V.2), corrected constant.** *With $k$ known isotropic endpoint blocks and $|V_i|^2\le(1+o(1))n_i$,
  $H\le(1+o(1))N\prod_{e\notin E}J_e=N^{(20-k)/10+o(1)}$. On equal-weight blocks at width $N^{-3/10+o(1)}$ the genuine
  tuples have $\alpha\le(20-k)/3$ and are A-points for $\delta\le(k-8)/(2(20-k))$.* [J] Rechecked with the record's
  normalisation (A-point: $\mathrm{dist}\le H^{-(1/4+\delta/2)}$): $1/\alpha=3/(20-k)$ and
  $\delta=2(3/(20-k)-1/4)=(k-8)/(2(20-k))$. For $k=9$, $\delta=1/22>1/30$; for $k=8$, $\delta=0$. **Nine is the
  threshold.** A Zariski-dense nine-block family would refute the **conjecture** $A_5(\delta)$ for every
  $\delta<1/22$ (hence every $\delta<1/30$), and would contradict Vojta's prediction $\alpha^{\rm ess}\ge4$. It
  contradicts no theorem.

**Computed** (referee). An exhaustive toy over $\mathbb F_{11}[T]$ and $\mathbb F_{13}[T]$ (split form $XY$, constant
cofactors, $\deg\delta\le2$) found 0 propagation failures and 0 violations, and found 480 non-exact five-block nonstar
models at $\deg\delta=2$ exactly. So the exponent 2 is sharp in the toy. Caveat (referee's): no toy tuple reached
$\deg\delta\le1$, so the predictions T2 and T3 are tested only vacuously. On the V.3 families the hypothesis fails by
$P^{3/2}$: defects are $P^{7/2}$, even the IV.5 identity fails (residual of degree 6), and the 191-digit instance has
$\log|\delta|/\log P\in[3.85,3.91]$, pre-asymptotic for $3.5$. This is consistent with the theorem: the five-block
curve contains the nonstar triangle $12,13,23$, so it **must** have defect $\ge P^{2-o(1)}$.

**Ruling on exact (b') GAP.** This is the referee's question, not a claim of Astra. Astra states that she has "no
unconditional upper bound of eight" and "no essentially dense critical family with five whole endpoint blocks"
(EVIDENCE paragraph). The referee's list of open items is accurate, and I adopt it as the "Open" column of the
table in §3(c). **No gap in any claim of Astra.**

### 1.4 Withdrawn or corrected in Part V (collected)

1. Coordinator's summary, "the two-sided question is Question B itself", and Astra's summary item 4, "bounded or
   polylogarithmic multipliers thus return to Gaussian Question B". **Withdrawn as worded.** The corrected form is
   the dictionary in §1.1: the ten-pair-profile, Liouville-endpoint slice, which $A_5(\delta<1/10)$ covers. Astra's
   own V.1 sentence ("the applicable statement is the Gaussian $A_5$ ... with the very same ratio point and
   exceptional set") is correct as a statement of coverage.
2. "$|V_i|^2/n_i\to1$" in V.3 holds only at rate $N^{-1/20}$ (defect $n^{7/8}$). Every statement about the V.3
   families must carry this rate.
3. Incidence counts: report isotropic **blocks**. "Ten incidences" in a mixed orientation equals 4 blocks. Endpoint
   counts are even. The counts 1 and 5 mean "at least".
4. The five-block family is one curve and does not escape closed sets.
5. "Width $\asymp N^{-3/10}$": the constant is the size of the fixed cofactor. All numerical instances are far from
   the critical regime and certify identities only.
6. The 191-digit instance is not an output of the specialisation lemma (13 divides both $t_{25}$ and $J_{25}$).
   Harmless.
7. The specialisation lemma's bad set must include the resultant primes of the non-isotropic remainders (cosmetic).

## 2. What is new relative to the record, and correct

Grep: the record has the full-isotropy endpoint and "exact star norms force isotropy" (free thread Part III), the
range formula (`astra_free_III`), "bounded mean log slack forces genuine points on a positive fraction of marks"
(`README_III.md` l.171), and non-isotropy of the III.11 and IV.2 models (numerical 0/25, 0/20). It has no theorem
in which **partial** isotropy forces anything, and no isotropic model at critical width.

* **V.4 rigidity** (new): partial isotropy (a nonstar triple) plus defect $\le P^{2-\varepsilon}$ forces exact norms.
  This is the first statement of its kind. **Peeling hierarchy** (referee, new): $P^1/P^2/P^3$ for one incidence /
  nonstar triple / 3-2-1 sets ($\ge8$ blocks, $K_4$).
* **Rank-one branch closed under isotropy** (referee, new): it closes the Part IV residual branch whenever two
  isotropic edges meet, at defect $\le P^{2-\varepsilon}$.
* **Four-star family** (new): 4 isotropic blocks at critical width, dense off any $Z$, defect $n^{7/8}$.
* **Five-block curve** (new): one curve only.
* **IV.2 pencil non-isotropy, proved** (new; previously numerical).
* **Transfer (V.2) and the nine-block threshold** (new): $\alpha\le(20-k)/3$, threshold $k=9$ under $A_5$.
* **Elementary new items**: the sandwich $N/C^2\le H\le NC^2$ with an explicit circle; $c=2\Leftrightarrow c=1$;
  $c^4\mid\prod|t_{ij}|$; two-sided width $\Rightarrow c_i=N^{o(1)}$; V.2's shear construction; the specialisation lemma.

Not new: the factorisation $V_i=G_iU_i$, the local multiplier rule, the range formula, propagation, and Liouville
averaging.

## 3. Rulings

### 3(a) Isotropy without size

**Ruling: isotropy without size is worthless as a lever, and this is now a theorem (V.2).** Full marked isotropy (all
20 endpoint or all 25 balanced incidences), fixed cofactors, index one and squarefree split parents allow
**arbitrarily small** width off any proper closed set. The mechanism is that $\Gamma(N)$ preserves every local
datum (lines mod $p$, minors, primitivity), and a shear inside $\Gamma(N)$ sends all columns to infinity in one
direction. What is violated is $|V_i|^2\sim n_i$: the multipliers $c_i$ become huge ($\sim10^{382}$-$10^{768}$ on the
referee's instance), and the genuine height $H$ of (V.1) decouples from $N$. Measured against $H$, the width is only
$\approx H^{-1/5}$, not an A-point at all. So "a width lower bound in terms of the parent alone" is false, and any
statement that uses isotropy must also use the archimedean size $|V_i|^2/n_i$.

### 3(b) Full isotropy with two-sided size: is it an exact reformulation of Question B?

**No.** Here is the exact reformulation.
* (Forward, proved.) Full isotropy plus two-sided width $N^{-3/10\pm o(1)}$ on equal-weight blocks gives
  $c_i=N^{o(1)}$. By (V.1) this gives a genuine primitive Gaussian 5-tuple at height $H=N^{1+o(1)}$, with the same
  ratio point, a ten-pair-cut profile (weight $W/10$ each, $o(W)$ elsewhere) and width $H^{-3/10+o(1)}$.
* (Backward, [J], resting on the record's dictionary `astra_free_IV/answer.md` l.866-868 and the unit-class setup of
  the expert note §3, unrefereed in this direction.) A genuine primitive 5-tuple with a ten-pair-cut profile at
  that width, within one unit class, gives exact stars $G_i$ with $c_i=1$, which are fully isotropic with two-sided
  width. Mixed classes give $c_i\in\{1,2\}$.
* **So the two-sided fully isotropic problem is, up to $N^{o(1)}$, exactly: do Zariski-dense genuine 5-tuples exist
  on the ten-pair-cut profile at the Liouville/Cilleruelo-Córdoba endpoint ($\delta=1/10-o(1)$, $\alpha=10/3$)?**
  [J] Existence would refute $A_5(\delta)$ for every $\delta<1/10$, hence the $A_5$ route to Question B. It would
  not refute $(L_5^\circ)$, since the profile is unbalanced. Non-existence off $Z$ is implied by $A_5(\delta)$ for
  any $\delta<1/10$, which is much weaker than Question B ($\delta<1/30$, the window $(15/4,4)$, balanced
  fifteen-cut profiles). The integer-matrix language adds nothing: with subpower multipliers it is literally the
  Gaussian factorisation. Kill criterion (i) of the referees applies.

### 3(c) Partial isotropy: blocks versus existence / forced exactness / open

Setting: endpoint orientation, equal-weight blocks $J_e=P^{1+o(1)}$, $P=N^{1/10}$, subpower cofactors, width
$N^{-3/10+o(1)}$, $|V_i|^2/n_i\to1$. "Dense" means off every fixed proper closed $Z$. $k$ counts isotropic
**blocks**. Transfer: $H\le N^{(20-k)/10}$, $\alpha\le(20-k)/3$.

| $k$ | existence (best known defect) | forced exact norms (proved) | transfer $\alpha\le$ | open |
|---|---|---|---|---|
| 0 | dense (IV.2, $P^{1+o(1)}$, mixed signs) | below $P^{1-\varepsilon}$: exact or rank one (IV.4) | $20/3$ | $\theta\le1/4$ |
| 1 | dense as a sub-block of the four-star ($P^{7/2}$); "exactly one" not constructed | below $P^{1-\varepsilon}$ (peeling $s=1$) | $19/3$ | exact-count family; defect in $[P,P^{7/2})$ |
| 2-4, star $K_{1,k}$ | dense (four-star, $P^{7/2}$, $\delta_{\rm centre}=0$) | below $P^{1-\varepsilon}$ only; rank one impossible below $P^{2-\varepsilon}$ | $18/3$-$16/3$ | split residual $E=[V_a,X]L/D$ below $P^{2-\varepsilon}$: neither built nor excluded |
| 2-4, nonstar (two disjoint edges, path, triangle) | triangle $12,13,23$ only inside the five-block curve (one curve) | nonstar triple ($k\ge3$): below $P^{2-\varepsilon}$; two disjoint edges: below $P^{1-\varepsilon}$ only | $17/3$-$16/3$ | any dense nonstar family, at any defect |
| 5 | one rational curve ($K_{1,4}+23$, $P^{7/2}$); no escape | below $P^{2-\varepsilon}$ (V.4) | $5$ | dense family at any defect; any model with defect in $[P^2,P^3)$ |
| 6-7 | none | below $P^{3-\varepsilon}$ for 3-2-1 sets ($K_4$; 125/210 and 110/120 of the sets) | $14/3$, $13/3$ | existence at all |
| 8 | none | below $P^{3-\varepsilon}$ (all) | $4$ ($\delta=0$; no contradiction) | existence at all |
| 9 | none | below $P^{3-\varepsilon}$ | $11/3$ ($\delta=1/22$): **dense would refute $A_5$** | existence; unconditional exclusion |
| 10 | full isotropy: forces $c_i=1$, exact, genuine (any defect $o(n)$) | always | $10/3$ = the genuine ten-pair problem (§3(b)) | = §3(b) |

Reading [J]: isotropy and norm accuracy are traded against each other along the table. The only constructions put
the isotropy on a star, and they pay with defect $n^{7/8}$, above the whole range $P^{3}$ in which any forcing
theorem lives. Every "forced exact" entry turns a model into a genuine ten-pair tuple at $\alpha=10/3$, whose dense
existence is itself open (§3(b)). So no entry excludes anything unconditionally.

### 3(d) What changes for Question B, the growth rate and the uniform bound (blunt)

* **Growth rate:** unchanged, $(2/3+\varepsilon)\log R/\log\log R$.
* **Uniform bound $M(C)$:** unchanged.
* **Question B:** unchanged. Part V answers Part IV's open question 1 ("isotropic integral models"). Up to 4 blocks
  are free on a dense set, 5 are free on a curve, 10 is genuineness, and size-free isotropy is worthless. That
  closes the hope that local isotropy by itself is a lever. The levers proved in Parts III-V are all of the form
  "size $\times$ isotropy $\Rightarrow$ genuine": exact norms; full isotropy with $|V|^2<2n$ or $c_i=N^{o(1)}$; defect
  $P^{1-\varepsilon}$ with one incidence; $P^{2-\varepsilon}$ with a nonstar triple; $P^{3-\varepsilon}$ with a 3-2-1 set;
  plus the one-way transfer with 9 blocks. **None of them gives an inequality for genuine tuples.** A Question B
  proof must still beat Liouville at $e$ for genuine balanced points, and nothing in Part V touches that.
* **The integer-matrix programme for Question B** (Parts III-V of this thread) has now mapped the fake side
  completely enough. Fakes fail exactly at "exact norms" or "isotropy with accurate size", and both of these collapse
  to genuineness. I recommend closing the programme as a route to Question B. It has become a classification of
  fakes, not an approach to the genuine inequality.

## 4. Recommended updates to the expert note (recommendations only; the record is read-only)

1. §5 item 12, after the "Part IV, Question B thread: the norm defect" bullet, add:
   > *Part V, Question B thread: local isotropy* (`round31/README_VL.md`). Call a block $J_e$ isotropic when
   > $J_e\mid|V_i|^2$ at its marked columns (this propagates along the marked line). (i) Without a size bound, full
   > isotropy with fixed cofactors and index one allows arbitrarily small width off any closed set ($\Gamma(N)$
   > shear); the genuine height $H$ decouples from $N$. (ii) Full isotropy with two-sided width $N^{-3/10\pm o(1)}$
   > forces subpower multipliers and is the genuine problem on the ten-pair-cut profile at the Cilleruelo-Córdoba
   > endpoint $\alpha=10/3$, covered by $A_5(\delta)$ for any $\delta<1/10$ and strictly weaker than Question B.
   > (iii) In the endpoint orientation with $|V|^2/n\to1$ (defect $n^{7/8}$), four isotropic blocks forming a star
   > occur densely, and five occur on one curve. (iv) At defect $\le P^{2-\varepsilon}$ ($P=N^{1/10}$) a nonstar
   > triple of isotropic blocks forces exact norms, and at $\le P^{3-\varepsilon}$ so do $\ge8$ blocks. (v) With $k$
   > isotropic blocks the model yields genuine tuples of exponent $\le(20-k)/3$; a dense nine-block family would
   > refute $A_5$. Nothing here bears on genuine tuples.
2. §3 (Question B), under the approximation reading: add "the ten-pair-cut slice at $\alpha=10/3$ (Cilleruelo-Córdoba
   endpoint): whether it is Zariski dense is open; non-density is implied by $A_5(\delta<1/10)$ and is strictly
   weaker than Question B".
3. Closed routes: add "local isotropy without archimedean size (V.2)" and "full isotropy with subpower multipliers as
   a new formulation (it is the Gaussian factorisation)".
4. Wherever isotropic models are cited: count isotropic **blocks**, not incidences, and state the defect rate.

## 5. Open sub-questions, ranked by honest probability of progress on Question B

All probabilities are the judge's estimates [J].

1. **Density of the ten-pair-cut Liouville slice.** Are genuine primitive 5-tuples with a ten-pair-cut profile and
   width $H^{-3/10+o(1)}$ Zariski dense in $T_5$? A "no" would be the first unconditional beat of Liouville on a
   dense set at $e$ for one profile (still weaker than Question B). A "yes" would refute the $A_5$ route. P ≈ 0.01.
   *First step:* extend the record's template-density bound (ten-pair template curves of fixed block degree sweep
   out dimension $\le2$) to block degree growing with the parameter, or find the first degree at which dimension 3
   is reached.
2. **Isotropic star at small defect.** Construct or exclude the split residual $E=[V_a,X]L(X)/D$ with 2-4 isotropic
   edges at $a$ and defect below $P^{2-\varepsilon}$. This fixes the minimal isotropy count at the $P^2$ level.
   P ≈ 0.003 (calibration only). *First step:* impose $DE=[V_a,X]L$ with symbolic $L$ on the four-star ansatz
   (degree-4 columns over $\mathbb Q(i)[T]$) and solve for $\deg\delta_j\le3$.
3. **A dense nonstar family** (a triangle or three edges without a common vertex), at any defect. Decide whether
   V.4's coefficient-matching conditions ($Bb_1=A+1/C-15$, $Db_1=A/C+C-15$, $h^2-b_1h-1/C=0$, with the $z_j$
   consistent) have Zariski-dense rational solutions. P ≈ 0.003. *First step:* compute the dimension and the
   rational parametrisability of that variety with a Gröbner basis over $\mathbb Q$.
4. **The nine-block threshold.** Prove unconditionally that no dense family has $\ge9$ isotropic blocks with
   $|V|^2/n\to1$. This would be an unconditional consequence of a conjecture's prediction, but by §3(c) it is
   essentially equivalent to excluding dense genuine tuples at $\alpha\le11/3$, i.e. as hard as an $A_5$-type
   statement. P ≈ 0.002. *First step:* none recommended before item 1.

## 6. Scores (0-10, on the goal)

| block | score | reason |
|---|---|---|
| full (V.1-V.2) | 0.5 | correct; partly already in the record; size-free isotropy and two-sided full isotropy both reduce to genuineness; the "Question B itself" wording withdrawn |
| partial (V.3) | 1 | correct, new constructions (dense four-star, one five-block curve); defect $n^{7/8}$; incidence counts overread |
| exact (V.4) | 1.5 | the one clean new theorem (partial isotropy plus defect $P^{2-\varepsilon}$ forces exactness), sharp in a toy, extended to $P^3$; the nine-block threshold is exact; characterisations only |
| **overall Part V** | **1** | correct and well calibrated; closes local isotropy as a stand-alone lever; no movement on Question B, the growth rate or the uniform bound |

No web sources were needed. Part V cites no external theorem beyond Gauss's lemma, CRT and Zariski density of
congruence subgroups (standard). The specialisation sieve is elementary, and no Tao-Ziegler input is used.
