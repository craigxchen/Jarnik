> **Round 31 has thirteen judge reports.** This file is the judge's report on the first pass (Part I). Later passes:
> `README_II.md`, `README_III.md`, `README_IV.md` (open thread, Part IV: the orbifold-BMY inequality (IV.2), the
> uniform bound in the function-field model, the characteristic-2 construction), `README_IVL.md` (Question B
> thread, Part IV: the norm defect), `README_V.md` (open thread, Part V: unbounded clusters over the fixed field
> $\mathbb F_2$ with an anisotropic norm; filter F3+), `README_VL.md` (Question B thread, Part V: local
> isotropy), `README_VI.md` (open thread, Part VI: unbounded clusters over $\mathbb F_3(t)$ with $A^2+B^2$;
> squarefree parents over $\mathbb F_2(t)$; the parity-code condition; filter F3++), `README_VII.md` (open thread,
> Part VII: unbounded clusters with every prime exponent at most two, also in a ramified model; exactly squarefree
> parents in odd characteristic open; filter F3+++) and `README_VIL.md` (Question B
> thread, Part VI: the Subspace theorem on fixed-field unit-sum blocks; $A_4(\delta)$ false for $\delta<1/6$; the
> ten-pair slice is equivalent to $A_5(\delta<1/10)$); `README_VIIL.md` (Question B thread, Part VII: five marks
> over a characteristic-0 function field; the sparse theorem; dense polynomial Pell families; the unrestricted
> question open); `README_VIII.md` (open thread, Part VIII, back on integer circles: the rank theorem
> $M\le1+\omega(N')+K_C$, a height-order inequality, angle-model fakes inside the low-height window); `README_VIIIL.md` (Question B
> thread, Part VIII: $(L_8)$ implies Rudin's conjecture and a divisor bound below Lenstra's exponent $1/4$; round
> 8's fixed-threshold Subspace hypothesis is false; conditional growth from relative-height counts). Referee reports: `ref_*.md` (I), `ref2_*.md` (II),
> `ref3_*.md` (III), `ref4_*.md` and `ref4L_*.md` (IV), `ref5_*.md` and `ref5L_*.md` (V), `ref6_*.md` and
> `ref6L_*.md` (VI), `ref7_*.md` and `ref7L_*.md` (VII; `ref7L_alpha_*.md` are three further checks of the bound
> $\alpha^{\rm ess}_e(-K)\le4$), `ref8_*.md` and `ref8L_*.md` (VIII). `langer_check.md` records the check of Langer's Theorem 0.1 in the original. The
> coordinator's messages are `astra_brief.txt`, `astra_*/reply*.txt` and `astra_*_reply4.txt`; the colleague's
> cumulative reports are `astra_free_VI/answer.md` and `astra_lemma_VI/answer.md` (Parts I-VI).

# Round 31 -- judge report: an outside audit (Astra) of the record

Sources: `astra_free/answer.md` (open-ended audit, §§1-7) and `astra_lemma/answer.md` (attack on Question B, §§1-6),
both frozen, with the brief `astra_brief.txt`. Seven referees each took one block: `ref_cover.md`, `ref_baker.md`
(astra_lemma §§1-2), `ref_square.md` (astra_lemma §§3-4, 6 and the parity/angle-model items), `ref_local.md`
(astra_lemma §5), `ref_divisor.md` (astra_free §§1-2, 4), `ref_sieve.md` (astra_free §§3, 5, 6) and `ref_ff.md`
(astra_free §7). Judge scripts are in `judge/`. **Every claim labelled [J] is the judge's own and is unrefereed.**

Scope. The directory also contains a second exchange (`astra_free_II/`, `astra_lemma_II/`, `ref2_grid.md`,
`ref2_template.md`). It is **not** covered by this report, which judges only the two frozen reports named above and
the seven referee reports listed.

Labels. "Proved" means a complete proof checked line by line by a referee, or an exact or exhaustive computation that
a referee reproduced with independent code. "Computed" means exact data reproduced independently. "Evidence" means a
search that is not exhaustive. No referee reported FALSE. No referee reported GAP on a headline claim. The judge
therefore had no disputed claim to rule on. The scope corrections below are the referees', and the judge read the
relevant passages of Astra's reports to confirm them.

## 0. Verdict in five lines

1. **Nothing changes for the integer problem.** The growth rate stays $(2/3+\varepsilon)\log R/\log\log R$. There is
   no uniform bound, and Questions A and B and route 5 keep their status. Astra says the same of her own work.
2. **The strongest new result is a proof of $A_5(1/31)$ off a proper closed $Z$, for every full-rank primitive
   five-tuple whose valuation lattice has index $I\ge7$** (astra_lemma §2). Two ingredients make it work: McKinnon–Roth
   on ramified isogeny covers, and a two-logarithm Baker–Wüstholz bound for large index. Two referees found no error.
   It is peripheral: generic tuples, the F1 fake profiles and all explicit families have index 1. [J] checked this on
   600 random Paley-profile 5-subtuples, all of index 1.
3. **Second-strongest: the sieve rigidity theorem** (astra_free §3). The all-height, completely multiplicative
   graft condition (round-20 (**)) forces $\phi_p\equiv\pm\arg\pi_p\pmod{\pi/4}$ (via Ratner), and $B_0=4$ is sharp.
   So item (b) of overview §5.5 moves from "heuristic collapse" to "proved collapse". Item (a) and the prime-only
   sharp sieve stay open.
4. **Two closures:**
   * the round-24 family $M_C$ has **no** nondegenerate rational member in either $\mathbb Q$-form, and each form is
     excluded at 2 and at 3 (Q24.2 settled negatively);
   * at $D=8$, $c=0$, squarefree $N$, the function-field model has **rank $\le4$**, and $M(8,0)=8$ with a unique
     extremal code (the round-30 follow-up, done).
   Also new: a converse to round 8's "integral reduced chords $\Rightarrow$ genuine", and an F1-obstructed divisor
   reformulation of growth.
5. **Barrier statements.** Astra's literal-slip point about the expert note's opening is correct and should be
   fixed. Her F2 point is a wording matter: the record already states route 5's supercriticality. Four further
   sentences of the expert note are now out of date (§3.3). **Round score on the goal: 2/10.**

## 1. Referee-corrected statements per block

### 1.1 cover: astra_lemma §1 and §2.1 (referee: no error; scope corrected)

**Proved**
* **§1 conventions.** For a primitive five-tuple of norm $H=R^2$, put $y_j=z_j/z_5$. Only split primes occur in $H$,
  and $y_j=\epsilon_j\prod_p g_p^{b_{jp}}$ with $g_p=\pi_p/\bar\pi_p$. The index $I=[\mathbb Z^4:L]$ of the column
  lattice is intrinsic and equals the gcd of the $4\times4$ minors. Further, $h(s)=\tfrac12\log D(s)$ and
  $|s-1|\ge D(s)^{-1/2}$. The $A_5(\delta)$ arc condition is $\Delta\le H^{-(1+2\delta)/4}$. The common phase is
  removed by ratios, not set to zero. (Routine; equivalent to `round20/gauss.md` Lemma 1.1.)
* $B^4=70$ for the $A_4$ root polytope, and the toric height is $\log H$. (Already round22 `toric.md` Prop 1.1.)
* **Rational lift.** If the columns lie in $L=P\mathbb Z^4$, then $x_k=\prod_p g_p^{h_{kp}}$ with $h=P^{-1}b$ is
  rational and $y=\epsilon f_P(x)$, with at most $4^4$ unit translates. No gcd or denominator cost is dropped. (New.)
* **Fixed-cover theorem, corrected wording.** $X_f$ is the twisted toric variety of the same fan in the finer lattice
  $N'=P^{-1}N$, with polytope $P^T R_5$ of normalised volume $70I$, so that $B_f^4=70I$. All fibre points of
  $x^P=\epsilon^{-1}$ lie in $(S^1)^4$; they are real roots-of-unity points. "Properness" is to be read as follows:
  *$f_P$ on $(S^1)^4$ is an $I$-sheeted covering map, so a real lift of $y$ near $e$ lies within $C_P\Delta$ of a
  unique fibre point.* Heights upstairs and downstairs are comparable, with constants depending only on
  $(P,\epsilon)$.
* **McKinnon–Roth as used.** Cor 4.2 (p. 27) gives $\beta_q(L)\ge\frac n{n+1}(\mathrm{Vol}(L)/\mathrm{mult}_qX)^{1/n}$
  on an irreducible $X$, which may be singular. The judge checked this against `papers_check/read_mr.md` lines 131-152.
  Thm 5.1 allows $x_v\in X(\bar k)$, an ample $\mathbb Q$-bundle and $S=\{\infty\}$. The finiteness used is (5.2),
  applied once per fibre target, giving a finite union. $S$ is fixed, so the growing-$S$ caution does not apply.
* For each bounded range $[d_0,I_*]$ of indices, $A_5$ at exponent $g<B_{d_0}=\frac45(70d_0)^{1/4}$ holds off a proper
  closed $Z$. Reason: there are finitely many HNF lattices with $256$ values of $\epsilon$ each, and the image of a
  proper closed set under the finite map $\bar\psi$ is proper closed. Rank 4 is essential.

**Gap:** none. **Scope correction:** the uniform statement "for all $I\ge7$" also needs §2.2 (block 1.2).

### 1.2 baker: astra_lemma §2.2, §2.3 and the conclusion of §2 (referee: correct; one citation caveat)

**Proved**
* **Inequality (2).** For every full-rank primitive tuple, with absolute and effective constants,
  $$\frac{-\log\Delta}{W}\ll I^{-1/4}\log(2I)+\frac{\log(2+W)}{W}.$$
  The proof takes a short dual vector $v\in L^*$ with $\|v\|_\infty\le I^{-1/4}$ (Minkowski; Astra's factor 2 is
  harmless). Its order $m$ divides the exponent of $L^*/\mathbb Z^4$, so $m\le I$. Then $u=\prod g_p^{\langle v,b_p\rangle}$
  is non-torsion with $2h(u)\le4\|v\|_\infty W$, and $(y^c)^4=u^{4m}$ with $|u^{4m}-1|\le4\|c\|_1\Delta$. Finally,
  Baker–Wüstholz in two logarithms ($u,-1$; degree 2) gives $-\log|u^n-1|\ll(1+h(u))\log 2n$. The bounds $I\le24(W/\log5)^4$
  and $\|c\|_1\le8I$ finish it.
  * *Citation caveat (referee):* the BW specialisation was checked through a secondary source (arXiv:1611.04664) and
    against the product-of-heights form in `round8/lw.md`, **not** against the Crelle full text. The judge accepts it.
    It is the standard two-logarithm case, and any Baker-type bound of shape $\ll h(u)\log n$ suffices.
  * By-product: $\Delta\le H^{-\kappa}$ forces $I\le I(\kappa)$, effectively.
* **Theorem (astra_lemma §2), as corrected.** For $d_0\ge1$ and $0<g<B_{d_0}$ there are a proper closed $Z\subset T$
  and $H_0$, depending only on $(d_0,g)$, with the following property. Every primitive five-tuple with
  $\mathrm{rank}\,L=4$, $I\ge d_0$, $H\ge H_0$ and $y\notin Z$ has $\Delta>H^{-1/g}$.
  * The threshold $I_*(g)$ is effective; $Z$ and $H_0$ are not.
  * Numerics: $B_5=3.4602>17/5>10/3$ and $B_6=3.6216<15/4$. Also $B_7=3.7639>124/33=3.7576$; [J] recomputed, with a
    margin of only $0.006$.
  * **Hence $A_5(1/31)$ holds for full-rank tuples with $I\ge7$, and the approximation exponent exceeds the Liouville
    value $10/3$ for $I\ge5$.** Index 7 is sharp for the method.
* **§2.3, deleting a fixed prime set $S_0$.** The theorem remains true with $I$ computed from the columns outside
  $S_0$; outside rank $<4$ counts as index $\infty$. The constants depend on $S_0$, as stated. Cases A, B and C were
  all checked exactly. The example with full index 1 and outside index 7 is realised.

**Gap:** none.

### 1.3 square: astra_lemma §§3, 4, 6, and the angle-model and identity items (referee: correct; one hypothesis added)

**Proved**
* **(4)/(5).** If $y_1,\dots,y_4$ are dependent in $U(\mathbb Q)/U(\mathbb Q)^2$ (torsion bit included), then
  $\Delta\ge2^{-1/2}H^{-1/4}$. Units are arbitrary and there is no exceptional set. This is exponent 4.
* **(6), balanced squareclass case:** exponent $5(1-\eta)/(1+\eta)$, with $H_\Sigma\ge(1-7\eta)w-O(1)$.
* **Reduction:** the one-unit-class $(L_5^\circ)$ needs only the indices $I\in\{1,3,5\}$. Even $I$ forces mod-2
  dependence. This is **conditional on §2 with $d_0=7$** (blocks 1.1-1.2, now refereed), so it now stands.
* **Zariski-dense balanced full-rank family in the squareclass class**, with an added hypothesis: *the singleton
  primes are fixed, or satisfy $\log p=o(\sqrt w)$, with even exponents of order $w/\log p$*. This shows (6) is not
  vacuous, and nothing more.
* **(7):** $\min d_{ij}\le(4W+2w_E)/7$. This gives exponent $7/2$ for squarefree norms with $\le9$ primes, and 4 for
  $\le6$ primes. Both constants are optimal for pair averaging (exact LP).
* **Ten-pair-cut angle model.** This is an index-1, Smith-form $(1,1,1,1)$ fake satisfying every monomial gap (axial
  and diagonal), with $\Delta\asymp H^{-3/10}$, and it avoids any proper closed $Z$. "Sharp" means sharp for
  gap-only inputs. (A routine instance of round8 Thm 6.1 on the round23 2|3 profile.)
* **Identity (8):** $\sum_iP_i/\prod_{j\ne i}t_{ij}=0$ on ten-pair-cut one-unit-class tuples. Checked exactly on 200
  tuples. No bound follows from it.
* The "exact loss" statements of §6 and the fifteen-cut index-1 profiles: trivial restatements.

**Gap:** none. All bounds in this block are consequences of monomial gaps, so fakes satisfy them (Lemma K1 class).

### 1.4 local: astra_lemma §5 (referee: theorem correct; closing remark corrected)

**Proved**
* The algebra of round24 `search.md` §5.1 ($c,D,u,v,G$, $v/u=4(m^2-1)/(m^2-9)$, the formula for $u(m)-u(-m)$) is
  reproduced in sympy. The exact square-class identity is
  $\mathrm{disc}\,G=u(2D)\,[2(m+1)c/((m-3)^5(m-1)^3)]^2$.
* **Descent.** A nondegenerate ray-generic member over $\mathbb Q$ forces $u\in\mathbb Q^{*2}$, rational $m$, and a
  split $G$. *Correction:* the step that makes $x_0$ rational must state why the 16 projected values are distinct (a
  coincidence would put three places in a degree-2 fibre). Alternatively, sum the two Galois-stable subsets
  separately. Averaging all 30 place images gives the wrong value.
* **Corrected theorem: neither $\mathbb Q$-form of $M_C$ has a nondegenerate ray-generic member; each form is excluded
  independently at 2 and at 3.**
  * Twisted form at 2: $A=-(m+3)c$ and $B=2D$ cannot both be squares (table mod 256).
  * Split form at 3: shown 3-adically.
  * The extra exclusions found by the referee: twisted at 3, and split at 2, where $m\equiv5\pmod8$ needs precision
    beyond $2^8$.
  * Astra's sentence "the two forms have different local obstructions" is **withdrawn**: there is no split/twisted
    contrast, and so no anisotropy signal.

### 1.5 divisor: astra_free §§1, 2, 4 (referee: all correct; significance F1-obstructed)

**Proved**
* **Reduction (1)-(2).** If growth fails for some fixed $C$, there are fixed $B,b>0$ and unbounded normalised $N$ with
  $D_B(N)\ge b(W/\log W)^2$. Here $D_B$ counts the divisors $16\sqrt N\le d\le\sqrt NW^B$ with
  $d-\lfloor\sqrt d\rfloor^2=t^2$, $(s,t)=1$ and $16t^2\sqrt N\le d$. There are also $\gg(W/\log W)^3$ square-product
  triples. Both the subdivision and the one-unit-class normalisation are needed, and the referee exhibited data
  showing each one is.
* The average bound $\sum_{X}^{2X}D_B\ll X^{3/4}\log\log X$ holds. It is pointwise-useless.
* The nonprimitive counterexample $N=L^2(x^2+1)$ checks out, so coprimality is essential.
* The pointwise core restriction holds. It gives no saving at squarefree parents.
* **Triangle theorem.** Any three edge divisors with separately chosen primitive near-real representations satisfy
  $\psi_i=\psi_j+\psi_k$ for some labelling. This is a single Liouville gap, which fakes also satisfy.
* **Five-vertex reconstruction** ($m\ge5$, or a unique diameter): primitive near-real representations on all edges
  give a genuine cluster at norm $N$ realising the profile up to column reversal and shift. There is an imprimitive
  variant. The four-vertex exception at $N=131490079981$ is exact.
* **Scalar reconstruction (§4).** Integer $t_{ij}$ on all edges incident with a straddling set $U$ ($|U|\ge2$,
  $m\ge3$) give genuine norm-$N$ points with these phase differences. **This is the first proof of the converse that
  `round8/referee_stability.md` §5.4 flagged as unproved.** A non-straddling counterexample was reproduced. At the
  barrier scale $O(m\log(2r)+W+m)$ tested edges suffice.
* §1 logic (scope of fake tests; conjunction and quantifier remarks): correct and routine.

**Computed (reproduced exactly)**
* The census of odd split-supported $N\le10^7$ at $B=4$: 18180/113/2 parents, and zero square-product triples.
* The 22-prime squarefree census: max 6 candidates; 404 incidences, all phase-closed, all of parallelogram shape
  $0,v_1,v_2,v_1+v_2$.
* The 13-prime prime-power census: 53/42, all closed.

**Referee answer adopted.** For squarefree $N$, $D_B(N)$ equals an angle-model count of signed divisor forms within
$W^{B/2}/4$ of their Liouville limit. Paley fakes have $\gg m^2$ distinct difference vectors meeting it. **So the
target $D_B=o((W/\log W)^2)$ is F1-obstructed.** It is also stronger than growth, and may be false.

### 1.6 sieve: astra_free §§3, 5, 6 (referee: correct; scope corrected)

**Proved**
* **Rigidity theorem.** Let the angles at $p\in S$ be grafted, extended completely additively to primitive ideals, and
  assume they satisfy (3) (= round-20 (**)) at **all** heights for primitive $\alpha,\beta$ with norms coprime to $S$.
  Then $\phi_p\equiv\pm\arg\pi_p\pmod{\pi/4}$, and conversely those angles satisfy (3).
  * Proof steps checked: Ratner's orbit-closure theorem on $(\mathrm{SL}_2(\mathbb R)/\Gamma(q))^2$ (the only
    intermediate algebras are $\Delta$ and the whole algebra); the rational-direction cost $\sqrt{h/Q}$ (attained, by
    the referee's independent search); and triangularity of the signed-valuation matrix, giving $Q_{0,p}=p$.
* **Operator bound.** An all-primitive-ideal dual large sieve $AX+B_0K$ with $B_0<4$ forces (3).
  **Sharpness:** the irrational one-prime graft at $p=13$ ($\phi=\arg(2+3i)/2$) passes both monomial-gap families
  and satisfies the sieve with $B_0=4$, yet fails (3). The $p=29$ three-fibre test is an illustration only.
* **Corrected scope statement.**
  * Item (b) of overview §5.5 (the multiplicative product-scale gap) becomes: *proved; it forces genuine angles mod
    $\pi/4$, so it is the problem itself (as composite Prop 3.1); $B_0=4$ is attained by an irrational graft;
    ineffective in height.*
  * It does **not** touch item (a) (no multiplicativity) or the prime-only sharp $(X+K)$ sieve.
  * It says nothing about a sieve checked only up to a prescribed height, or about template-dependent thresholds.
  * It gives no growth.
* **(5) suffices for growth, and $(L_5^\circ)\Rightarrow(5)$.** This is a routine variant of round23 `split.md` Thm 1.
  No reverse implication is known. Formally (5) is weaker; it is not shown to be easier.
* **§6.** For primary generators, $A_4(\delta)$ fails for $\delta<1/14$ and $(L_4^\circ)$ fails in one unit class. The
  Jacobian is checked exactly over $\mathbb Q(i)$ ($\det=1674i/1349553725$) and mod 101. This re-proves the record's
  Theorem 4.5 with a new family, and repairs the unit-class point the record had judged harmless.

### 1.7 ff: astra_free §7 (referee: correct; scope stated)

**Proved**
* The round-30 census is reproduced independently: 64 six-row rank-5 orbits, of which **21** survive Lemma 1.1. They
  are identical to the judge's `classes_D8_c0.pkl`. The referee's Lemma-1.1 test is complete (LP plus a bounded DFS),
  so **the linear tool is exhausted at $D=8$**.
* New nonlinear tools are correct: the quartic identity, the proportional-row lemma, the four-collinear-row (Rolle)
  lemma, and the unit-ideal and positive-product certificates.
* **All 21 classes are excluded.** The referee reproduced this in four gauges with Groebner bases over $\mathbb Q$.
* **Theorem, corrected scope:** *for squarefree $N$ over $\mathbb C(t)$ (hence $\mathbb Q(i)(t)$), $D=8$, $c=0$, every
  cluster has rank $\le4$ and size $\le8$; the 8-row sign code is unique up to symmetry and is attained by the cube
  example.*
  * The theorem does not cover $\mathbb F_q(t)$ or non-squarefree $N$.
  * $M(8,0)=8$ is now a theorem. It rests on pairwise balance plus the rank exclusion.
  * The judge's $M(7,0)\le5$ and $M(9,0)\le5$ (round 30, [J]) **remain unrefereed**.

### 1.8 Withdrawn or corrected in the round
* "The two $\mathbb Q$-forms of $M_C$ have different local obstructions" (astra_lemma §5): withdrawn.
* The density step of the squareclass family: the hypothesis $\log p=o(\sqrt w)$ is added.
* The uniformity of §2.1 is now stated as resting on §2.2.
* The $x_0$-rationality step of astra_lemma §5 is repaired.
* Brief paraphrase "every determinant critical": the overstatement is in the brief, not the record (§3.3).

## 2. What is genuinely new relative to the record, and correct

The judge checked these against the record by grep.

1. **Ramified isogeny covers in McKinnon–Roth** (astra_lemma §2.1). The record had only dismissed étale covers (MR
   Cor 8.9; `papers_check/read_mr.md` line 185). Restricting to the thin set of points that lift raises the volume to
   $70I$.
   * Combined with the **two-logarithm dual-vector Baker bound** (2), this gives the first unconditional essential
     exponent above $15/4$ at $e$ on a Zariski-dense set of points.
   * The set is thin: it is the union over index $\ge7$ of translated cover images.
   * The record's literature statement is about $\alpha^{\rm ess}$ on *all* points off $Z$. That statement is not
     contradicted; it is refined.
2. **The parity bounds (4)-(6)** on squareclass-dependent tuples, and the bookkeeping that reduces the remaining
   indices to $\{1,3,5\}$. These are new, but Liouville-class.
3. **Sieve rigidity** (astra_free §3), with sharp $B_0=4$. This is new: Ratner and commensurability appear nowhere in
   the record. It closes item (b) of overview §5.5.
4. **Emptiness of both $\mathbb Q$-forms of $M_C$.** This upgrades round 24's ineffective Faltings finiteness, and
   Q24.2 is closed.
5. **Rank $\le4$ at $D=8$ in the function-field model**, and the uniqueness of the extremal code. Round-30
   follow-up §7 is done.
6. **Scalar and five-vertex reconstruction theorems.** These prove the round-8 converse "integral reduced chords
   $\Leftrightarrow$ genuine" that was flagged as unproved. Also new: the triangle theorem, and the parent-only
   divisor reformulation $D_B$. All of these are equivalences or fake-satisfied facts, not bounding inputs.
7. **Exact data.** The divisor censuses; the four-vertex exceptional parent $131490079981$; the $p=13$ graft; and a
   new primary-generator $A_4$ counterexample family.

Routine or already in the record:
* $B^4=70$; the §1 height conventions;
* $(L_5^\circ)\Rightarrow(5)$ and (5) $\Rightarrow$ growth (a variant of split Thm 1);
* §6 re-proving Thm 4.5;
* the ten-pair-cut angle model (round8 Thm 6.1 on round23's profile);
* the saturated-balance census (round 30).

## 3. What changes for the integer problem

### 3.1 Growth rate, uniform bound, routes
* **Growth rate:** unchanged at $(2/3+\varepsilon)\log R/\log\log R$ (Theorem 1.1).
* **Uniform bound $M(C)<\infty$:** unchanged and open. No lower bound for $M_\infty(C)$ moves.
* **Question A ($L_8$):** untouched. No block addresses it.
* **Question B ($A_5$, $(L_5^\circ)$):** status **open**. What it gains:
  * The exceptional-set statement $A_5(1/31)$ holds off a proper $Z$ for full-rank tuples with $I\ge7$.
  * Exponent 4 holds for squareclass-dependent tuples with no $Z$, and exponent 5 when they are also balanced.
  * The one-unit-class $(L_5^\circ)$ is reduced to $I\in\{1,3,5\}$.
  * **Index 1 is the whole problem.** It is the generic case: the referee's synthetic share of $I>1$ falls to 0% at
    14-15 primes, every explicit family has $I=1$, and [J] all 600 sampled Paley 5-subtuples have $I=1$ (`judge/idx.py`).
    At $I=1$ the cover has degree 1 and McKinnon–Roth gives only $2.31$.
* **Route 5:** untouched. No block addresses its zero estimate ZE$^{\rm all}_H$.

### 3.2 F1 and F2
* **F1 (fake circles).** No result of the round escapes F1 as a growth input. Sorted by how each relates to F1:
  * *Satisfied by fakes:* the divisor target $D_B$ (referee: the fakes realise $\gg(W/\log W)^2$), the triangle
    theorem, and the parity bounds.
  * *Characterisations of genuineness, not bounds:* the reconstruction theorems and sieve rigidity.
  * *Restricted to index $\ge7$:* the cover theorem, which uses actual algebraic points (McKinnon–Roth, Baker). F1
    fake profiles have index 1 [J], so it cannot reach them.
  * Astra's logical remark is correct: a model that satisfies two lists separately does not cap their conjunction.
    The record's main caps (Paley fakes, `round8/lw.md` Thm 6.1) are joint constructions, so the remark leaves them
    standing.
* **F2 (criticality).** Astra says "every determinant is critical" is contradicted by route 5. The record does not
  make that unqualified claim. Expert note line 387 and §5 item 3 (lines 535-538) already say that route 5's LP is
  supercritical and stops at a zero estimate. Lemma K1 (line 274) concerns inequalities, and route 5 produces vanishing,
  not an inequality. **Ruling: no substantive overstatement.** A clarifying qualifier is still worthwhile (below).

### 3.3 Sentences of `expert_note/expert_note.md` that should be corrected

Recommendations only; the record is read-only for this judge.

1. **Lines 40-42**: "No bound $o(\log R/\log\log R)$ is known for every $R$, either unconditionally or from any
   standard conjecture. The only standard-conjecture routes give the uniform bound directly."
   This is a literal contradiction (Astra, confirmed by the divisor referee), since a uniform bound is an $o$-bound.
   Replace with: *"No bound $o(\log R/\log\log R)$ is known unconditionally. Conditionally, the only known routes go
   through Vojta's Main Conjecture and give the uniform bound directly; no intermediate growth bound is known from a
   weaker standard conjecture."*
2. **Lines 68-70**: "P3-all$^*$, the sharp-constant large sieves and the low-height window … remain uncapped".
   Add: *"The completely multiplicative, all-ideal, all-height form of the sharp sieve (round-20 (**)) is not
   uncapped but collapses: it forces genuine angles mod $\pi/4$, with $B_0=4$ attained by an irrational graft (round 31,
   ineffective in height). The prime-only and non-multiplicative forms remain uncapped."*
3. **Line 274** ("Lemma K1: every recorded method is exactly critical") and **line 532**. Add the qualifier: *"every
   recorded method that yields an inequality at the actual point; route 5's forced vanishing is supercritical (item 3
   of §5) but yields no inequality without its open zero estimate."*
4. **Lines 298-306 and 533** ("no unconditional theorem with $\alpha^{\rm ess}>\tau_x$ … was found"; McKinnon–Roth
   "weaker than Liouville"). Add: *"Restricted to the thin, Zariski-dense set of points whose valuation lattice has
   index $I\ge5$ (resp. $\ge7$), McKinnon–Roth on the degree-$I$ isogeny cover plus a two-logarithm Baker bound gives
   exponent $>10/3$ (resp. $>15/4$, i.e. $A_5(1/31)$) off a proper closed set (round 31). At index 1 nothing changes."*
5. **Lines 371-373 and 409-410** ("Only finitely many curves of its complex component are defined over $\mathbb Q$
   (… Faltings; ineffective …)", "only finitely many over $\mathbb Q$ on the known component"). Replace with: *"No
   nondegenerate ray-generic member of the known component is defined over $\mathbb Q$, in either form; each form is
   excluded at 2 and at 3 (round 31). Other components (Q24.1) are open."*
6. **Lines 633-634** ("$M(8,0)=8$ … (**judge's computation, unrefereed**) … what excludes them, if anything, is
   non-linear"). Replace with: *"At $D=8$ all 21 classes are excluded by non-linear certificates (round 31, refereed):
   rank $\le4$ and $M(8,0)=8$ with a unique extremal code, for squarefree $N$ over $\mathbb C(t)$. $M(7,0)\le5$ and
   $M(9,0)\le5$ remain unrefereed."*

## 4. Open sub-questions created by the round, ranked by honest probability of a growth-rate improvement

All probabilities are the judge's estimates [J].

1. **Index 1 in Question B with an input that fakes violate.** The cover method has degree 1 there, and parity is
   void, so a genuinely new idea is needed. One concrete test: does identity (8) constrain the actual Gaussian
   arguments beyond the monomial gaps? P(growth improvement) ≈ 0.01.
2. **The prime-only sharp sieve (item (a)) by a Ratner-type rigidity argument**, together with an effective version
   (a height of the first violation in terms of the graft). This is the only analytic form left that matters.
   P ≈ 0.01.
3. **Refute or prove the primitive divisor count $D_B=o((W/\log W)^2)$.** It is F1-obstructed, so only a proof using
   actual integrality could help. A refutation, which is likely feasible, would cost nothing for growth. P ≈ 0.005.
4. **Indices 3 and 5.** Combine the cover ($B_3=3.05$, $B_5=3.46$) with parity or a second cover, to finish the
   reduction to index 1. This is thin-set bookkeeping. P ≈ 0.002.
5. **Other components of the 15/4 Pell class (Q24.1).** P ≈ 0.001.
6. **The function-field pipeline at $D=12$, and refereeing $M(7,0)$ and $M(9,0)$.** P < 0.001. The line stays closed.

## 5. Scores (0-10)

| block | score | reason |
|---|---|---|
| cover (lemma §§1, 2.1) | 4/10 | Correct and new: McKinnon–Roth on ramified covers is the round's best idea. Peripheral: it acts on the index $\ge7$ thin set, while the core of Question B is index 1. |
| baker (lemma §§2.2-2.3) | 4/10 | A correct, effective two-logarithm device that makes the index theorem uniform. Standard technique, and it relies on a secondary-source BW citation. |
| square (lemma §§3, 4, 6) | 2/10 | All correct. The bounds are Liouville-class on thin congruence classes, and the reduction to $\{1,3,5\}$ leaves the generic case. |
| local (lemma §5) | 3/10 | Correct and new: closes Q24.2 in both forms. The record had already rated this component irrelevant to Question B. |
| divisor (free §§1, 2, 4) | 3/10 | Clean reformulation, plus the first proof of the round-8 converse. Its proposed target is F1-obstructed. |
| sieve (free §§3, 5, 6) | 4/10 | A genuine new theorem (Ratner rigidity, sharp $B_0=4$) that corrects the barrier map. It closes a form that was not the relevant one. |
| ff (free §7) | 3/10 | Completes the round-30 follow-up correctly (rank $\le4$, $M(8,0)=8$). It does not transfer to integers. |
| **round, on the goal** | **2/10** | A competent, error-free outside audit. Several closures and one new technique, all peripheral. No change to the growth rate, the uniform bound, Questions A or B, or route 5. |
