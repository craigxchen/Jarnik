# Kill report: arakelov.md (round 20)

Target: `round20/arakelov.md` (status: reduced to caps). Scripts and outputs (independent code; nothing of the
proposal's `arakelov_checks/` is imported): `round20/kill_arakelov_checks/` = `kgauss.py`, `kclusters.py`,
`kdict.py` → `kdict_out.txt`, `kadelic.py` → `kadelic_out.txt`, `kmulti.py` → `kmulti_out.txt`, `k2adic.py` →
`k2adic_out.txt`. Notes: `round20/kill_arakelov_notes.md`.

**Verdict: KILLED as a lead.** No growth improvement and no uniform bound are claimed, and nothing survives that
could become one.
* Every input the lens uses lies in a capped class: checklist item 1 (row 2), item 3 (row 8), and item 5 / Lim.
* The residual "needed statement" is not a new statement. Its precise form ("equivalently ... $(L_8)$" or Vojta $D=0$
  on $\mathrm{Bl}_YX_7$ / $\overline M_{0,8}$) is **at least as strong as the target**: $(L_8)$ alone gives the uniform
  theorem (overview §4.2), and $V_8$ / Vojta on $\mathrm{Bl}_YX_7$ give it conditionally (overview Thms 4.1, 4.3). Its
  vague form ("global coherence, weaker than lattice gaps at every prime, not fixed-$S$") names no inequality.

The **mathematics of the write-up largely survives re-derivation and recomputation**: Theorem A is exact on all 26
configurations, the §4.2 table is reproduced to three decimals, and the §4.3 per-row counts are reproduced. But:
* it contains two arithmetic slips (§3);
* there is a real gap in the proof of Prop B.2: the 2-adic place is omitted, and the true 2-adic local minimum is
  about $2\times$ what the proof assumes (§2.3);
* Prop B.2 has a scope gap for $M\le X<2M$;
* it cuts at $M$ where it should cut at $X$;
* "equivalently" is overclaimed (§4).

None of these opens a lane; each one, once corrected, only tightens the reduction to caps.

## 0. What the proposal claims

1. **Theorem A** (dictionary, labelled proved and exact). On the regular model $\mathcal X_N$ of $uv=Nt^2$ over
   $\mathbb Z[i]$:
   * components equal levels;
   * $(P_i\cdot P_j)_{\rm fin}=\log\mathrm{Norm}\,c_{ij}$;
   * the vertical correction of $D_0$ is $Q_p$;
   * the cut identity holds.
2. **Theorem B** (proved). The $(2/3+\varepsilon)$ finite inequality is the adelic $M$-point (Fekete) inequality
   $\sum_v\epsilon_v\le0$.
3. **Lemmas 2.1–2.3 and Prop B.2** (proved). Genus-0 Arakelov inequalities lie in the class ADC, and ADC is capped
   by the row-2 fakes.
4. **Prop 3.3(b)** (proved). Zhang, Northcott and arithmetic bigness in the $(H,E)$-plane would need
   $\tau^*(M)>2$.
5. **§6.1** (evidence). The multi-point vector $\Gamma^{>M}_T$ ($|T|\ge3$) vanishes on the record clusters and
   on fake completions.
6. **Theorem C** (proved). Heights on $\overline M_{0,M}$ obtained from local data and product formulas reduce to
   pair product formulas, so pseudo-adelic completions of fakes satisfy them.
7. **Needed statement.** Global coherence, "equivalently" $(L_8)$ / Vojta $D=0$.

## 1. Recomputation of the dictionary (A)

**Re-derivation.** I re-derived each part by hand. Each holds for normalised clusters with $N$ odd and only split
primes, which is exactly what the proposal assumes.

* **(1)** $v_\pi(z_j)=a_{jp}$. The chart $V_{a+1}$ has $\alpha=z/\pi^a$, a unit, and $\beta=\bar z/\pi^{e-a-1}$, which
  vanishes mod $\pi$.
* **(2)** Different levels give $v_\pi(z_i-z_j)=\min(a_i,a_j)=v_\pi(g_{ij})$. Equal levels give
  $i_\pi=v_\pi(\alpha_i-\alpha_j)$. At $\varpi\nmid N$, including $(1+i)$ because $N$ is odd and the $(u,v)$-model is
  smooth there, $g_{ij}$ is a unit.
* **(3)** Write $y_\tau=x_\tau-x_{\tau-1}$. The chain system gives $y_\tau=-S_\tau(c)$. Summation by parts on the
  path Laplacian gives $-\Phi^2=\sum y_\tau^2$. The fibre $\sum C_k$ is in the kernel, so the normalisation $x_0=0$ is
  harmless. At $\bar\pi$ one has $S'_\tau=M/2-S_{e-\tau+1}$, so the sum is the same.
* **(4)**
  $\sum_{i<j}[\min(a_i,a_j)+\min(e-a_i,e-a_j)]=e\binom M2-\sum|a_i-a_j|=e\frac{M(M-2)}4+Q_p$,
  using $\sum|a_i-a_j|=\sum_\tau S_\tau(M-S_\tau)=eM^2/4-Q_p$.

**Exact check with my own code** (`kdict.py`). The configurations are rebuilt from the same inputs: the referee's
10-cluster and its two 9-subclusters, the five record clusters, nine best windows, and nine random subsets (my random
subsets differ from the proposal's because the enumeration order differs).
* (1): the component is tested from chart valuations.
* (2): $i_\pi$ is computed from the chart coordinate ($z/\pi^a$ for $a<e$; $\bar z$ for $a=e$) at both $\pi$ and
  $\bar\pi$, and compared with $v(c_{ij})$. I also check that $\mathrm{Norm}\,g_{ij}$ is supported on primes of $N$.
* (3): $\Phi$ is solved from the intersection matrix of the chain in exact rationals, for five random zero-sum charges
  per prime and for $2D_0$ (compared with $4Q_p$).
* (4): checked as an exact big-integer identity.

**All 26 pass** (`kdict_out.txt`: `ALL OK: True`). **Theorem A survives.**

It is, however, a re-labelling of known quantities: profile, reduced chord and cut slack. It carries no new
inequality. The proposal says so itself ("the global equality itself is just the product formula for the
discriminant").

## 2. Re-derivation of proof steps (B), (C), (D)

### 2.1 Theorem B: survives as a reformulation

* **Normalisation, re-derived.**
  * With $|x|_\infty=|x|^2$, the arc gives $-\sum_{i<j}\log|z_i-z_j|_\infty\ge-\binom M2\log(C^2R)$.
  * The places over $p\mid N$ contribute $W(\frac{M^2}4-\frac M2)+(Q_p+2\,\mathrm{coll})\log p$.
  * With $W=2\log R$, the difference is $\binom M2\log R-(\frac{M^2}2-M)\log R=\frac M2\log R$. Halving gives exactly
    the right side $\frac M4\log R+\binom M2\log C$ of `two_thirds/proof.md` (5.1).
  * The weights $\log q$ at split and inert $q$, and $\frac K2\log2$ at $(1+i)$, match (5.1) term by term.
* **"Free levels".** The proposal does not say whether the minimum over free levels equals the constrained bracket
  $\Psi^{(e)}$, which has $n_0,n_e\ge1$. It does. If $n_0=0$, move one point from the lowest occupied level $t_0$ to
  level 0. The cut sum changes by $t_0(1-M)<0$, and the collision terms do not increase. So the unconstrained optimum
  already has $n_0,n_e\ge1$.
* **Correction (minor, constant level).** Theorem B states that the $(2/3+\varepsilon)$ inequality "**is**" the adelic
  Fekete inequality. That is true only after the genuine $(1+i)$ minimum is replaced by the crude
  $\binom M2\log2$. The local points of $x^2+y^2=N$ at $(1+i)$ form a torsor of the norm-one torus over
  $\mathbb Z_2$.
  * They fall into $1,2,4,4,8,8,16,16,\dots$ classes mod $(1+i)^a$, i.e. about $2^{a/2}$. The unit count is $2^{a-1}$.
    (`k2adic.py`: identical for $N\equiv5,13,17,29,2193\bmod 2^{12}$.)
  * The genuine 2-adic minimum $\sum_aE(M,\mathrm{cl}_a)$ is $85$, $2809$, $48336$ at $M=10,50,200$, against
    $\binom M2=45,1225,19900$: a ratio of $1.89$, $2.29$, $2.43$.
  * So the true adelic Fekete inequality is stronger than (5.1) by $\Theta(M^2)\log2$, the same order as Remark B.1's
    archimedean gain. Constants do not matter here, but the identification "is" is inexact.

### 2.2 Lemmas 2.1–2.3: survive

* **Lemma 2.1.** $\bar K(s)^*$ is generated by constants and linear forms.
* **Lemma 2.2.** The Berkovich convex hull is determined by pair valuations.
* **Lemma 2.3.** In genus 0 the Néron–Tate height vanishes, and a symmetric form on $\{\sum c=0\}$ is determined by
  polarisation from its values at $e_i-e_j$.

All three are standard and correct. None is new.

### 2.3 Prop B.2 (ADC capped): the conclusion survives; the proof has two gaps

The ADC class is defined generously. Clause (ii) admits **any** local lower bound at **each** place, and the proof
then runs through cases (a)–(e).

* **Gap 1: the place $(1+i)$ is omitted (minor-to-major for the proof, minor for the conclusion).** ADC (ii) at
  $(1+i)$ contains the genuine 2-adic pigeonhole of §2.1, about $2.4\binom M2\log2$ at $M=200$.
  * The cited witness is `round5/construct.md` at $(27+o(1))M\log M$. Its residue data are defined only at **odd**
    prime powers (`construct.md` Def. 1.2: "$Q$ = odd prime powers"). Its pair gap (R) carries a single factor
    $\sqrt2$, i.e. one $(1+i)$ contact per pair.
  * The round-5 referee already flagged this: "$q=2$ should be included" (`round5/README.md`, construct.md
    corrections).
  * Cases (a)–(e) never mention $(1+i)$. So the proof does not show that the 27-witness satisfies the 2-adic ADC
    inputs.
  * *Repair:* cite `round5/lower.md`, whose witnesses satisfy 2-adic and twisted collisions at $X=2M$ according to its
    referee. That gives $(48+o(1))M\log M$, still at the rate $\log R/\log\log R$.
  * The cap therefore stands, with a different witness and constant.
* **Gap 2: the range $M\le X<2M$ (minor).** Case (c) treats only $q>X$, and case (b) treats only the designed depths.
  At $q\le X$ with $q^{a-1}\le X<q^a$, the first undesigned depth has $(q-1)q^{a-1}$ classes, which can be $<M$.
  * Example: $q=5$, $M=X=24$. Depth 2 has $20$ classes and forces $E(24,20)=4$ collisions mod $25$, which the design
    does not supply.
  * For $X\ge2M$ the count is $\ge q^a(1-1/q)>X/2\ge M$ (split) or $>q^a>X$ (inert), so the gap closes.
  * The statement should read $X\ge2M$, which is where the witnesses live in any case.
* **Scope (the author's own caveat, made precise).** The genuine-places version of ADC uses each pair product
  formula only as the inequality "finite contacts $\le$ archimedean budget". That is exactly the pair case of axiom
  (R). So **ADC is row 2 of overview §5.4, re-expressed in Arakelov language, and nothing more.**

### 2.4 Prop 3.3(b): survives, and is a restatement of row 8

Re-derivation of the argument:
1. Write $L=aH-bE$ and $h_E\ge\frac12\log R'-\log C-O(1)$. Then $h_L\le(a-\frac b2)\log R'+O(1)$.
2. A contradiction off a fixed closed set needs $h_L\ge\varepsilon h_H-O(1)$, i.e. bigness of $(a-\varepsilon)H-bE$,
   i.e. $b/(a-\varepsilon)<\tau^*$.
3. Combining, $b/\tau^*<a-\varepsilon<b/2$, which forces $\tau^*>2$.

The remaining case also fails. For $b/a\in[\tau^*,2]$, $L$ is not pseudo-effective, so the essential minimum is
$-\infty$ and Zhang/Northcott give nothing.

The argument is correct. It is the τ\* barrier of `threshold/upper.md` (row 8; checklist item 3) restated for
adelically metrised bundles. The proposal explicitly does not treat classes involving exceptional divisors over
$\Delta_{ij}$ (the determinant method, item 6). **Capped, not new.**

### 2.5 Theorem C: survives as stated, and is near-tautological

* **(1)** Keel: $\mathrm{Pic}$ is generated by boundary divisors. Correct.
* **(2)** A relation $\sum m_TD_T=\mathrm{div}f$ makes $f$ a unit on $M_{0,M}$. The units of a hyperplane-arrangement
  complement are monomials in the linear forms, i.e. products of $z_i-z_j$ in the Kapranov chart. Correct.
* **(3)** "Every inequality that follows from local inequalities and product formulas holds for every object that
  satisfies local inequalities and product formulas." This is a tautology once "pseudo-adelic configuration" is
  defined.
  * The cap content sits entirely in the claim that A_X fakes have completions. That claim rests on
    Prop B.2(e) (axiom (R), pair case), so it inherits the 2-adic gap of §2.3. Repair: use the `lower.md` witnesses.
  * Completions use **formal places of real weight**, which are not objects of the record's fake class (§5.1 of the
    overview: angle models). With genuine places only, Theorem C(3) degrades, by the author's own *Scope*, to the
    per-pair inequality. That is again row 2.
  * So Theorem C proves no cap beyond row 2. What it leaves uncapped (exact integrality of all
    $\mathrm{Norm}(z_i-z_j)$ at once) is, used in full, the problem itself (`round10/analytic.md` Thm D; overview
    §5.5, "additive data ... used in full they are the problem itself").

## 3. Multi-point term (E) and Theorem C

**§4.2 table recomputed** (`kadelic.py`, own code). The local minima are:
* pigeonhole at $q\nmid N$, split or inert;
* exhaustive minimisation over level vectors at $p\mid N$;
* $\binom M2\log2$ at $(1+i)$;
* Fekete on the chord segment at $\infty$, with Lobatto $\delta_M$.

All 15 arc rows agree with the proposal to three decimals: slack, $\infty$, $q\le M$, $p\mid N$, $q>M$, and the Fekete
gain ($\delta_{10}=0.70171$, $\delta_{50}=0.54949$, $\delta_{200}=0.51530$). Evidence only.

One caveat on the labelling. The $\infty$ column uses the segment of the **chord**, but the points lie on the arc,
whose capacity $R\sin(\theta/4)$ exceeds the chord-segment capacity $R\sin(\theta/4)\cos(\theta/4)$. So $\epsilon_\infty$
there is not a valid lower bound for arc points. The effect is $O(\theta^2)$, which is negligible here.

**§4.3 recounted** (`kmulti.py`, own code). For split $q\nmid N$ and points of equal norm,
$z_i\equiv z_j\bmod\varpi$ if and only if $q\mid z_i-z_j$ in $\mathbb Z[i]$, so I count classes mod $q$ directly.
* All per-row counts agree. For example, the 10-cluster has 4 pairs (at $83,89,109,197$) and 0 triples, and the single
  triple is at $q=7$ in `w6_2576450045`.
* **The stated totals are wrong.** The pairs total **29**, not 28 as stated in §4.3/§0. The expectation is
  **25.67**, not 25.5. The triple expectation is $0.023$, which is consistent with "0.02".
* No contact at depth $\ge2$ ($q^2\mid z_i-z_j$) occurs at any $q>M$.

**The cut-off should be $X$, not $M$ (Prop 6.1(b)).** An $A_X$ fake with $X\ge M$ carries designed residues at every
$q\le X$, so $\Gamma^{(M,X]}$ is fake data, not only $\Gamma^{\le M}$. With the witness value $X=2M$, the only
observed large-prime triple ($q=7\le12$, $M=6$) falls inside the designed range. So **$\Gamma^{>2M}_T=0$ for every
$|T|\ge3$ on all 15 arc clusters, with no exception**. The non-separation claim is therefore stronger than stated.

**It is also a theorem, not evidence.** By Theorem C(3), any inequality on $\Gamma$ that comes from local data plus
product formulas holds on pseudo-adelic completions. A completion may in addition place any ultrametrically
consistent multi-point pattern at private places, provided the pair budgets cover it. So $\Gamma^{>X}$ can never be
the separating input of a local-plus-product-formula argument. The proposal's §6.1 conclusion ("cannot be the
separating input") is therefore correct, but for structural reasons, and the 15-cluster evidence adds nothing.

**Kill test (i): fakes satisfying the inputs at barrier scale.**
* All of ADC: the `round5/lower.md` witnesses, $(48+o(1))M\log M$ with 2-adic collisions at $X=2M$. The
  `construct.md` witnesses at $27M\log M$ cover the odd places only (§2.3).
* The Zhang/bigness class: no fake is needed, because it is a variety cap (row 8).
* The Theorem-C class on $\overline M_{0,M}$: the same witnesses, completed at formal places, or at genuine places
  with the pair product formula as an inequality.

The only input that fakes fail is global rationality/integrality. Imposed in full this is the lattice-gap condition
at every prime, which is the problem itself (`round10/analytic.md` Thm D).

## 4. Needed statement vs target / Vojta

The proposal states the needed statement in two forms and calls them "equivalent". **They are not equivalent, and
neither is a usable target.**

* **Form 1** reads: "an inequality that uses global coherence of the local configurations across all places, weaker
  than lattice gaps at every prime, not via a fixed-$S$ exceptional set."
  * This is a description of a *class* of possible inputs, not a statement. No inequality is written down, so
    nothing can be tested or capped.
  * Its two boundary conditions are already the record's: round 10 Thm D (full lattice gaps are the problem), and
    lit18 Props A, B together with the project constraint (no fixed-$S$ exceptional sets).
* **Form 2** reads: "a proved $D=0$ Vojta-type height gap on $\mathrm{Bl}_YX_7$ or $\overline M_{0,8}$, i.e.
  $(L_8)$." This form is **at least as strong as the target**:
  * $(L_8)$ alone implies the uniform theorem (overview §4.2: "conversely $(L_8)$ alone gives the uniform theorem").
  * Vojta $D=0$ on $\mathrm{Bl}_YX_7$ gives it as well (overview Thm 4.1), and so does $V_8$ (Thm 4.3).
  * The uniform theorem implies every growth improvement. So Form 2 is a hypothesis stronger than what is to be
    proved.
  * Its status after rounds 18–19 is open, untestable (census Thm E: minors $\ge10^{49}$), not provable by auxiliary
    polynomials on $\overline M_{0,k}$ (item 3), and not reachable through unit pencils under abc (item 8).
* **"Equivalently" is an overclaim.** Form 1 could in principle contain growth-only statements strictly weaker than
  $(L_8)$. The lens produced none. Theorem C shows that no such statement comes from heights on $\overline M_{0,M}$
  built from local data and product formulas.

So the residual is the record's conditional lane, unchanged (kill criterion (iii)).

## 5. Barrier checklist, item by item

For each item: where the lens lands, and my reduction.

1. **Fake-satisfied inputs.** ADC is row 2, re-expressed (§2.3). Theorem A turns pair/cut/profile quantities into
   intersection numbers; Theorem B is the pigeonhole/chain bracket. Néron–Tate data on covers are pair data plus
   residue/character data. **Capped**, via the `lower.md` witnesses at $(48+o(1))M\log M$ once the 2-adic place is
   included.
2. **Averaged statements.** Quantitative equidistribution is an energy, i.e. a pair statement. Galois-orbit
   equidistribution does not apply (bounded degree, large height). **Not used / pair.**
3. **Auxiliary polynomials and bigness.** Prop 3.3(b) is the τ\* barrier (row 8) in adelic language. Slopes on $X_M$
   are Liouville. **Capped.**
4. **Single linear form; SLW.** Not used.
5. **$|S|$-only thresholds.** Mordell–Weil packing gives counts $e^{O(\rho)}$, with $\rho$ up to $M-1$ (fakes have
   full rank; overview Thm 3.12: bounded rank is *equivalent* to the uniform theorem). **Lim**, same shape as item 5.
6. **$(R_H)$ / CIRCLE$(q)$.** Classes involving diagonal exceptional divisors, and finite-place vanishing filtrations,
   are deferred to the determinant method. **Not reopened.**
7. **$(L_5)/(L_6)$ mechanisms at $k=8$; testing $(L_8)$.** The residual *is* $(L_8)$/Vojta. No test is possible
   (census Thm E). **Stopped.**
8. **The abc unit-pencil lane.** Lemma 2.1 shows that all multiplicative relations on $\mathbb P^1$ are pair
   relations. The three-term abc relation is squarefree-vacuous (round 8). **Closed.**
9. **Dirichlet exponent** (`round4/hardness.md`). Per pair the lens gives Liouville on the torus, which is attained.
   Jointly it gives exponent exactly $1/2$ for the diagonal line. **Untouched.**
10. **Project constraints.**
    * No fixed-$S$ exceptional set is used.
    * Gaussian gcds are kept exactly ($\gamma_{ij}=\log\mathrm{Norm}\,c_{ij}$, verified).
    * No phase is set to zero.
    * $N$-dependent metrics are explicitly excluded.

    **Respected.**

## 6. What survives

**Nothing survives as a lead toward growth or uniformity.** Killed: the residual needed statement is at least as
strong as the target (§4), and every inequality the lens can produce is capped (§5).

What survives as correct mathematics. All of it is reformulation or evidence and gives no growth input.
* **Theorem A** (exact; independently re-verified on 26 configurations). The exponent profile is the set of components
  of $\mathcal X_N$ met by the sections, $(P_i\cdot P_j)_{\rm fin}=\log\mathrm{Norm}\,c_{ij}$, the vertical correction
  of $D_0$ is $Q_p$, and the cut identity holds.
* **Theorem B** as a reformulation: the $(2/3+\varepsilon)$ finite inequality is the adelic Fekete inequality with
  the crude $(1+i)$ bound. The genuine 2-adic local minimum is larger by a factor of about $2$–$2.4$ in contacts.
  This is a constant only.
* **Prop 3.3(b)**: correct, a restatement of row 8.
* **Prop B.2 and Theorem C**: correct in conclusion. B.2 needs the `lower.md` witness for the $(1+i)$ place and
  $X\ge2M$. C proves no cap beyond row 2.
* **Evidence.** The §4.2 slack table is reproduced exactly. The §4.3 per-row counts are reproduced, with corrected
  totals of 29 pairs (expected 25.67) and 1 triple (expected 0.023). Moreover $\Gamma^{>2M}_T=0$ for all $|T|\ge3$ on
  all 15 arc clusters.

**Corrections the write-up needs:**
1. §0 and §4.3: "pairs 28 ... vs 25.5" should read "29 ... vs 25.7".
2. Prop B.2: add the $(1+i)$ place. The 27-witness has odd-prime residues only, so cite `lower.md`
   ($48M\log M$), and state $X\ge2M$.
3. Theorem B: "is" should become "is, with the crude $(1+i)$ bound".
4. Prop 6.1(b) and §6.1: cut at $X$, not $M$. The non-separation follows from Theorem C, not from the evidence.
5. Needed statement: drop "Equivalently". Form 2 implies the uniform theorem; Form 1 is not a statement.
