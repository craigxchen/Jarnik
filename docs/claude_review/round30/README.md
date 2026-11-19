> **Added after round 31 (Part IV, `../round31/README_IV.md`).** The *uniform* bound in this model is now a theorem:
> $M(D,0)\le24$ for every $D$ and $\sup_DM(D,c)<\infty$ for every $c$, over $\mathbb C(t)$, $\mathbb R(t)$,
> $\mathbb Q(i)(t)$ (Langer's orbifold BMY on the section surface). Statements below that the uniform question is
> open, or that it has the shape of open truncated function-field Vojta, apply only to the sharp statement (R$_4$)
> ($\limsup M=8$), which stays open. The upper bound $4(D+2)$ in the small-cases table becomes $\min(4(D+2),24)$.
> The conditional "unbounded $\mathbb Q(i)(t)$-clusters would refute $A_5(\delta)$" has a false premise.

# Round 30 -- judge report: the function-field model

Lanes: `search.md`/`kill_search.md`, `prove.md`/`kill_prove.md`, `bridge.md`/`kill_bridge.md`. Judge checks are in
`judge_checks/` with notes in `judge_notes.md`. **Every claim labelled [J] is the judge's own and is unrefereed.**

Model and notation as in the round brief: $c=0$ unless stated, $K=\lfloor D/2\rfloor+c$, $L=D-1-K$. "Proved"
means a complete proof, or an exact or exhaustive computation that has been reproduced independently. "Certified"
means a Krawczyk or exact certificate. "Evidence" means numerics or a non-exhaustive search.

## 0. Verdict in five lines

1. **No rank-5 cluster at $c=0$ was found for $D\le12$.** This is search evidence, not proof. The first clause of the
   stop rule holds.
2. **There is no general rank bound.** Rank bounds were proved only on two very special loci: the cube surface
   (largest cluster 8, rank 4) and the horizontal line ($M\le2$). Read literally, the second clause of the stop rule
   is therefore not met. In substance, both bounds are of the kind the rule meant to exclude (§4).
3. **The best tool of the round is the saturated-lattice balance lemma** (search Lemma 1.1 = prove Lemma 1.2). It
   gives Prop K1 and the $D=8$ equidistant non-existence. An exhaustive search with it (§3.2) gives the judge's
   results **$M(7,0)\le5$, $M(8,0)=8$ and $M(9,0)\le5$** [J, unrefereed]. It does **not** kill rank 5: 21 classes of
   6-point rank-5 codes survive at $D=8$. So it settles sizes in small cases, not ranks.
4. **(S4) is reconciled.** Vojta's truncated inequality for the SNC pair $(\mathrm{Bl}_eX,\pi^*\partial+E)$ gives
   $r(L+1)\le|S|-1+\varepsilon h$, so (R$_4$) has exactly the Vojta shape (bridge is right; coordinator and prove used
   a weaker consequence).
5. **Nothing changes for the integer problem.** The best cube template constant is $C_\infty=145877.66$, far above
   the golden $5.8977$. No record improves, and Question B gains no A-points.

**Decision: close the function-field line for internal work.** One bounded exception is allowed: the finite $D=8$
check of §7 (21 explicit codes).

## 1. Killer-corrected statements per lane

### 1.1 search (killer: survives with corrections; one false claim)

**Proved**
* **Lemma 1.1 (saturated balance).** Let $V$ be the real solution space of the moment equations, and let
  $0\ne a\in V\cap\mathbb Z^D$. Then $\|a\|_1\ge L+1$ and $\min(a^+,a^-)\ge L/2$. This needs $N$ squarefree. The
  prove lane proved the same lemma independently.
* **Prop 2.1.** Half-cube codes need $4\mid D$.
* **Thm 3.1.** On the nondegenerate cube locus, the cluster through $Q$ has exactly 8 points and rank 4. One
  degeneracy line in the proof is misstated, but the conclusion stands. The prove lane's Thm 3.1 strengthens this to
  every base point.
* **Horizontal line.** $M\le2$ for $D\ge7$, by Plotkin with even $d=2L$.
* **Exact balanced clique numbers:** $M(8,0)\le14$ and $M(9,0)\le6$.
* **Rank-$\le4$ code maxima for $6\le D\le16$.** These come from an exhaustive enumeration solved with a floating
  MILP; the killer filled in $D=11,13,15$. An 8-point rank-$\le4$ code exists only for $D\in\{8,12,16\}$.
* **Prop K1 (killer).** At $D=10$ and $D=14$, every cluster with $\ge6$ points has rank $\ge5$. The 141 "codes"
  form one class, which Lemma 1.1 kills at $D=10$; the two classes at $D=14$ die the same way.

**Certified**
* $M(9,0)\ge4$ and $M(10,0)\ge5$.
* A non-composed 8-point rank-4 half-cube$(2,1)$ cluster at $D=12$, over $\mathbb R$. It is not shown to be a
  $\mathbb Q(i)$ template.
* Rank-5 six-point clusters at $c=1$ for $D=10$ and $D=12$, beyond the covering range.
* (killer) Rank-4 five-point clusters on the vertical-line locus at $D=10$.

**Withdrawn**
* "Vertical lines give rank $\le3$". It was a parameter count, and it is refuted. [J] I rechecked killer solution 2
  numerically: it has 5 zero-sum sets, including $\emptyset$, of rank 4.
* "141 codes never realised".
* "16/46 types never realised". At most 11 remain unrealised under the killer's solver, and two of the lane's
  unrealised types were certified.

**Evidence:** $M(7,0)=4$, $M(9,0)=4$, and no rank 5 at $D=8$.

### 1.2 prove (killer: survives; "balance does not separate" withdrawn)

**Proved**
* (S1)–(S3), with (S3) corrected: the $N_v$ version is false and the Blaschke version is right.
* Lemma 1.2 (as search Lemma 1.1).
* Props 2.1–2.3: minimal pairs are real-polynomial fibres, the symplectic three-point identity, and the ordered
  phases of equidistant clusters.
* Prop 3.3.
* **Thm 3.1 (cube surface).** Every nondegenerate real member has largest cluster exactly 8 over all 256 base
  points, and that cluster has rank 4. The 48 special real members (one $|y|$-set up to symmetry) have profile
  $(8,4,\dots)$. The killer recomputed this with resultants and rigorous root isolation. The scope is the symmetric
  surface only.
* **Thm 3.2.** $M\le2$ on horizontal lines for $D>6c+6$.
* **(killer)** There is no 6-point equidistant cluster at $D=8$: each of the 3 code classes contains a weight-2
  integer vector in the span of its differences.

**Certified (killer)**
* Nondegenerate 5-point rank-4 equidistant clusters at $D=8$. The fibre profile is $(6,5,5,4,\dots)$, not the
  lane's $(5,5,4,4,\dots)$.

**Withdrawn**
* Prop 2.5 as used: "balance does not separate real clusters from fakes". Only *pairwise* balance is faked. The
  Theorem-B fakes and their doubles contain integer kernel vectors of weight $D/4<L+1$.

**Heuristic only:** the dimension bound $\dim Y\le3$ for deformations, and the failure analyses of the Wronskian,
Mason/gcd and Lehmann–Riedl–Tanimoto (LRT) tools. The failure analyses are sound as explanations.

### 1.3 bridge (killer: survives, minor corrections)

**Proved**
* **Template theorem 1.2.** This is a routine formalisation of (S6). A primitive $\mathbb Q(i)(t)$-cluster with $D$
  even and $c=0$ specialises at rational $t_0\to\infty$ to integer clusters, with $C(t_0)\to\Gamma/\sqrt{|g(t_0)|}$.
  The best constant is $C_\infty=\Gamma/\sqrt{\prod_\ell g^*_\ell}$, attained along an arithmetic progression. When
  $K<D/2$, $C\to0$.
* **Lemma 3.2.** The conjugation-symmetric squarefree 8-point cluster carries the cube code. The 5.754 record cluster
  has exactly the cube code (exact factorisation).
* **The 56 five-subsets of a cube cluster.** All have $\alpha=4$, so none is an A-point. None lies on a quadric cone
  or a cubic threefold (exact). They sweep 3-folds; the killer proved this with an exact Jacobian. The 24 subsets that
  contain two conjugate pairs are dense in subtori $u_au_b=u_cu_d$.

**Exact**
* $C_\infty=145877.658\ldots$ at $y=(-\frac94,-\frac23,\frac14,\frac83)$, with $t_0\in\frac12\mathbb Z$.
* A 74-digit integer example, checked exactly.
* Sub-cluster constants: $30.88$ (5 points), $432.87$ (6), $7946.5$ (7).

**Evidence**
* That this is the minimum over the cube family (28 of 38 classes computed).
* Membership in $W$ and $\iota W$.

**Corrections**
* At $D=7$ only Plotkin's bound 8 was quoted as proved. In fact the balanced clique bound 7 (search) is proved, and
  [J] $M(7,0)\le5$ (§3.1).
* $3\sqrt[3]2$ is the characteristic-0 Corvaja–Zannier constant.
* The $\varepsilon<1/4$ threshold depends on the normalisation of $h$.

**Conditional, and correct:** truncated function-field Vojta for $(\mathrm{Bl}_eX_r,\pi^*\partial+E)$ gives
$r\mu\le|S|-1+\varepsilon h$.

## 2. Small cases ($c=0$, squarefree $N$)

| $D$ | $L$ | lower bound | upper bound | rank facts | status |
|---|---|---|---|---|---|
| 2–6 | | 4, 3, 6, 3, 5 | the same | | exact (round 29); [J] filter reproduces $D=5,6$ |
| 7 | 3 | 4 (exact) | **5 [J]** (was 7: balanced clique) | the only 5-point type is $\mathcal T_7$; rank $\le4$ | $=4$ is evidence (3000 starts) |
| 8 | 3 | 8 (cube, exact) | **8 [J]**, so **$M(8,0)=8$ [J]** (was 14) | no rank-5 found; equidistant 6-point impossible (proved) | rank-5 open |
| 9 | 4 | 4 (certified) | **5 [J]** (was 6) | rank $\le4$ [J]: no 6-point code survives | $=4$ is evidence |
| 10 | 4 | 5 (certified) | 12 (Plotkin $A(11,6)$) | **K1:** $\ge6$ points needs rank $\ge5$ | $M=5$ unless rank 5 |
| 12 | 5 | 8 (certified, half-cube$(2,1)$) | Plotkin | no rank-5 extension among 918 sets (evidence) | |
| 14 | 6 | — | — | **K1** | |
| $8\mid D$ | | 8 (cube family, exact over $\mathbb Q(i)(t)$) | $4(D+2)$ | rank-$\le4$ max 8 only at $4\mid D$ ($D\le16$) | |

At $c=1$:
* $M(7,1)\in[8,14]$: the lower bound is exact, the upper bound is the proved balanced clique.
* Rank-5 six-point clusters are certified at $D=10$ and $D=12$; they lie outside the covering range $D\le8$.
* Random rank-5 types fail at $D=14,16$. This is thin evidence.

## 3. Reconciliation and judge checks [J]

### 3.1 Contradictions between lanes

**(S4): "(R$_4$) is sharper than Vojta" (coordinator, prove) against "exactly the Vojta shape" (bridge).**
Resolved in favour of bridge. Let $X$ be a smooth toric compactification of $\mathbb G_m^r$ with boundary $\partial$,
so that $K_X+\partial=0$. Let $\pi:X'=\mathrm{Bl}_eX\to X$ be the blow-up at the identity $e$. Then
$K_{X'}=\pi^*K_X+(r-1)E$.
* Since $e\notin\partial$, the divisor $\pi^*\partial+E$ has simple normal crossings, and
  $K_{X'}+\pi^*\partial+E=rE$.
* A rank-$r$ cluster gives a curve $\mathbb P^1\to X'$. It meets $\pi^*\partial$ only above the $\le2D$ roots of $N$,
  and meets $E$ with total contact $h_E\ge m=L+1$ at $t=\infty$, plus $1$ for each further point.
* Truncated Vojta for this pair gives $r\,h_E\le N^{(1)}_{\pi^*\partial}+N^{(1)}_E-2+\varepsilon h$, that is
  $$rm+(r-1)\#\{\text{other contacts}\}\le|S|-1+\varepsilon h .$$
* The version without $E$ in the divisor, $(r-1)m\le|S|-2+\varepsilon h$, follows from this one by subtracting
  $m\ge1$. It is a *weaker consequence of the same conjecture*, not a competing statement.
* At $c=0$ ($m\approx D/2$) the strong form forbids $r=5$ when $\varepsilon h<D/2+1$ and allows $r=4$. At
  $\varepsilon=0$ it is round-29 Theorem D, $r(L+1)+1\le2D$. The cube family (rank 4, $4m=2D>2D-1$) is consistent
  with this: it must use the $\varepsilon h$ slack, so it lies in the exceptional set for small $\varepsilon$, and it
  is not dominating.
* **Conclusion:** (R$_4$) is precisely the prediction of function-field Vojta outside the exceptional set. It is not
  stronger. Function-field Vojta for blow-ups of tori is open (memory-quote; no source was found this round), so this
  is calibration, not evidence.

**"Balance is faked" (prove Prop 2.5) against kill_prove.** Resolved for kill_prove. I re-derived the witness:
$\mathbf 1_H=\frac14(\mathbf 1+\chi_{110}+\chi_{101}+\chi_{011})$, $H=\{000,111\}$, has weight $2<4$. Being
non-negative, it also fails Gordan. Pairwise balance is faked; saturated balance is not.

**The $D=7$ upper bound.** The values on record are 8 (Plotkin, kill_bridge) and 7 (search, balanced clique, proved).
The judge check below gives 5.

**The vertical locus.** The search lane claims rank $\le3$ there. kill_search has a certified rank-4 counterexample,
and the claim is withdrawn.

### 3.2 Judge computation: exhaustive saturated-balance codes [J, unrefereed]

**Method.** `judge_checks/satdfs.py` and `sym8.py` enumerate every code in $\{\pm1\}^D$ ($c=0$) that satisfies all of
the following:
* pairwise distance $\ge L+1$;
* pairwise balance $\min(p,n)\ge L/2$;
* at every node, the saturated condition of Lemma 1.1.

The saturated condition is tested exactly: an integer nullspace in python-flint, then an int64 test of every integer
vector $a\in[-2,2]^D$ that violates $\|a\|_1\ge L+1$ or $\min(a^\pm)\ge L/2$. This is a subset of Lemma 1.1, so every
prune is a proof. Each run is exhaustive up to $S_D\times\{\pm1\}$: every code contains a vertex that can be moved to
the representative $(+^w,-^{D-w})$, $w\le D/2$.

| $D$ | balanced clique (pairs only) | **saturated-filter maximum** | known lower bound |
|---|---|---|---|
| 5 | 4 | **3** | 3 (exact) |
| 6 | 7 | **5** | 5 (exact) |
| 7 | 7 | **5** | 4 (exact) |
| 8 | 14 | **8** | 8 (cube, exact) |
| 9 | 6 | **5** | 4 (certified) |

The filter reproduces the true values at $D=5,6$. **Consequences [J]:**
* **$M(7,0)\le5$** (was 7).
* **$M(8,0)=8$ exactly** (was $8\le M\le14$). The upper bound is the filter maximum; the lower bound is the exact cube
  cluster.
* **$M(9,0)\le5$** (was 6).

These bounds hold for squarefree $N$ and rest only on Lemma 1.1, which both lanes and both killers proved. The search
was run once, with one implementation (runtimes: 8 min at $D=8$, 13.5 min at $D=9$). Before these are quoted as
theorems, a referee must rerun them with independent code.

At $D=8$, **21 classes** of 6-point rank-5 codes survive the filter (`judge_checks/classes_D8_c0.pkl`). Every
surviving code has at most 8 points, so rank 5 at $D=8$ cannot raise $M$, but the filter does not exclude rank 5.
At $D=9$ no 6-point code survives at all, so every cluster there has rank $\le4$ [J].
At $D=7$, $M\in\{4,5\}$; the survivors at size 5 are one $S_7$-orbit, which matches search's single type
$\mathcal T_7$.

**A rank-5 code that passes every linear filter (D=8) [J].** The rows are
$$\begin{smallmatrix}+&+&+&+&+&-&-&-\\+&+&+&-&-&+&+&-\\+&+&-&+&-&+&-&+\\+&-&+&+&-&-&+&+\\+&-&-&-&+&+&+&+\\-&+&+&-&+&+&-&+\end{smallmatrix}$$
* 13 pairs have $(p,n)=(2,2)$ and two pairs have $(3,3)$.
* The rank is 5, and the span is saturated (index 1). Its LLL basis consists of weight-4 vectors with
  $a^+=a^-=2$.
* No integer combination with coefficients in $[-4,4]$ violates Lemma 1.1, and Gordan's condition is feasible.
* Realisation attempt: 400 Levenberg–Marquardt starts in the 14-dimensional gauge-fixed space (15 equations). Every
  converged solution is degenerate ($\rho_k=\rho_l$); the best nondegenerate residual is $1.8\cdot10^{-4}$. This is
  evidence of non-realisability, which is expected, since the expected dimension is $16-2-15=-1$.
* So the saturated-lattice route proposed by kill_prove **cannot by itself prove (R$_4$) at $D=8$**. The step left
  over is non-linear.

## 4. Verdict on (R$_4$) and on the pre-registered stop rule

**(R$_4$) is open.** This round neither proves nor refutes it.

For:
* No rank-5 cluster has been found at $c=0$ for $D\le12$. This is evidence. The searches covered all 46 surviving
  five-point types at $D=8$ (plus one extra set), the $D=12$ half-cube family, the cube locus, and the $D=8$
  equidistant codes, where its absence is now proved.
* (R$_4$) has exactly the Vojta shape (§3.1).

Against:
* Rank 5 *does* occur at $c=1$ beyond the covering range ($D=10,12$, certified), so the dimension count is not a
  hard wall.

What the round does settle:
* Proved reductions: at $D=10$ and $D=14$ (Prop K1), $M=5$ resp. $M\le5$ unless rank $\ge5$ occurs.
* The linear tools (saturated balance plus Gordan) do not exclude rank 5 at $D=8$: 21 classes survive. They do bound the size: $M(8,0)=8$ [J, §3.2].

**Stop rule.** The rule fires when (a) no rank-5 cluster is found at $c=0$ for $D\le12$ and (b) no rank bound is
proved on any positive-dimensional locus.
* (a) holds.
* (b) fails on the letter. The cube surface and the horizontal line carry proved rank bounds.
* In substance both bounds are vacuous for the question:
  * The cube surface is the one known rank-4 family. "Rank exactly 4 on it" means it does not extend; it says nothing
    off the surface.
  * On the horizontal line the excess structure degenerates to Prouhet–Tarry–Escott (PTE) moments, and $M\le2$ there.
  * Neither proof technique (an excess combinatorial $m=2$ equation; Chebyshev rows) generalises to an open set of
    root configurations.
* The two arguments offered for staying open both collapsed under refereeing:
  * search's "real-algebraic obstruction" lead turned out to be linear algebra (Prop K1);
  * prove's "balance is faked" turned out wrong, but the corrected filter does not suffice either [J].

**Judge ruling: the line is closed for internal work.** The rule is treated as having fired in spirit, and the
literal non-firing is recorded here explicitly. The one bounded exception is §7.

## 5. What is genuinely new and correct

1. **The saturated-lattice balance lemma** (search Lemma 1.1 = prove Lemma 1.2). It is the round's only tool that
   kills configurations, and it is elementary (the argument principle). Its consequences:
   * Prop K1;
   * no equidistant 6-point cluster at $D=8$;
   * the Theorem-B fakes are not genuine;
   * [J] $M(7,0)\le5$, $M(8,0)=8$, $M(9,0)\le5$, and rank $\le4$ at $D=9$ (exhaustive, unrefereed).
2. **Exact or certified objects:**
   * the non-composed $D=12$ 8-point half-cube$(2,1)$ cluster;
   * $c=1$ rank-5 clusters at $D=10$ and $D=12$;
   * certified 5-point equidistant rank-4 clusters at $D=8$;
   * vertical-locus rank-4 clusters at $D=10$;
   * $M(9,0)\ge4$ and $M(10,0)\ge5$.
3. **Cube surface theorems.** Search Thm 3.1 and prove Thm 3.1: the cube family is isolated in rank. The 8-point
   structure does not extend, and the 48 special members change nothing.
4. **Half-cube arithmetic (Prop 2.1)** and the complete rank-$\le4$ code tables for $D\le16$. Conditional on (R$_4$):
   $\limsup M(D,0)=8$ can only come from $4\mid D$, and $M\le6$ at $D=10,14$; at those two $D$ the bound improves to
   $\le5$ by K1.
5. **The template theorem with exact local factors**, and the value $C_\infty^{\rm cube}=145877.66$.
6. **The (S4) correction** to the SNC pair $\pi^*\partial+E$.

None of these is a new barrier type. Novelty overall: a known barrier in a new model, plus exact data.

## 6. What changes for the integer problem

* **Templates.**
  * Function-field clusters over $\mathbb Q(i)(t)$ with $D$ even and $c=0$ give integer clusters with an explicit
    limiting constant (bridge Thm 1.2, proved).
  * The cube gives $M_\infty(C)\ge8$ only for $C>145877.66$. The golden family already gives $C>5.89771$, so **no
    lower bound for $M_\infty(C)$ improves**.
  * The $D=12$ cluster is certified over $\mathbb R$ only, so it is not yet a template.
  * Odd $D$ (that is, $K<D/2$) gives $C\to0$ templates. That route needs 6 or more points at $K=(D-1)/2$, and none is
    known.
* **Records.** The 5.754 record cluster carries the cube code, but it was not shown to be a cube value. The 5.674
  (rank 7) cluster and the 10-point cluster are not specialisations of any known $c=0$ function-field cluster.
* **Question B (exceptional set).**
  * Cube 5-subsets have $\alpha=4$ exactly, so they give no A-points.
  * They avoid all round-27 cones and cubic threefolds (exact).
  * They sweep 3-folds, and 24 of them are dense in subtori $u_au_b=u_cu_d$ (proved).
  * So the cube family lies inside the expected exceptional structure without producing a violation. **Nothing
    changes for $A_5$ or $(L_5^o)$.**
* **Conditional statements.**
  * Unbounded $\mathbb Q(i)(t)$-clusters would refute $A_5(\delta)$ for $\delta<1/30$ and $V_E(5)$ (via J4).
  * They would refute Cilleruelo–Granville only with bounded constants, and nothing in the model controls those.

## 7. Follow-up (at most one; optional, bounded)

**Decide the 21 surviving rank-5 classes at $D=8$.**

*Why it is worth it.* If none of the 21 is realisable, then every cluster at $D=8$, $c=0$ has rank $\le4$. That would
be the first full-$D$ rank theorem in the model, with $M(8,0)=8$ already in hand [J]. It would also show whether the
pattern "saturated filter plus elimination of a finite list" can be run at $D=12$, the first case where both the
half-cube and rank-5 candidates compete.

*Concrete first step.*
1. A referee independently reruns `judge_checks/sym8.py` at $D=7,8,9$.
2. For each class in `classes_D8_c0.pkl`, write the 15 moment equations $\sum_k\delta_k\,\mathrm{Im}\rho_k^j=0$
   ($j\le3$, 5 independent differences) in the 16 real unknowns.
3. Fix the real affine gauge ($\sum x_k=0$ and one $y_k=1$).
4. Prove that the resulting system has no solution with $y_k>0$ and pairwise distinct roots. Use Groebner bases or
   resultants over $\mathbb Q$ with the nondegeneracy saturation $\prod_{k<l}|\rho_k-\rho_l|^2\ne0$. The system is
   over-determined by one equation, as the expected dimension is $-1$.

*Evidence so far [J].*
* 60 Levenberg–Marquardt starts per class.
* For 10 of the 21 classes (classes 0–9), the attempt produced either no nondegenerate run or a best nondegenerate
  residual $\ge6\cdot10^{-5}$. The other 11 classes did not run before the session ended.
* The single code analysed in full (§3.2) gave 400 starts, all degenerating.

*Honest probabilities.*
* All 21 classes proved non-realisable within one round: about 0.4. Elimination with 14 unknowns and degree 3 may
  not finish without homotopy software, which is unavailable here.
* That this leads to a proof of (R$_4$) for all $D$: below 0.03.
* That it changes anything for the integer problem: about 0.01.

If the follow-up is not taken, nothing is lost: the line is closed.

## 8. Scores (0–10)

| lane | score | reason |
|---|---|---|
| search | 5/10 | The most new exact and certified data: the $D=12$ non-composed cluster, $c=1$ rank 5, and the full rank-$\le4$ tables. It shares the key lemma. One false "proved" claim (vertical lines), and its lead was mis-identified. |
| prove | 4/10 | Correct cube-surface and horizontal-line theorems with full proofs, and correct (S1)–(S3) corrections. Its central negative claim (Prop 2.5, "balance is faked") is false for the saturated form. |
| bridge | 4/10 | A correct template theorem with exact constants, and a correct (S4) correction (§3.1). It fixes the integer-side calibration: no record and no Question-B change. Little novelty. |
| **round, on the goal** | **2/10** | No rank bound in general. No change for Cilleruelo–Granville, $M_\infty(C)$ or Question B. The function-field line is closed. Small-case values sharpen ($M(7,0)\le5$, $M(8,0)=8$, $M(9,0)\le5$ [J]), but they do not touch the uniform question. |
