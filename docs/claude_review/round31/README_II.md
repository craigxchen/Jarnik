# Round 31, Part II -- judge report: Astra's second exchange

Sources (frozen):
* `astra_free_II/answer.md`, from the heading "# Part II" (§§II.1-II.5). These answer the coordinator's
  `astra_free_II/reply2.txt`.
* `astra_lemma_II/answer.md`, from "# Part II" (§§II.1-II.7). These answer `astra_lemma_II/reply2.txt`.

Eight referees each took one block:

| block | report | source sections |
|---|---|---|
| min | `ref2_min.md` | free II.1 |
| template | `ref2_template.md` | free II.2 |
| grid | `ref2_grid.md` | free II.3 |
| index | `ref2_index.md` | free II.4 |
| tri | `ref2_tri.md` | free II.5 |
| height | `ref2_height.md` | lemma II.1-II.3 |
| window | `ref2_window.md` | lemma II.4-II.5 |
| plane | `ref2_plane.md` | lemma II.6-II.7 |

Part I is judged in `README.md` in this directory, and this report relies on it where Part II uses Part I statements.
Judge scripts are in `judge_II/`. **Every claim labelled [J] is the judge's own and is unrefereed.**

Labels:
* "Proved": a referee checked the proof line by line, or reproduced an exact or exhaustive computation with
  independent code.
* "Computed": exact data, independently reproduced.
* "Evidence": numerics or a search that is not exhaustive.

No referee reported FALSE. No referee reported GAP on any claim, headline or otherwise. The judge therefore had no
disputed claim to rule on. He did read the passages behind the two coordinator questions singled out for §3: free II.4
lines 545-605, and lemma Part II lines 232-345 together with the II.6/II.7 statements.

## 0. Verdict in five lines

1. **Nothing changes for the integer problem.** The growth rate stays $(2/3+\varepsilon)\log R/\log\log R$. There is
   no uniform bound. Question A and route 5 are untouched, and Question B is open. Astra says the same of her own work.
2. **The most useful result is the index clarification (free II.4).**
   * In a shortened genuine one-unit cluster, all but $O(m^3)$ ordered five-tuples have index one. The constant is
     absolute, via Baker–Wüstholz.
   * A subcluster of $\ge m/145-O(1)$ points has *every* five-tuple of index one.
   * Hence **"index-one $(L_5^\circ)$ suffices for the uniform bound" is established.** It is formally weaker than
     $(L_5^\circ)$ and no easier in substance.
   * Part I's high-index cover theorem and its $\{1,3,5\}$ reduction are thereby irrelevant to clusters.
3. **The most interesting mathematics is the fixed-plane theorem (lemma II.6-II.7).** Ten-pair-cut tuples whose block
   norms lie on a fixed rational plane $U$ satisfy $T\gg_U(\max n_e)^{1/6}$, and so have exponent $\ge38/11$, or
   $40/11$ for comparable norms. The constants depend on $U$ and are ineffective.
   * [J] It is not the record's first gain over $10/3$ on genuine index-one tuples. Balanced Liouville already gives
     $15/4$, and Part I (7) gives $7/2$ with at most nine primes.
   * [J] It is the first gain over $10/3$ inside the ten-pair-cut profile class. That is where the gap-only angle
     model attains $10/3$, so it is also the first exponent bound in the record that a gap-admissible fake demonstrably
     fails on a nonempty family (§3.3).
   * It is thin, non-uniform in $U$, and below $15/4$.
4. **Everything else is correct and peripheral.**
   * The $D=8$ minimal-vector count is 12, attained only by the cube.
   * The finite specialization templates have clique number 3.
   * The free-nonvanishing grid designs are all subcritical (a closed negative).
   * The triangle quadratic-form identity is already in the record (round27 `stab.md` Prop 1.1(iii)).
   * The six-norm short-interval lemma is new.
   * There are no polynomial norm pencils, and covers cannot lift index-one points.
5. **Part II score on the goal: 2/10.** It is an error-free exchange that closes three questions the coordinator asked
   (2(a)/(b), 4, and lemma Q1-Q2), sharpens Question B to "index one", and moves no rate.

## 1. Referee-corrected statements per block

### 1.1 min: free II.1 (referee: all correct; one verification caveat)

**Proved** (function-field model, $D=8$, eight distinct roots in the upper half plane, $c=0$)
* **Reduction.** A weight-4 moment-kernel vector has masses $(2,2)$, by the round-30 mass lemma. The moment points
  $h_i=(y_i,x_iy_i,3x_i^2y_i-y_i^3)$ are distinct. Each 4-support carries at most one pairing. Vectors up to sign
  correspond to distinct supports, and so to distinct norm divisors.
* **No four moment points are collinear** (a monic quartic with no cubic term, plus Vieta). At fixed height, no three
  are collinear. This is a shorter proof of the Part I four-collinear lemma.
* **Main count.** *There are at most 12 weight-four ternary kernel vectors up to sign (24 oriented), with equality
  exactly for configurations carrying the eight-row cube code.*
  * Census by rank: 70, 1799, 20640, 94190, 70168 admissible flats.
  * The three planar affine classes (28/18/4 records) are excluded. One class reduces to $2b^2+1=0$.
  * The rational cube example has exactly 12 such vectors up to sign, and none of weight 1-3.
* **Repeated-factor variant.** With midpoints $2e_j-e_i-e_k$ allowed, the bound is still 12, with the same equality
  case. The hexagon rank-4 class is excluded by hand, and all 59 planar shapes are excluded.
  * *Verification caveat:* the referee's independent sympy saturation excluded 52 of the 59. The other 7 (out-indices
    18, 22, 26, 31, 37, 48, 51 in `ref2_min/classes_mid12.txt`) rest on Astra's certificates
    (`ii_min_all_verify.py`). The referee reran these and checked their positive-product logic.
  * This is a gap in independent reproduction, not in the proof. [J] No Singular, Macaulay2 or msolve is available
    here either, so the judge could not close it.
* **First gap.** A configuration with 11 directions has 12.
* **Minimum-14 example.** Roots $x+i$ with $x\in\{1,5,8,10,12,16,22,23\}$ give a kernel-lattice minimum of 14, with 13
  shortest directions up to sign.

**Corrections**
* Astra's "28 unordered pairs" should read: 24 pairs have weight 4 (12 supports, each twice) and 4 pairs have weight 8.
* Two slips in the coordinator's question are corrected:
  * the count is 12 divisors, not "28 edges up to sign";
  * weight $D/2$ is only a *lower* bound for the lattice minimum (the minimum-14 example), not "exactly the minimum".
* The "completeness" of the affine invariant is unproved. It is not needed, because the referee's $S_8$
  classification gives the same classes.

**Gap:** none.

### 1.2 template: free II.2 (referee: all correct)

**Proved / computed**
* **Configurations C1 and C2.** These have pair-sum coincidences, with 7 resp. 9 ternary weight-4 directions (9 resp.
  11 counting midpoints). Fibre histograms are 120/32/8/12 and 104/32/16/8. The referee checked these from the
  definition, without assuming the mass lemma.
* **Specialization cost.** With parent $N_0$, $d_v/(16t_v^2\sqrt{N_0})\to1/(16b_v^2)$. Taking $N=\mathrm{lcm}\,d_v$
  gains $\sqrt G$.
* **CRT progression.** Five (C1) resp. seven (C2) of the $d_v$ simultaneously satisfy every condition of the divisor
  test, on infinitely many parents. For each fixed $B$ this holds *eventually*.
* **Clique bounds.** C1: no four lattice points realize all five labels. C2: the seven ratios have no four-clique.
  Minimum-shell clique maxima are 4, 3, 3 (cube, C1, C2), and fibre maxima are 8, 4, 5.
* **Ellipse obstruction.** An axis-parallel ellipse with real centre gives at most two weight-$m$ directions.
* **Substitution.** $t\mapsto t^k+a$ gives no new $\mathbb Q(i)[t]$ constant-offset divisors.

**Corrections (scope)**
* "For each fixed $B>0$" holds only beyond $B$-dependent, astronomically large $n$. At the computed parents (773-1436
  digits) the least admissible $B$ is 1.0 to 2.4, and $d/\sqrt N$ stays near $\sqrt L\approx1.7\cdot10^{10}$.
* The referee's C2 contraction $L$ is 4329 times larger than Astra's, from extra gains at 3, 13 and 37.
* With all 7 or 9 directions the cutoff fails ($L\approx10^{13}$ to $10^{16}$). This is why only subsets are used.

**Gap:** none. These are fixed finite counts with no growth in the number of directions.

### 1.3 grid: free II.3 (referee: correct, one overstatement)

**Proved**
* **Normalization.** $\kappa_{\rm raw}=\sum_p(e_pV+L_p-U_p)\log p/(WV)$, with critical value $1/2$ and diagonal order
  exactly $V$.
* **Finite Vandermonde.** $\kappa=(\binom s2+\binom{n-s}2)/\binom n2$. The half cut gives $(n-2)/(2(n-1))$, and fair
  bits give $1/2$.
* **Cyclic six-mark design.** $\det=-\prod_{k=1}^5\sum_j\xi^{kj}z_j$ (the sign is $-1$). It is nonvanishing at angular
  width $<\pi/3$. Also $\kappa_{\rm raw}=1-W^{-1}\sum_p(\max-\min)\log p$, which is $\le\frac12-\log16/W$ on actual
  profiles.
  * **Correction:** $\kappa=0$ on *every* nonconstant cut, not only on three-three cuts. The fair expectation is $1/32$.
* **Every-total-degree design.** $V-1$ pair factors and one factor $u+v-s_{1,d-1}-s_{2,0}$. Fair-cut
  $\kappa=((V-1)/2+1/8)/V$.
* **Walsh-model table** for $d=1..4$. Raw: 61/155, 464/1085, 127/310, 1259/3255. Collected: 67/155, 513/1085,
  522/1085, 524/1085. These are model values only.
* **Lower-simplex obstruction.** The minimum over nonsingular minors of $\kappa$ is $\le\frac12-\log16/W$.
* **Fair-cut cap.** $E(m_S+m_{S^c})\le V/2$.
  * **Not new:** this is the DKSS multiplicity Schwartz–Zippel lemma, already in the record as threshold
    `construct.md` Thm B / Cor C3, identified as DKSS in round27 `kill_stab.md`. Astra does not cite it.
* **Laurent bookkeeping.** The forced exponent is $e_p(T+2Q)-\mathrm{width}_p$, so the design is supercritical iff
  $\sum\mathrm{width}_p\log p<vW/2$. This restates the record's mean-width criterion.

**Referee's synthesis, adopted.** Every design whose nonvanishing is free makes $F_G$ a product of linear forms with
zero coefficient sum. Each such factor contributes at most $W/2-\log16$ on an actual profile, so **no
free-nonvanishing design is supercritical.** A supercritical design needs an irreducible factor of degree $\ge2$, and
its nonvanishing at the actual point is exactly route 5's open zero estimate ZE$^{\rm all}_H$.

**Gap:** none.

### 1.4 index: free II.4 (referee: all correct)

**Proved**
* **Real affine independence of all rows**, from pair separation alone ($|x-y|\le x+y-2xy$, mean $\le1/2$). The fact
  is known (round5 `referee_construct.md`:233; round10 Lemma 3.2), but this proof is new and applies to gap-admissible
  angle models.
* **Growth-scale bounds.** $e_p\log p<4W/m$ (= round8 `stability.md` Thm 3.3(2)). Every five-tuple has a bounded
  minor $(2D)^4$ at $W\le bm\log m$.
* **Two free rows.** There are $O_q(m^3)$ rank losses mod a fixed odd $q$, with no hypothesis on $W$. The drifting
  arc centre $(h+k)\theta_0$ is handled correctly.
* **Parity.** No five-tuple loses rank mod 2.
* **Counting theorem.** On a normalized profile with $W\le bm\log m$, all but $O_b(m^3)$ ordered five-tuples have index
  one.
* **Baker–Wüstholz removal.** On genuine shortened clusters ($\Delta\le\frac12N^{-1/4}$), BW bounds the index
  absolutely, via Part I lemma §2.2 (refereed correct in `ref_baker.md`). So $O(m^3)$ holds with an absolute constant
  at any growth rate. The BW statement was confirmed only through secondary sources, the same caveat as Part I.
* **Index one suffices.** Index-one $(L_5^\circ)$ suffices for the uniform bound, and index-one (5) suffices for growth.
  See §3.2.
* **Calibration.** $H\le N^{31/33}$ suffices for $A_5(1/31)$. The typical $15/16<31/33$, with slack $1/528$.
* **Selection theorem.** Every genuine shortened one-unit cluster has a subcluster of $\ge m/145-K_Q$ points in which
  every five-tuple has index one. Every inequality was checked.
* **Paley double core.** Every proper restriction has index 1, and the whole profile has index $M/2$. Reproduced for
  $M\le32$.

**Computed:** the Gini certificate 315/256, 239/128, 633/256. It is not used by the final proof.

**Caveat (referee):** no genuine record cluster satisfies the shortening hypothesis; they are about 10× too long.
On them the index-one conclusion fails (indices 0, 3 and 4 occur on the 7-point record). The theorem is therefore
untestable on genuine data, though not contradicted by it.

**Gap:** none.

### 1.5 tri: free II.5 (referee: correct; not new)

**Proved**
* **Fake analogue.** In the angle model, near-Liouville edges have distinct labels $\prod p^{|v_p|}$, including prime
  powers. So fakes match $\gg m^2$ distinct labels. This is F1 confirmed for the divisor count.
* **Triangle labels.** The labels are $ab, bc, ca$ with $abc\mid N$ and $\gcd(a,b,c)=1$. Also $d_{12}d_{31}/d_{23}$ is
  a square.
* **(II.1)** $c=(au^2+bv^2\pm2suv)/w^2$ with $c\mid N/(ab)$, together with its reverse parametrization. Both were
  verified on 9473 genuine triples.
  * **Not new:** this is the binary-quadratic-form identity of round27 `stab.md` Prop 1.1(iii) (discriminant
    $-4t^2$). Astra does not cite it.
* **Profile condition.** The exact profile condition is the $k=3$ chain condition.
* **Vieta descent.** $C'=2(A+B)-C$. The descent is strictly decreasing, and the seeds satisfy $3AB\le4U$.
  * The descent changes the parent: the reflected triple $(1417,61,4553)$ is a genuine triangle at
    $393547661\ne131490079981$. So it cannot organize a fixed-parent count.
* **Example.** The worked example at $N=131490079981$ is reproduced.

**Correction (coordinator's question).** (II.1) is the three-point case of the *full Gaussian anchored data*
(ptolemy.md §1), not of Question A's data. It is precisely the norm condition $|P_2|^2=n_{12}n_2$, which Question A
forgets.

**Gap:** none. The "completion count" has no estimate. The only bound it yields is worse than the target.

### 1.6 height: lemma II.1-II.3 (referee: correct; one constant)

**Proved**
* **II.1 cutoff.** The sufficient cutoff is $H_S\le C^{-g}N^{g/4}$.
  * At index $\ge7$ this means $H_S\le C^{-94/25}N^{47/50}$.
  * At index $\ge9$ every $H_S\le N$ qualifies ($630>625$).
  * $H_S\asymp N$ does **not** suffice, which answers the coordinator's Q1 negatively as posed.
  * Squareclass dependence forces $|G_S|^2\le4C^4$.
  * Bounded-gcd squareclass example: $u_j=t+c_j+i$, $c=(0,1,3,7)$.
* **(II.2) concentration.** All but $O_{C,\varepsilon}(m^4)$ ordered five-tuples have $\log H_S/W$ within
  $\varepsilon$ of $15/16$. The mean $15/16$ is round23 `split.md` (ii); the concentration is new.
* **Corrected grid conclusion.** Index-$\ge7$ tuples number at most $(D+O_C(1))m^4$. At index $\ge9$ the bound applies
  to all tuples, and squareclass-dependent tuples number $O_C(m^4)$.
* **Growth-scale index count.** At $W\le Am\log m$, all but $O_{A,C}(m^4)$ tuples have index one.
  * **Correction:** the genuine-point gap is $\arcsin(e^{-w/2}/2)$, since $|u-1|=2|\sin|$. This is harmless for
    $\ell\ge2$.
  * **Superseded** by free II.4's $O(m^3)$, which needs no growth hypothesis.
* **Squarefree index.** Every full-rank squarefree index lies in $\{1,2,3\}$, from the maximal 0/1 determinant of
  order 4. This is classical and new to the record.
* **Covers.** Covers that are étale above $T$ are translated isogenies. A rational lift needs every prime column in the
  image sublattice, which excludes $I=1$ for unit translates. This also corrects the coordinator's
  ramification-divisibility argument, which holds only on one sheet.
* **Non-thinness.** Squarefree ten-pair profiles are not thin in $T^4$, by quantitative Hilbert irreducibility (Cohen
  1981; Serre, *Lectures on Mordell–Weil* §13.1). This was verified via secondary sources only, and it imposes no arc
  condition.
* **Cubic cover.** $v^3=J/b$ with $J=\prod t_{ij}/H^3$ lifts rationally when $\prod t$ is $b$ times a cube. The target
  has multiplicity 3, so the McKinnon–Roth ratio stays $210/3=70$: **no gain**.
* **Uniqueness.** $J^{k/3}$ is the only formal block-wise $k$-th-power pattern (Smith invariants $1^{14},3$).

**Gap:** none.

### 1.7 window: lemma II.4-II.5 (referee: correct; one exponent supplied)

**Proved**
* **Six-norm lemma.** If six integers in $[X,X+c_0X^{1/6}]$ have product $s^2+t^2$ with $t\ne0$, then
  $|t|\ge c_1\sqrt X$ (absolute constants).
  * [J] `judge_II/sixnorm.py` rechecks $\deg(F-Q^2)=2$ symbolically. A random sample ($64\le X<6000$, window
    $2X^{1/6}$) finds no representation with $|t|<2\sqrt X$.
  * The sharpness of the exponent $1/6$ is unresolved.
* **Application.** Via (II.5), $|u_i-u_j|=2|t_{ij}|/\sqrt{D_{ij}}$.
  * Exponent 4 ($\Delta\gg H^{-1/4}$) holds when the other four norms are $\asymp X$.
  * **Correction:** with no hypothesis on the other four norms, Astra's "pointwise gain" is, by the referee's
    averaging argument, the explicit uniform exponent **$g\ge17/5$**.
* **Congruence.** (II.6)-(II.7) give $b\le9T^4K^6+4T^6$ and the exponent $\min(1/6,(1-6\delta)/4)$. Only nine residues
  are needed.
* **Normalization and specialization.** The normalization (monic) and the specialization to an exponent-$10/3$
  template with bounded gcd are a routine formalization of the record's (B4).
* **Literal wording.** The wording permits disjoint-edge conjugate overlap and excludes the other overlaps.
* **(II.8) moments.** At $d=1$ this is the record's `growth/cotangent.md` statement (450 LM runs, all degenerate),
  cited correctly by Astra.
* **No pencil.** No polynomial norm pencil exists in any degree: an identity $n_e=F_e(G)$ forces every $F_e$ to have no
  real zeros.
* **Reconstruction from four base blocks.** $C_i=U_i+h_iD_i$, with the gcd and disjointness conditions necessary and
  sufficient.
* **Generalized circles.** Linear blocks cannot have all ten roots on a generalized circle.

**Computed:** the mod-7 exclusion. There is no good-reduction $\mathbb F_{49}$ solution of the first five moment
equations (95,760 stars, 0 solutions).
* *Weight:* this is exactly the random expectation ($1.4\cdot10^{-4}$ expected; the KMAX = 3 control finds 1638
  against a predicted 789).
* So it is **not** evidence of rigidity. It shows only that a $\mathbb Q(i)$-rational $d=1$ solution has bad
  reduction at 7.

**Gap:** none. The full polynomial ten-pair system is open even at $d=1$.

### 1.8 plane: lemma II.6-II.7 (referee: correct; three cosmetic corrections)

**Proved** (one-unit-class ten-pair model, pairwise coprime block norms, norm vector in a fixed rational subspace $U$
of dimension $\le2$)
* **Reduction to Pell.** With $U$ and the signed $t$ fixed, a sextic constant on infinitely many rays is handled as
  follows:
  * factors of degree $\ge3$ go by Thue;
  * a linear factor contradicts (II.10);
  * two non-proportional quadratics are impossible;
  * definite quadratics give finitely many;
  * what remains is $F=\kappa Q^3$ with $Q$ indefinite and irreducible over $\mathbb Q$.
* **Conic step.** $n_e=a_eu+b_e/u$. The real $y_e=\pm\sqrt{b_e/a_e}$ have constant star sums and constant reciprocal
  star sums.
* **Real five-star lemma.** There are no nonzero real $y$ on the edges of $K_5$, distinct up to sign, with constant
  star sums and constant reciprocal star sums. So the Pell conic is excluded.
  * **Correction:** the critical-value determinant is $-b^2(9b+16)(24b+1)/16$ (sign corrected; same zeros).
  * Reality is essential, since nondegenerate complex solutions exist.
* **II.7 cases.** Proportional forms give $B\ll_UT^6$. In the Plücker step, either $B\ll T^2$ or $M=u\wedge(pv+qw)$.
  The transfer of (8) to (P3) gives $B\ll T^4$ unless $V=0$; both steps Astra flagged were checked symbolically.
* **Genus-five curve.** It is a smooth complete intersection of three diagonal quadrics in $\mathbb P^4$; the genus
  was recomputed two ways. Faltings then gives finitely many projective cofactor patterns, depending only on $U$. This
  is where ineffectivity enters.
* **Roth/Liouville size step.** For a fixed pattern, $|\lambda|\gg B^{1/2}$. Roth is needed only for squared cubic
  factors and simple roots of degree $\ge3$; Liouville covers the quadratic cases.
* **Theorem.** $T=\max|t_e|\gg_U(\max n_e)^{1/6}$, and $\gg_U(\max n_e)^{1/4}$ when the ten forms on $U$ are pairwise
  non-proportional. The constants are uniform in the tuple, depend on $U$, and are ineffective only in the stratum
  where all Plücker forms vanish and $V=0$.
* **Exponents.** The approximation exponent is $\ge38/11$ on fixed planes, $40/11$ for comparable norms, and $H^{-2/7}$
  in the proportional case.
  * **Correction:** this needs *private blocks*, so that $H=\prod n_e$ is the primitive norm. That is automatic in
    Part I's model but not restated in the II.7 headline.
* **Arithmetic progressions.** These give exponent 4 in every edge assignment (10,080 anchor triples and 9! slope
  assignments, 0 degenerate configurations). The constants there are absolute and effective.
* **Prime test case.** $34913+2100r$ ($0\le r\le9$) are primes $\equiv1\pmod4$, giving a squarefree index-one test case.
* **Green–Tao–Ziegler** gives arbitrarily large comparable prime families $4x+4ry+1$. It is used only for
  nonemptiness.

**Gap:** none.

### 1.9 Withdrawn or corrected in Part II
* free II.3: "$\kappa=0$ on three-three cuts" is strengthened to "on every nonconstant cut".
* free II.3: the fair-cut cap is credited to DKSS and to the record's Thm B / C3.
* free II.5: (II.1) is credited to round27 `stab.md` Prop 1.1(iii).
* free II.1: the "28 pairs" sentence is corrected. The coordinator's "weight $D/2$ is exactly the minimum" is
  withdrawn (it was the coordinator's slip, not Astra's or the record's).
* free II.2: "for each fixed $B$" holds only beyond astronomically large, $B$-dependent parents.
* lemma II.2: the genuine-gap constant is corrected by a factor 1/2. The $O(m^4)$ count is superseded by $O(m^3)$.
* lemma II.4: the unquantified "pointwise gain" becomes $g\ge17/5$.
* lemma II.6: the determinant sign is corrected.
* lemma II.7: the private-block hypothesis is added.
* Part I §2's reduction of $(L_5^\circ)$ to indices $\{1,3,5\}$ is not wrong, but it is **superseded** for clusters by
  free II.4.

## 2. What is genuinely new relative to the record, and correct

The judge checked these against the record by grep, and the referees did the same.

1. **The five-tuple index theory of clusters** (free II.4): the $O(m^3)$ count, the two-free-rows argument, the
   selection theorem $m/145$, and the Paley double-core observation. grep finds no prior five-tuple index discussion
   apart from Part I. This is the round's most consequential clarification.
2. **The fixed-plane rigidity theorem** (lemma II.6-II.7): the real five-star lemma, the genus-five cofactor curve,
   and the arithmetic-progression exponent 4. This is new and mathematically the strongest piece. See §3.3 for its
   F1 status.
3. **The six-norm short-interval lemma** and its congruence variant (lemma II.4). New and elementary.
4. **The exact $D=8$ minimal-vector count** (12, cube-only), with the first-gap and minimum-14 facts. The record
   listed this as not obtained (round30 `prove.md` §5).
5. **Closed negatives:**
   * free-nonvanishing grid determinants are subcritical (free II.3);
   * covers étale above $T$ cannot lift index-one points, and the cubic cover gives no volume gain (lemma II.3);
   * no polynomial norm pencils (lemma II.5);
   * the $H_S$ cutoff $N^{47/50}$ (lemma II.1).
6. **Finite templates C1/C2** with exact CRT progressions; the ellipse and substitution obstructions.

Not new, or routine:
* the fair-cut cap (DKSS = threshold Thm B / C3);
* (II.1) of free II.5 (round27 `stab.md`);
* (II.8) at $d=1$ (`growth/cotangent.md`);
* the mean $15/16$ (round23 `split.md`);
* affine independence (round5, round10);
* $e_p\log p<4W/m$ (round8);
* the Laurent mean-width criterion;
* the mod-2 Sidon fact (Part I).

## 3. What changes for the integer problem

### 3.1 Growth rate, uniform bound, Question A, route 5
* **Growth rate:** unchanged at $(2/3+\varepsilon)\log R/\log\log R$.
* **Uniform bound:** open and unchanged. No lower bound for $M_\infty(C)$ moves; the C1/C2 templates have contraction
  $\sim10^{10}$.
* **Question A ($L_8$):** untouched. Free II.5 only identifies (II.1) as exactly the norm condition Question A forgets.
* **Route 5:** untouched.
* **The grid-determinant idea (coordinator Q4) is answered negatively for every free-nonvanishing design.** Such
  designs factor into zero-sum linear forms, each of which is subcritical on actual profiles. Any supercritical grid
  design needs an irreducible factor of degree $\ge2$, whose nonvanishing at the actual point *is* an instance of
  ZE$^{\rm all}_H$. The grid re-poses route 5's missing zero estimate; it does not bypass it. **Record this as
  closed.**

### 3.2 Question B: is "index-one $(L_5^\circ)$ suffices" established?
**Yes.** The judge read free II.4 lines 545-605 and the referee's §9 chain, and rules it established.
* The hypothesis is exactly the record's $(L_5^\circ)$: fixed proper closed $Z\subset T_5$; absolute $\eta,\gamma,w_0$;
  one unit class; primitive; all fifteen cuts in $[(1-\eta)w,(1+\eta)w]$; conclusion $H_\Sigma\ge\gamma w$. The only
  change is that the hypothesis is restricted to tuples of index one.
* The chain:
  1. Subdivide the arc (a constant factor in $C$) and keep one unit class.
  2. Part I §5 (refereed) gives $\gg m^5$ balanced tuples with $H_\Sigma=O(W/m)$.
  3. Remove the higher-index tuples: there are $O(m^3)$ of them, with an *absolute* constant, by BW. This is what
     uniformity needs, and it involves no fixed-$S$ exceptional set.
  4. Remove the zeros of $Z$ by Schwartz–Zippel ($\le\deg\cdot m^4$).
  5. Obtain the contradiction.
* The standing cautions are respected: the gcd is divided out, with $t_{ij}$ and $d_{ij}$ unchanged; the common phase
  drift is carried through the mod-$q$ step, not set to zero; and no threshold depends on a template.
* The citation caveat is the same as Part I: the BW two-logarithm statement was confirmed via secondary sources only.
  Any Baker-type bound of shape $\ll h(u)\log n$ suffices.
* **Substance.** Index-one $(L_5^\circ)$ is formally weaker but not easier. Every hard model lives there: the Paley
  fakes, the ten-pair-cut angle model at $10/3$, and the balanced squarefree profiles. The precise statement of
  Question B for clusters is now **the index-one case of $(L_5^\circ)$**. Part I's index-$\ge7$ cover theorem and the
  $\{1,3,5\}$ reduction contribute nothing to the uniform bound or to growth.

### 3.3 Are the fixed-plane and short-interval gains the record's first unconditional gains over $10/3$ on genuine index-one tuples?
**Literally, no** [J]. Two earlier results already beat $10/3$ on genuine index-one tuples:
* Balanced profiles get $15/4$ from Liouville with toric sections (round 23 `width.md`; expert note line 286). Balanced
  squarefree 15-cut profiles have index 1 (Part I §4).
* Part I (7) gives $7/2$ on all squarefree primitive norms with $\le9$ primes, and $4$ with $\le6$ primes.

Both are *gap-only* bounds, which fakes satisfy.

**In the sense that matters, yes, with qualifications** [J, unrefereed]:
* Part I §6 (refereed) shows that on the **ten-pair-cut profile** the gap-only angle model attains exactly $10/3$.
* That construction uses the ten block norms only through their sizes (primes $p_E\in[X,2X]$, with free angles). The
  judge reads it as applying verbatim on any fixed plane with comparable norms, for example the arithmetic-progression
  plane $p+rq$ with $q\le p/9$ (Green–Tao supplies the primes). Not refereed.
* So on these families gap-admissible fakes reach $\Delta\asymp H^{-3/10}$. Genuine tuples have
  $\Delta\gg_UH^{-11/40}$ (comparable norms; exponent 4 on the arithmetic-progression plane with absolute constants),
  or $\Delta\gg H^{-1/4}$ under the six-norm window.
* As far as the judge's grep of the record shows, these are the **first unconditional approximation-exponent bounds
  that beat the fake value on a nonempty family of genuine index-one tuples**. They use integrality of the cofactors
  $t_{ij}$, identity (8), and the shared Gaussian blocks, none of which fakes have.
* For the window family the fake side is heuristic only, since six primes in $[X,X+cX^{1/6}]$ are not known to exist
  unconditionally.

**Why this does not move Question B:**
1. **Wrong regime.** The fixed-plane hypothesis is a codimension-8 linear condition on the ten norms, and nothing in
   the growth regime supplies it. Every tuple lies on *some* plane of height $\approx X$, so only sequences on one
   fixed plane are constrained.
2. **Wrong uniformity.** The constants depend on $U$ and are ineffective (Faltings, Roth) in the hard stratum. A plane
   of height $X^{o(1)}$ is not covered.
3. **Wrong profile.** The ten-pair profile has no singleton cuts, so it lies outside the balanced $(L_5)$ domain.
4. **Too weak.** $40/11=3.636<15/4$, and $A_5(\delta)$ with $\delta<1/30$ needs more than $15/4$.

The short-interval window has relative width $X^{-5/6}$, which nothing in the problem provides.

**Verdict:** a real conceptual first (genuine beats fake, on a thin family) with zero effect on the rates.

### 3.4 F1 and F2
* **F1:** confirmed on the divisor side (free II.5: fakes match $\gg m^2$ distinct labels). It is escaped only on the
  thin families of §3.3.
* **F2:** route 5 supercriticality is untouched. The grid analysis confirms that every *free* design is exactly or
  sub-critical.

### 3.5 Suggested additions to `expert_note/expert_note.md` (recommendations only; the record is read-only)
1. **Question B:**
   * *"For shortened genuine one-unit clusters, all but $O(m^3)$ five-tuples have index one (absolute constant, via
     Baker–Wüstholz), and a subcluster of $\ge m/145-O(1)$ points has every five-tuple of index one (round 31). Hence
     $(L_5^\circ)$ restricted to index-one tuples suffices for the uniform bound."*
   * Downgrade the index-$\ge7$ cover theorem to "irrelevant to clusters".
2. **Approximation exponents:**
   * *"On the ten-pair-cut profile, where gap-only inputs give exactly $10/3$, genuine tuples with block norms on a
     fixed rational plane $U$ have exponent $\ge38/11$ ($40/11$ for comparable norms; $4$ on generic planes, including
     arithmetic progressions), with constants depending ineffectively on $U$ (round 31)."*
3. **Function-field model:** *"$D=8$: at most 12 constant-offset divisors up to sign, cube only. The lattice minimum
   can exceed $D/2$."*
4. **Grid determinants:** list them as closed. Free nonvanishing implies subcritical. Credit the cap to DKSS / Thm B.

## 4. Open sub-questions, ranked by honest probability of a growth-rate improvement

All probabilities are the judge's estimates [J].

1. **Fixed-plane rigidity without a fixed plane.** Extend the star-sum collapse to planes $U$ of height $X^{o(1)}$, or
   to comparable ten-pair norms with no plane hypothesis, or to the 15-cut profile with singleton blocks (where (8) and
   Question B actually live). It is the only F1-escaping exponent mechanism in the record. P ≈ 0.01.
   * *First step:* make the $V\ne0$ branch of II.7 ($B\ll T^4$) effective, with explicit dependence on the height of
     $U$, and test whether it survives at height $X^{\varepsilon}$.
2. **Index-one $(L_5^\circ)$ via identity (8) and cofactor integrality on balanced profiles.** This is the
   reformulated core of Question B. P ≈ 0.008.
   * *First step:* write the analogue of (8) for the fifteen-cut profile, and ask whether a power lower bound on
     $\max|t|$ follows under a comparability hypothesis alone, as in II.4's window.
3. **Decide the $d=1$ ten-pair polynomial system rigorously**, via the four-base-block reconstruction: six cubic
   resultants modulo the real affine group, saturated against gcd and collisions. A solution would be a genuine
   $10/3$ template, which is negative information. Non-existence would be the first rigidity at index one in the
   function-field model. P ≈ 0.003.
   * *First step:* Groebner basis or certified homotopy on the reconstruction system (needs Singular/M2/msolve, not
     installed here).
4. **Supercritical grid designs with an irreducible factor of degree $\ge2$**, posed explicitly as instances of
   ZE$^{\rm all}_H$ on sparse cuts. P ≈ 0.003.
   * *First step:* find the smallest design whose forced width is supercritical on a model sparse profile, and state
     its zero estimate.
5. **The common-parent completion count** (free II.5): bound
   $\sum\mathbf 1[(au^2+bv^2\pm2suv)/w^2\mid N/d]$ by $o((W/\log W)^3)$ using the integer $s$. P ≈ 0.002.
   * *First step:* compute the count on the record's parents to see whether it is even plausibly $o$.
6. **The cubic cover with a finer local constant** at its multiplicity-3 target point. This would need an
   approximation constant at a singular point beyond McKinnon–Roth Cor 4.2. P ≈ 0.002.
   * *First step:* compute the Seshadri constant and the volume of the blow-up at the singular point for the cubic
     hypersurface.
7. **Minimal-vector counts at $D=10,12$ that are sublinear in $D$**, and templates with a growing number of
   specializable directions under controlled CRT moduli. P ≈ 0.001.
   * *First step:* run the `ref2_min/flats.cpp` census at $D=10$, ternary case only.
8. **Close the independent-verification gap** on the 7 planar shapes of the repeated-factor $D=8$ variant (Singular or
   M2). Bookkeeping only. P < 0.001.

## 5. Scores (0-10)

| block | score | reason |
|---|---|---|
| min (free II.1) | 3/10 | Exact, complete answer to Q2(a)/(b) at $D=8$, plus two slips corrected. Says nothing for growing $D$. |
| template (free II.2) | 2/10 | Correct finite constructions with exact CRT progressions; fixed counts; astronomically non-uniform in $B$. |
| grid (free II.3) | 2/10 | A clean closed negative to Q4. The key cap is the record's DKSS Thm B, uncited. |
| index (free II.4) | 4/10 | The most consequential clarification: Question B for clusters is exactly index-one $(L_5^\circ)$, with absolute constants. Formally sharp, substantively a relabelling of the hard case. |
| tri (free II.5) | 1/10 | Correct, but its main identity is already in the record (uncited), and the descent changes the parent. |
| height (lemma II.1-II.3) | 2/10 | Correct clarifications (the $N^{47/50}$ cutoff, no lift at index one, cubic cover with no gain); the index count is superseded. |
| window (lemma II.4-II.5) | 3/10 | A new, elementary six-norm lemma (exponent 4, uniformly 17/5). The polynomial test stays open; the mod-7 evidence is empty. |
| plane (lemma II.6-II.7) | 4/10 | High mathematical quality, and the first genuine-beats-fake exponent on index-one tuples [J]. Thin, non-uniform, ineffective, below $15/4$. |
| **Part II, on the goal** | **2/10** | Error-free; three coordinator questions closed; Question B sharpened to index one; one conceptual first. No change to the growth rate, the uniform bound, Questions A/B or route 5. |
