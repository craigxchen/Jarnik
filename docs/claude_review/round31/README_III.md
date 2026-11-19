# Round 31, Part III -- judge report: Astra's third replies

Sources: Part III of `astra_free_III/answer.md` (heading "# Part III", §§III.1-III.11, scripts `iii_*`) and Part III
of `astra_lemma_III/answer.md` (§§III.1-III.10, scripts `p3_*`). The coordinator's questions are in `reply3.txt` in
each directory. Eight referees each took one block:

| key | report | Astra sections |
|---|---|---|
| density | `ref3_density.md` | free III.9 (and the parts of III.1, III.7 it uses) |
| quartet | `ref3_quartet.md` | free III.1-III.3 |
| rank | `ref3_rank.md` | free III.4-III.5 |
| model | `ref3_model.md` | free III.6, III.11 |
| proj | `ref3_proj.md` | free III.7, III.8, III.10 |
| four | `ref3_four.md` | lemma III.1-III.4 (plus summary items citing III.9) |
| linear | `ref3_linear.md` | lemma III.5, III.7, III.8, III.9, III.10 (barred list) |
| sweep | `ref3_sweep.md` | lemma III.6, III.10 (template degree) |

Judge scripts are in `judge_III/` (`thresholds.py`, `sector.py`, `sector2.py`). The Part I judge report is
`README.md`; no Part II judge report (`README_II.md`) exists at the time of writing, so nothing here relies on one.
**Every claim labelled [J] is the judge's own and is unrefereed.**

Labels. "Proved" means a complete proof checked line by line by a referee, or an exact computation a referee
reproduced with independent code. "Computed" means exact data reproduced independently. "Unreproduced" means Astra's
computation that no referee reran. No referee found a false theorem. One referee (density) reported FALSE on the
*non-triviality* of a headline (§1.1, ruled on below), and one (linear) recorded the still-open decision of the strict
linear system as GAP; that is an open problem, not a gap in a claimed proof (§1.7).

## 0. Verdict in five lines

1. **Nothing changes for the integer problem.** The growth rate stays $(2/3+\varepsilon)\log R/\log\log R$. There is
   no uniform bound, and Questions A and B and route 5 keep their status. Astra says the same (summary item 11).
2. **The density theorem III.9 is proved** (with a one-word repair: the sample must be in uniformly random order). It
   is an *almost-all* statement in a density-zero family of parents. **Its headline numbers are not meaningful:** an
   elementary pair first-moment count (density referee, rechecked by the judge) gives $\le2$ points at $D=4$ and $\le1$
   for $D>4$, against Astra's 17 and 13. Quartet exposure adds something only for $1<D\lesssim2.08$ (e.g. 37 at $D=3/2$).
   It does not see the worst parent and does not separate genuine circles from fakes in any useful sense (§3.2).
3. **The integer countermodel III.11 is correct**: integrality of a primitive $2\times5$ star matrix, *all* ten
   divisibilities, bounded cofactors, index one and norms to relative error $O(N^{-1/10})$ give width $\asymp
   N^{-3/10}$, so they cannot beat the five-point exponent 10/3. It is projectively a point of the record's C1 family
   and says nothing about Question A. III.6 (growing $m$) keeps only half the log divisibility, the trivial CRT ceiling.
4. **The rest is route closure, all correct:** unit-equation rank must be linear in the extension count (route (a)
   closed); fixed-cofactor elimination is dominant (closed); quartet cross-ratio consistency is automatic (closed);
   quartet recovery and global projective rigidity hold only for $m\gtrsim10^4$; the strict linear-template system
   is still undecided, with tangent dimension $\le3$.
5. **Score for Part III on the goal: 1/10.** Competent, almost error-free, and almost entirely peripheral.

## 1. Referee-corrected statements per block

### 1.1 density: free III.9 (referee: proof correct; headline superseded; FALSE on non-triviality)

**Proved**
* **III.9, corrected statement.** Fix $D>1$ and $C>0$. Draw an *ordered* sequence of $r$ distinct primes uniformly
  from the split primes in $[r^D,2r^D]$ (equivalently, put the drawn set in uniformly random order) and let $N$ be their
  product. There are $M_D$ and $c=c(D,M_D,C)>0$ such that for $r\ge r_0(D,C)$, outside a fraction
  $\exp(-c\,r\log r)=N^{-c/D+o(1)}$ of parents, no arc of length $C\sqrt R$ carries more than $M_D$ points of a normalised
  cluster with normalised norm $N$. Here (III.7), $(D-1)F_M>3D/(2(M-1))$ with $F_M=\frac{(M-1)(M-4)}{8(M-2)(M-3)}$,
  gives $M_D=17$ at $D=4$, $13$ for large $D$, $37$ at $D=3/2$, and a finite value for every $D>1$. The exceptional set
  depends on $C$ (and on $D,M$). It is not uniform in $C$.
* **Coverage (corrected).** Covered: $x^2+y^2=N$ for $N$ in the good set; and $x^2+y^2=n$ with $n=N|G|^2$, for clusters
  all of whose points are divisible by one Gaussian integer $G$ (powers of 2, inert squares, a split factor common to
  all points), at the same $C$. **Not covered:** any cluster on $n=Nq$, $Nq^2,\dots$ in which an extra split prime $q$ is
  a non-constant cut; any $n$ whose normalised norm is not a family member. The family has density zero.
* **Bounded-exponent extension.** Exponents $1\le e_p\le E$ chosen adversarially after the primes are drawn, under
  $(D-1)F_M/E>3D/(2(M-1))$. Proved, on the random support only.
* **Pair first-moment bound (referee's, §9 of `ref3_density.md`; checked by the judge).** With
  $\rho_M=\lfloor M^2/4\rfloor/\binom M2$: if $D>4\rho_M$ then, outside a fraction $\exp(-c\,r\log r)$, no arc of length
  $C\sqrt R$ carries $M$ points. Hence **$D>4$: at most 1 point; $8/3<D\le4$: at most 2; every $D>2$: $O_D(1)$.**

**Computed.** All (III.7) thresholds, by exact rationals (referee; [J] reproduced in `judge_III/thresholds.py`:
$D=4$: $0.371875>0.352941$ at $M=18$, $0.371429<0.375$ at $M=17$; maxima 133, 37, 25, 17, 13 at
$D=1.1,1.5,2,4,10^6$; pair-count maxima 20, 10, 4, 2, 2, 1 at $D=2.1,2.2,2.5,3,4,5$). Exhaustive enumeration of all
$4\cdot2^r$ points on random parents ($r=10..18$): maximum 1 point per arc except a single 2 at $r=10$, $D=3$, $C=100$.

**Withdrawn / corrected**
* "Order the sampled primes" read as increasing order: the conditional law is then not uniform and $1/\binom{K-r+s}{s}$
  is unjustified. Replaced by a uniformly random order.
* The headline "17 at $D\ge4$, 13 for large $D$" as the content of the theorem: superseded by 2 and 1.

**Judge's ruling on the FALSE verdict (non-triviality).** I read III.9 and §9 of the referee report. Astra does not
claim non-triviality in words, so nothing she wrote is false; the referee's verdict is about the meaning of the
headline, and on that I uphold it. I checked the pair-count proof line by line [J]:
(i) for two points with separating set $S$, $|S|=k$, $\omega=P_i\bar P_j/g_{ij}$ has norm $d=\prod_Sp$ and
$z_i\bar z_j=\varepsilon N\omega/\bar\omega$, so $2\arg\omega\in\frac\pi2\mathbb Z+[-\Delta,\Delta]$;
(ii) $\omega$ is neither real nor diagonal, since $d$ is odd, squarefree and $>1$;
(iii) the number of such $\omega$ with $|\omega|^2\le X$ in the eight thin sectors is $\le c\,X\Delta$ with no boundary
term, because the off-axis coordinate is a nonzero integer ([J] `sector.py`: count $\le0.50\cdot8\tan\delta\,(X+\sqrt X)$
over 60 random $(X,\delta)$, and zero once $X\Delta^2<4$);
(iv) each $\omega$ fixes the prime set of $S$, and after exposing the other positions (random order) the $S$-positions
carry a uniform $k$-subset of $K-r+k$ primes, so the probability is $\le cC(2P)^kP^{-r/4}/\binom{K-r+k}{k}\le
cC(6kD\log r)^kr^{-Dr/4}$;
(v) $\sum_{i<j}k_{ij}=\sum_pu_p\le r\lfloor M^2/4\rfloor$, so some pair has $k\le\rho_Mr$; summing over $2^r$ position
sets gives $\exp(-(D/4-\rho_M-o(1))r\log r)$. Correct. **So for $D\ge4$ the statement of III.9 is weaker than an
elementary fact, and the only range in which the quartet argument is needed is $1<D\le2$.** (Just above $D=2$ the
pair bound is very large, so the quartet bound is still better on $2<D\lesssim2.08$; at $D=2.1$ the pair count already
wins, 20 against 24.)

### 1.2 quartet: free III.1-III.3 (referee: no error; framing corrected)

**Proved**
* **(III.1)** The three-product Plücker identity $Ax-By+Cz=0$ with $A,B,C$ pairwise coprime (surviving exponent
  $h_3-h_2$ on the matching of the two lowest values), $ABC\mid N$, the lcm formula from the $d_{ij}$, and
  $|t_{ij}|\le\sqrt{d_{ij}}\,\Delta/2$, under the one-unit convention. *Not new*: this is `growth/ptolemy.md` Thm A(2),(4)
  and Thm B. Only the measured-divisor lcm formula is new, and minor.
* **Recovery lemma (III.2).** If $2U^2\max(u,v)<\gcd(B,B')$ then $uxz'=vx'z$; with coefficients coprime to the parent,
  $A=A'$, $C=C'$. New.
* **III.2 assertion 1.** On the narrow-bin Paley parent, a random prime assignment passes even one quartet with
  probability $O(m^42^{-cm})$, effectively $c=1/32$; the bound is below 1 only for $m$ in the low thousands.
* **III.2 assertion 2, corrected wording.** *Admissible assignments differ pairwise on more than $\lfloor r/64\rfloor$
  primes, so at most $r!/\lfloor r/64\rfloor!$ assignments admit an integer star matrix with the prescribed
  divisibilities and cofactor bounds.* This bounds a code; it does not show the code is empty.
* **III.2 assertion 3.** At most one of the $2^m$ within-pair orderings is admissible.
* Genuine clusters on the narrow-bin Paley profile would be admissible ($|tt'|\le\lfloor P/4\rfloor$), for one unit
  class and $\Delta\le\frac12N^{-1/4}$.
* Every assignment on narrow bins admits a gap-admissible angle fake. *Not new*: `round16/barrier.md` Thm 2.3.
* **III.3.** $\operatorname{lcm}_kQ_k\mid B(ij|xy)$; eliminating a common mark returns only the direct quartet
  congruence (squarefree case proved; prime-power case: 0 violations in 20,000 random configurations, **evidence**).
  A negative result.

**Computed** (reproduced with independent code). $m=7,11,19$: 0 surviving quartets (0/35, 0/330, 0/3876). First parent
$100049..100333$: 156/45045 pass for $k=2$ (16/12 factorable), 0/560560 for $k=3$. Planted fourteen-prime parent:
exactly 1 of 560,560 partitions passes ($-1145443A+1927968B-782525C=0$, $U=62500076862$), and no five marks pass under
any assignment. Planted at $m=7$; no asymptotic content.

**Withdrawn / corrected.** The opening "angle fakes fail the integer test" holds vacuously (angle fakes carry no integer
star matrix); it is not evidence of an F1-violating input.

### 1.3 rank: free III.4-III.5 (referee: no error; scope stated)

**Proved**
* **III.4 exact edge version.** If $s$ third vertices of a fixed edge $ij$ have divisor splits in one rank-$r_0$
  multiplicative coset, then $s\le8r_0+4$. No size conditions. Uses only pair separation, **so it holds in angle fakes**.
* **(III.4), corrected statement.** *Three anchors with cofactors $\le m^B$ and $(F_x,1-F_x)$ in a rank-$r_0$ coset give
  $s\le(4/D)r_0+O_B(\sqrt{r_0})+O_B(1)=(16/3+o(1))r_0+O_B(1)$*; the additive constant depends on $B$.
* One rational coordinate: $s\le(8+o_B(1))r_0+O(1)$, also with varying small multipliers and moving offsets.
* **Consequence.** At a growth failure the needed unit-equation rank is $\gg m$; ESS
  ($\exp((6n)^{3n}(r+1))$, Ann. Math. 155 (2002)) is then useless. ESS is a benchmark only, not used in a proof. This
  **closes the coordinator's route (a)** for anchored cross-ratio and triangle equations.
* **III.5.** Factorial density threshold; at $W=3m\log m$ some cofactor colour has $m^{5/4-o(1)}$ edges, forcing $C_{10}$
  but not $C_4$ (Bondy-Simonovits); excluding monochromatic $C_{10}$ would force $W\ge\frac{16}{5}m\log m$ (constant
  5/8). Purely hypothetical: nothing supports a $C_{10}$ exclusion.
* Parity: fixed-cofactor graphs are bipartite; $c=(a+b\pm2s)/h^2\mid N/d$ with $a\equiv\mp s\pmod h$.
* **Balanced four-cycle.** Genuine one-unit primitive 4-point clusters for arbitrarily large $N$ with all 7 cuts of
  weight $(\log N)/7+O(1)$, width $\asymp N^{-2/7}$, cofactor 2 on a $K_{2,2}$. Balance is only asymptotic: the $O(1)$ is
  18-31 nats at the certified instance. It refutes a 4-point $C_4$-exclusion lemma.

**Computed.** Paley $m=7,11,19$: projected anchor rows have full affine rank $m-4$ for all 35/165/969 triples.

**Not checked by the referee:** the twisted-height translation, the higher-word Walsh argument and the Auerbach
$\sqrt m$ bound in the linked notes.

### 1.4 model: free III.6, III.11 (referee: no error; scope corrected)

**Proved**
* **III.11, corrected statement.** *There are arbitrarily large split parents with ten pair cuts of equal weight up to
  $o(1)$, valuation index one, and one primitive integer $2\times5$ star matrix such that every $g_{ij}$ divides its
  minor with cofactor exactly $\pm(u_j-u_i)\in[40,520]$, $|V_i|^2=n_i(1+c_iK^{-2}+\dots)$ with fixed $c_i\ne0$ (relative
  error $O(N^{-1/10})$, absolute defect $\asymp n_i^{4/5}$), and $\operatorname{diam}\{2\arg V_i\}\asymp N^{-3/10}$.* At the
  distinguished primes, 0 of 25 marked incidences are isotropic (not "at least one non-isotropic line"). Two instances
  checked exactly (17 and 143 digits); the kernel generator equals Astra's $(f,h)$; Smith form $(1,1,1,1)$;
  $\operatorname{diam}\cdot N^{3/10}=1040$, minimum separation $\cdot N^{3/10}=80$.
* **Why it stops.** Pair-cut-only profiles give strict half separation only for $m=5,6$ ($m^2-9m+16<0$). The six-mark
  sheared example has all 15 divisibilities, radii exact to relative $2\cdot10^{-34}$ at the instance, width
  $\asymp N^{-4/15}$.
* **Exact norms force genuineness.** Primitive columns with $|V_i|^2=n_i$ and all minor divisibilities are
  $u_i\prod_p\pi_p^{a_{ip}}$, i.e. a genuine cluster. (The record's star/point dictionary in reverse.)
* **III.6.** For fixed $L$ and growing $m$: $m$ primitive columns on a Paley parent with $W=O_L(m\log m)$, width
  $\le\frac12N^{-1/4}$, norms and normalised minors within $m^{-L}$ of the targets, exact valuations on a set $S$, and
  $\log g_{ij,S}=\frac12\log g_{ij}+O((\log m)^2)$. Correct but nearly trivial: half the log divisibility is the CRT
  ceiling for columns $(X,1)$. Assertion 4 (no quartet admits full divisibilities with cofactors $\le3m^3+1$ at a
  suitable assignment) follows from III.2 assertion 1, which the quartet referee checked; so it stands.
* Bounded mean log slack in a fully isotropic integral model forces genuine points on a positive fraction of marks.

**Computed.** III.6 at $m=19,43,67$ ($P=10^{15}$ at $m=67$): exact valuations on $S$, $g_S$ divides every minor, full
$g_{ij}$ divides none (0/171, 0/903, 0/2211). Pre-asymptotic for small $m$, as the requirement $P\ge m^{4L+60}$ predicts.

**Withdrawn / corrected.** "Asymptotically exact star radii" must carry the rate: relative $O(N^{-1/10})$, far coarser
than the angular scale $N^{-3/10}$.

### 1.5 proj: free III.7, III.8, III.10 (referee: no error; thresholds made explicit)

**Proved**
* **III.7 local uniqueness.** $B>2U^2$ forces a unique $[x:y:z]$, so the quartet cross-ratio is determined by the six
  $d_{ij}$ and the cofactor bounds.
* **III.7 main theorem, with the explicit threshold.** *Under pair separation, for $m\ge16{,}390$ (Astra's constants;
  about $10^4$ with sharp counting), the parent and the labelled $d_{ij}$ determine all rational star directions up to
  one element of $\mathrm{PGL}_2(\mathbb Q)$.* Astra never states the number.
* Robustness: at most $m/10$ deleted partners per mark (needs $m\gtrsim4.1\cdot10^4$), or $o(m^2)$ deletions with
  uniqueness on $m-o(m)$ marks.
* Almost every parent from a pool of $K\ge m^{1+\eta}$ primes excludes every assignment on every Paley quartet
  (precursor of III.9; Paley-specific).
* **III.8 (III.6)** $R_{ij;kl}R_{ij;lh}=R_{ij;kh}$ when $Q>2U_1U_2U_3$; five-set consistency; gluing over all regular
  quartets for $m\gtrsim5\cdot10^5$. The gluing write-up is sketchy; the referee's induction (quartets meeting the chart
  triple in 2, 1, 0 marks) repairs it. A no-go: cross-ratio consistency adds nothing beyond single quartet tests.
* **III.10.** $t_{ij}=c\,s_is_jh_{ij}$; integral realisation by saturation; the primitive-column criterion; one quadratic
  form fixed by three directions. A reformulation of the residual problem, with no progress on it.

**Gap in scope (referee, accepted).** On the record's 8-point clusters (golden $n=4..22$), reconstruction from the
$d_{ij}$ alone determines at most one mark beyond the three anchors, and those clusters violate the pair-separation
hypothesis ($\min\log d/W=0.36$-$0.47$). III.7 says nothing at $m=8$.

### 1.6 four: lemma III.1-III.4 (referee: no error)

**Proved**
* III.1 bookkeeping corrections to Part I: faithful.
* **III.2.** For fixed signed $t$ the continuous norm vectors form an irreducible fourfold of degree 12 (4 linear
  equations, the Plücker quadric, one sextic). Degree confirmed by an exact rational slice (squarefree resultant of
  degree 12) and by monodromy. **One feasible $t$ projects dominantly onto ratio space**, so fixed-$t$ elimination leaves
  the zero ideal: the coordinator's fixed-$t$ Schwartz-Zippel route is closed.
* **III.3.** Exact Gaussian realisation test for odd private integral norms; the example $t\equiv1$ (real-feasible,
  empty arithmetic closure, by parity mod 2).
* **III.4.** On any rational norm plane, fixed signed $t$ gives $O_\varepsilon((2T)^\varepsilon)$ genuine norm vectors,
  uniformly in the plane height, via Bombieri-Schmidt (Invent. Math. 88 (1987) 69-81; verified via secondary sources
  only). Corrected: the ratio-class count is at most the norm-vector count. Four fixed marks admit at most four pure
  ten-pair completions per cofactor product. Fifteen cuts: the doubleton norm closure is the 7-dimensional cone over
  $\mathrm{Gr}(2,5)$, so no ten-pair sextic survives (negative). Genus-four cube cover (Riemann-Hurwitz $2g-2=6$),
  Faltings finiteness for fixed anchors only.
* Interpolation bound $O_\varepsilon(1+(AK)^{1/4}L^{(1+\varepsilon)/4})$: routine; the growth problem supplies no bound on
  $A,K,L$.

**Plausible only:** the final identification of the "one rational end / Pell" bridge with the strict linear system.

The positive counts concern the pure private ten-pair model, which II.1's Walsh bound already makes $O(m^4)$ on any
growing grid.

### 1.7 linear: lemma III.5, III.7-III.10 (referee: no false claim; two eliminations unreproduced)

**Proved**
* **III.8.** The strict linear-block scheme modulo affine maps has Zariski tangent dimension $\le3$ in characteristic 0
  and every $p>5$, and $\le2$ where $B\operatorname{diag}(x_e)Z$ has rank 4. (Expected dimension $-2$.) New.
* **III.10 (barred list).** No strict barred completion at rank 5, at most 1 at rank 4, at most 9 at rank 3; rank $\le2$
  impossible. The $\mathbb F_{101}$ rank-4 example reproduced (same jet up to normalisation).
* **III.5.** No strict solution with one or two distinct real parts. Pivot certificates: strict scheme of dimension
  $\le4$ over the anchor base (weaker than III.8).
* Real sign graph of upper-half-plane roots is empty, complete or a pentagon.
* III.9 set-up, Hankel test, plane count $O(P^{10}\log P)$, cube cover of genus 4; III.10 ratio polynomials of degree
  $O(6^{24d_0})$. Conditional or peripheral.
* **III.7 relabelling acts freely, corrected:** the case $d=4$ with a swap is skipped and must be added (trivial).
  *Conditional on the reflection slice below.*

**Computed.** Strict twisted census over $\mathbb F_{p^2}$ empty for $p=7,11,19,23,31$ (two independent algorithms,
with controls). New [referee]: at $p=7$ even four moments have no strict twisted point; at $p=11$ they have 220.
*Corrected meaning:* a strict $\mathbb Q(i)$-solution has a collision (possibly a root becoming real) mod each of
these $p$ in every $p$-integral chart. **Nothing follows over $\mathbb C$, $\mathbb R$, or for $\mathbb Q(i)$-existence.**

**Unreproduced (evidence only):** the two-height degree-32 elimination ($p=17,23,103,107$) and the reflection
degree-42 elimination with the golden-ratio exclusion. Singular/msolve were unavailable.

**Ruling on the referee's GAP.** The decision of the strict linear system (solution or emptiness) is open over
$\mathbb Q(i)$, $\mathbb R$ and $\mathbb C$, and Astra says so. This is not a gap in a claimed proof. Either answer is
peripheral: a solution gives one curve (by III.6 at most two-dimensional families), emptiness kills only block degree
$d=1$.

### 1.8 sweep: lemma III.6, III.10 templates (referee: no error; scope corrected)

**Proved**
* **III.6.** For fixed block degree $d$, the union of all strict degree-$d$ ten-pair template ratio curves has Zariski
  closure of dimension $\le2$ in $\mathbb G_m^4$ (char 0). The margin is exactly the $+1$ from the source derivative,
  checked exactly on the $k=4$ family (14 poles vs order 13; projected 12 vs 9). Corollary: leading cofactor shapes lie
  in a closed set of dimension $\le1$.
* **III.10** ten-pair: ratio polynomial of degree $O(6^{24d_0})$, $o(m)$ for $d_0\le c\log m$, $c<1/(24\log6)$.
* **III.10 fifteen cuts, corrected statement.** *Monic disjoint-support templates with contact $\ge\lceil D/2\rceil$ and
  $\max d_A\le d_0=o(\log m)$, evaluated at **real** source values, lie on a ratio polynomial of degree
  $O((d_0+1)^{15}8^{30d_0+4})=m^{o(1)}$, hence make up $o(m^5)$ tuples on any $m$-mark grid, uniformly in the
  coefficients* (already $o(m)$ for $d_0\le c\log m$, $c<1/(30\log8)$). The 4x4 log-determinant inequality sharpens
  `toric` Prop 4.3 by one, and that sharpening is needed at even equal degree.
* Relation to hypothesis (5) of Part I §5: supplied only on the template subclass. Not covered: block degree
  $\gtrsim\log m$, special finite $x_0$, large coefficients, multivariable families, isolated tuples. Hypothesis (5) is
  not advanced.

### 1.9 Withdrawn or corrected in Part III (collected)

1. III.9: "order the sampled primes" means a uniformly random order; the headline 17/13 is superseded by 2/1 (pair
   count); the exceptional set depends on $C$.
2. III.11: "asymptotically exact radii" means relative error $O(N^{-1/10})$; 0/25 marked incidences isotropic.
3. III.2 assertion 2: a code-size bound, not emptiness; "angle fakes fail" is vacuous.
4. (III.4): additive constant $O_B(1)$.
5. III.7 / III.8: explicit thresholds $m\ge16{,}390$ (or $\approx10^4$), $\gtrsim4.1\cdot10^4$, $\gtrsim5\cdot10^5$.
6. Lemma III.7 relabelling: add the $d=4$ swap case; conditional on the unreproduced reflection elimination.
7. Lemma III.4: ratio-class count $\le$ norm-vector count. Lemma III.10: real source values only.
8. Census: bad reduction at 7, 11, 19, 23, 31 only; no geometric conclusion.

## 2. What is genuinely new relative to the record, and correct

* **Recovery and rigidity from separation divisors** (free III.2 recovery lemma, III.7): for $m\gtrsim10^4$ the labelled
  $d_{ij}$ plus the height bound determine the $M_{0,m}(\mathbb Q)$ point. A clean converse-type rigidity. The Plücker
  identity itself is old (`growth/ptolemy.md` Thms A, B).
* **Linear lower bound on unit-equation rank** (free III.4, $s\le8r_0+4$): turns the record's "exponential in the rank,
  vacuous" into a proof that route (a) cannot work for anchored equations. F1-transparent.
* **Random-parent almost-all bound** (free III.9): new to the record (grep: no random-parent statement of this shape;
  Bourgain-Rudnick's almost-all separation concerns density-one $n$). For $D>2$ the phenomenon is elementary
  (pair count); the quartet-exposure argument is genuinely needed only for $1<D\le2$.
* **Integer critical-arc model at five marks** (free III.11): the split-only blocks and the index-$D$ lattice change that
  give critical archimedean geometry with full integrality. Projectively a point of the C1 family
  (`growth/ptolemy.md` Thm C1, `threshold/L8.md` §4).
* **Fixed-cofactor geometry** (lemma III.2-III.4): degree-12 fourfold, dominance onto ratio space, the
  $\mathrm{Gr}(2,5)$ cone at fifteen cuts, a plane-uniform Bombieri-Schmidt count.
* **Strict linear system** (lemma III.8, III.10): tangent dimension $\le3$, barred completions $\le1/\le9$; censuses at
  $p\le31$.
* **Template sweep** (lemma III.6, III.10): closure dimension $\le2$ (was $\le3$), effective in the degree.
* Minor: five-set cross-ratio consistency (no-go), balanced four-cycle example, parity of cofactor graphs, III.3
  no-amplification, $t\equiv1$ parity example.

Not new: the Plücker identity (ptolemy Thms A, B), "every assignment admits an angle fake" (round16 Thm 2.3), the
star/point dictionary, (L_5) false via C1.

## 3. What changes for the integer problem

### 3.1 Growth rate, uniform bound, Questions A and B, route 5

* **Growth rate:** unchanged, $(2/3+\varepsilon)\log R/\log\log R$. The constant 3 in $W\ge3m\log m-O(m)$ is untouched.
* **Uniform bound $M(C)$:** unchanged. III.9 is not a uniform bound for any individual $R$.
* **Question A ($(L_8)$):** unchanged. No Part III result has content at $m=8$: III.7 needs $m\gtrsim10^4$; III.11 is a
  five-mark object and its construction provably stops at $m=6$; III.6 keeps only half the divisibility.
* **Question B ($(L_5^\circ)$, $A_5(\delta)$):** unchanged. III.11 shows that integrality plus all divisibilities
  (without exact norms or isotropy) cannot give $A_5(\delta)$; template-density refutation is now blocked in dimension
  3 as well as 4; the strict linear template remains undecided.
* **Route 5:** unchanged; no Part III item touches the zero estimate.
* **Routes closed (record these):** (a) unit-equation counting with fixed anchors (rank must be $\gtrsim3s/16$);
  fixed-cofactor elimination (dominant); quartet cross-ratio consistency (automatic); pure ten-pair counts
  (already $O(m^4)$-exceptional); random-parent density statements as a route to the worst case.

### 3.2 Ruling (a): the density theorem III.9

1. **Is it proved?** Yes, with the random-ordering repair, for each fixed $D>1$ and $C$, with exceptional fraction
   $\exp(-c\,r\log r)=N^{-c/D+o(1)}$, $c=c(D,M,C)$ ("exponentially small in $r$", only *polynomially* small in $N$).
2. **Which circles and arcs?** Exactly those of §1.1 "Coverage": the circle $N$ for good $N$, and clusters on
   $n=N|G|^2$ sharing a common Gaussian factor $G$, on arcs of length $C\sqrt R$ for the fixed $C$. Not clusters in which
   an additional split prime is a non-constant cut. The parents form a density-zero family of integers with
   $\omega(N)\approx\log N/(D\log\log N)$, i.e. growth-scale numbers of primes, all of one dyadic size.
3. **Is it the record's first statement that genuine circles satisfy and angle fakes do not at the growth scale?**
   **No, in any useful sense.** [J] (i) It is a fixed-$M$ statement, not one at the growth scale in $m$; it is at growth
   scale only in $\omega(N)$. (ii) Any almost-all statement over parents is trivially violated by fakes, because a fake
   chooses its angles adversarially while the theorem averages over the genuine angles; the referee's pair count
   already "separates" in this sense at $M=2$ for $D>4$ (a 2-mark angle model exists on every parent). The input used
   is the count of Gaussian integers in thin sectors with prescribed support, an equidistribution average of exactly
   the kind the record says cannot reach a worst parent (`round10/analytic.md` Thm B, Prop C). (iii) The record's
   fakes are built on narrow bins $[P,Pe^{1/(4r)}]$ (`round16/barrier.md` Thm 2.3). Random dyadic draws land in such a
   bin with probability about $(c/r)^r=e^{-r\log r+O(r)}$, which is far *smaller* than III.9's exceptional fraction
   $e^{-c r\log r}$ with $c\approx0.02$ at $D=4$, $M=18$ [J, from $(D-1)F_{18}-3D/34\approx0.019$]. So the
   exceptional set may contain every parent on which the record builds fakes, and the theorem cannot distinguish
   genuine from fake there. Whether angle fakes with growing $m$ exist on *typical* dyadic parents is not established
   in the record; without that, no separation is shown even on typical parents. (iv) Precedent: Bourgain-Rudnick
   almost-all separation is already a genuine-only statement (density-one $n$). III.9 is new only in applying to a
   growth-scale-$\omega$ family.
4. **Wording for the record** (recommended, so that nobody mistakes it for a worst-case bound):
   > *Random parents (round 31, Part III, almost-all; not a worst-case statement).* Let $N$ be the product of $r$
   > distinct split primes drawn uniformly, in random order, from $[r^D,2r^D]$, and fix $C>0$. Outside a set of parents
   > of relative size $\exp(-c(D,C)\,r\log r)=N^{-c'/D}$, normalised clusters of norm $N$ (and their multiples by a common
   > Gaussian factor) have at most $M_D$ points on any arc of length $C\sqrt R$, where $M_D=1$ for $D>4$, $M_D\le2$ for
   > $8/3<D\le4$, $M_D=O(1)$ for $D>2$ (elementary pair count), and $M_D\le37$ at $D=3/2$, $M_D<\infty$ for every $D>1$
   > (quartet exposure). This is an average over a density-zero family. It says nothing about any individual $N$, about
   > parents with an extra non-constant split prime, about the narrow-bin parents on which fakes are built, or about the
   > growth rate.

### 3.3 Ruling (b): the integer countermodels III.6 and III.11

1. **Five points.** **Yes, III.11 shows that integrality of the star matrix with all divisibilities is insufficient at
   five points** to beat the exponent 10/3, when the remaining data are: primitivity, bounded cofactors, index one,
   the chain condition (vacuous here, pairwise coprime blocks), norms to relative error $O(N^{-1/10})$, and a narrow
   cone. The referee verified that after anchoring it is an $(L_5)$ datum with constant cofactors ($H_\Sigma=O(1)$).
   *Qualifications:* 10/3 is already the proved value at five points; $(L_5)$ was already known false via C1 and
   III.11 is a point of C1 projectively; the model is non-isotropic at every marked prime tested. What it adds is that
   the C1 counterexample survives full integrality and critical archimedean geometry. The missing property is
   **exact individual norms**, and Astra proves that exact norms with all divisibilities force a genuine cluster. So
   the gap between "integral model" and "circle" is exactly the problem, not a new lever.
2. **Growing $m$.** III.6 does **not** show insufficiency of integrality at the growth scale: it keeps only half the log
   divisibility (the trivial CRT ceiling), and with full divisibility the case is open. No integral fake with
   $m\gg W/\log W$ and full divisibility exists in the record.
3. **Question A.** Nothing. $(L_8)$ concerns eight marks; III.11 is fixed at five and pair-cut models stop at six;
   III.6 is half-divisible. Consistent with the record's expectation that $(L_k)$ can hold only for $k\ge8$.

### 3.4 Sentences of the record that should be updated (recommendations; the record is read-only)

1. `expert_note/expert_note.md` §1 after the F1 paragraph: add the §3.2 item 4 wording, with "almost-all; not a
   worst-case statement".
2. Wherever unit-equation counting is listed as a possible route: add "closed for anchored equations: the coset rank
   must be at least $(s-4)/8$ (round 31, III.4); the obstruction is F1-transparent".
3. Where (L_5)/C1 is discussed (expert note §2 item list): add "C1 has full integer realisations at critical
   archimedean scale with all divisibilities and bounded cofactors (round 31, III.11); only exact norms are missing".
4. Where `toric` Prop 4.3 is cited: the effective log-determinant inequality $(k-1)\mu_0+1\le B\cdot C$ and closure
   dimension $\le2$ for ten-pair templates (round 31, lemma III.6).

## 4. Open sub-questions, ranked by honest probability of a worst-case growth-rate improvement

All probabilities are the judge's estimates [J].

1. **A worst-case quartet obstruction.** Show that on *every* large parent (not only narrow-bin Paley parents) some
   regular quartet admits no small integer Plücker solution with the prescribed divisibilities (Astra's III.8 target).
   P ≈ 0.01. *First step:* on the record's golden family members and on the planted exceptional parents, measure the
   fraction of passing quartets as $m$ grows, and try to prove that "all regular quartets pass" forces a narrow-bin
   arithmetic condition on the primes.
2. **The III.10 residual scaling problem.** At a large exceptional profile, show that the unique reconstructed
   configuration admits no rational scales $s_i$ making all $c\,s_is_jh_{ij}$ small nonzero integers, plus primitivity and
   the single quadratic form. P ≈ 0.005. *First step:* run it on the seven-mark example of the planted fourteen-prime
   parent (shared passing three-per-class partitions), as Astra proposes.
3. **Cofactor pigeonhole, route (b).** Test the density of the shared complementary-parent conditions
   $c\mid N/d$, $a\equiv\mp s\pmod h$ across many common neighbours of one base edge at an actual large parent; cycle
   exclusion alone cannot work (the four-cycle example). P ≈ 0.002. *First step:* compute these conditions on the
   golden 8-point clusters and on random parents with planted clusters.
4. **An integral fake at growing $m$ with full divisibility** (or a five-mark model with norms finer than
   $N^{-3/10}$, or isotropic at a positive proportion of primes). This would be a new barrier, so its direct value for
   an improvement is low, but it would locate the obstruction precisely. P ≈ 0.002. *First step:* extend the III.6 CRT
   construction to keep the full $g_{ij}$ on a growing subset of pairs and see where the CRT period exceeds the target
   window.
5. **Decide the strict linear template system.** P < 0.001. *First step:* the eight-variable rank-4 chart of lemma
   III.10 plus the rank-3 pencil in msolve/Singular over two or three primes and over $\mathbb Q$.
6. **Sharper random-parent statements** (e.g. 1 point for every $D>0$ via a second moment, or other families). P ≈ 0.
   *First step:* none recommended; these cannot see the worst parent.

## 5. Scores (0-10)

| block | score | reason |
|---|---|---|
| density (free III.9) | 2/10 | Correct (one-word repair), new to the record, but for $D\ge4$ weaker than an elementary pair count, and almost-all in a density-zero family. |
| quartet (free III.1-III.3) | 2/10 | All correct. Code-size bounds at one special parent; the identity is old. A single admissible labeling would suffice for growth. |
| rank (free III.4-III.5) | 3/10 | Clean linear rank bound that rigorously closes route (a); F1-transparent; the $C_{10}$ calibration is hypothetical. |
| model (free III.6, III.11) | 3/10 | A correct, carefully verified integer realisation of C1 at critical scale; confirms what the record expected; III.6 is near-trivial. |
| proj (free III.7, III.8, III.10) | 3/10 | Correct and new rigidity, but only for $m\gtrsim10^4$; III.8 a no-go; III.10 a reformulation. |
| four (lemma III.1-III.4) | 3/10 | Correct, mostly new; main value negative (fixed-$t$ route closed). |
| linear (lemma III.5, III.7-III.10) | 2/10 | Correct partial bounds and censuses; the decisive question undecided, and either answer is peripheral. |
| sweep (lemma III.6, III.10) | 3/10 | Correct dimension drop and effective degree; covers only low-degree templates. |
| **Part III, on the goal** | **1/10** | Competent and nearly error-free, with useful route closures. No change to the growth rate, the uniform bound, Questions A or B, or route 5. The most visible theorem (III.9) is an almost-all statement whose headline is superseded by an elementary count. |
