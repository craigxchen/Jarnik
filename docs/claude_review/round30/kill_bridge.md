# Round 30, adversarial kill of lane bridge

Scripts: `round30/kill_bridge_checks/` (independent code; none of the lane's scripts is imported). Labels:
**re-proved**, **recomputed (exact)**, **evidence**, **corrected**, **not checked**.

## 0. Verdict

**Survives, with minor corrections. Status as claimed: obstruction plus calibration. No rank-5 cluster, no rank bound,
no record improved.** Every decisive number I recomputed matched. One numerical claim (the swept dimension 3) is now
**proved** by an exact Jacobian computation. The errors are a mislabelled bound at $D=7$, a mis-attributed constant in
the literature paragraph, and a "dense" claim that the summary states without its numerical label. None of them
changes a conclusion. Novelty: the template theorem is a routine specialisation/fixed-divisor argument and Lemma 3.2 is
elementary. The genuinely new content is computational ($C_{\rm cube}=145877.658\ldots$, the code identification of
the 5.754 cluster) plus the correction to (S4). Overall this is a repackaged known barrier. The stop rule applies.

## 1. Theorem 1.2 (templates): re-derived

* **Lemma 1.1, re-proved.** Put $\delta=z_j-z_i$. Then $z_i\bar\delta+\bar z_i\delta+\delta\bar\delta=0$, and the
  $t^{3D/2}$ coefficient gives $\gamma+\bar\gamma=0$. Exact check (`kb_cube.py`): all 27 nonzero $\gamma_{ij}$ of the
  cube at $y_J=(-\frac94,-\frac23,\frac14,\frac83)$ are purely imaginary, and
  $\Gamma=|\gamma_{0,\mathbf 1}|=153335/1296=118.314\ldots$.
* **(1), re-proved.** $C(t_0)=\max|z_i(t_0)-z_j(t_0)|/\sqrt{|g(t_0)|\,|z(t_0)|}\cdot(1+o(1))$. This is the
  best possible for given $t_0$: any admissible divisor of the values divides $g$, and a scaling
  $w\mapsto\lambda w$ multiplies $C$ by $\sqrt{|\lambda|}$.
* **(2), re-proved (split case).** For $\ell\notin P$ the residue argument works as stated. Primitivity is used to
  find a $z_j$ that omits the unique root near $t_0$, so it is needed for each sub-cluster separately. That holds for
  all subsets of size $\ge5$ of the cube (every column has weight 4). Boundedness and attainment of $g_\ell^*$ on a
  ball are correct: near each local root image the minimum is eventually constant.
* **(3) and (4), re-proved.** (3) is CRT. (4) uses the standard fact that the gcd of values of coprime binary forms
  divides their resultant ideal. One gap, harmless: (4) controls only $|t_0|\to\infty$. The case of $t_0$ bounded with
  $q\to\infty$ also gives $C\to\infty$, provided the differences $z_i-z_j$ have no common real zero. That proviso
  holds generically and is not stated.
* **Odd $D$ / $K<D/2$ gives $C\to0$, re-proved.** $|g|$ is bounded below on integers, so $C=O(|t_0|^{K-D/2})$.
* **Novelty.** This is the standard specialisation argument plus the fixed divisor of values of a polynomial family.
  It is correct, and it is a formalisation of (S6) and of `families.md` Lemma P rather than a new idea.

## 2. Cube constants: recomputed independently

All of the following are my own code: Gaussian Euclid, exact values, and local maxima found by $\ell$-adic sampling
rather than the lane's nearest-candidate argument.

| claim | my result | script |
|---|---|---|
| cluster at $y_J$ has exactly 8 points; all 256 transversals; $N$ squarefree; 16 distinct non-real roots | max clique 8, **unique**; $N$ squarefree, real coefficients | `kb_cube.py` |
| pair degrees 4, except $\deg(z_P-z_{P^c})=3$ | confirmed | `kb_cube.py` |
| $P$ (bad primes) | 23 primes, $2,\dots,14389$ | `kb_local.py` |
| $g_2^*=2^{-25/2}$, $g_3^*=3^{-8}$, $g_5^*=25$, all other $g^*_\ell=1$ | identical. For $\ell\ge7$ the test "some residue mod $\ell$ has $v_\ell>0$" is complete and finds none | `kb_local2.py` |
| $C_\infty(y_J)=145877.658383\ldots$ | $145877.65838315$ | `kb_local2.py` |
| same value for the 3 $\varepsilon$-assignments | confirmed (via $y$-permutations) | `kb_class.py` |
| $C(13999/2)=145877.668212\ldots$; $R^2$ has 74 digits, squarefree, 19 primes all $\equiv1\ (4)$ | confirmed exactly; first point matches `example.txt` up to a unit | `kb_gauss.py`, `kb_records.py` |
| $C\to145877.658383$ at the 200- and 488-digit $t_0$ | confirmed (needs $\ge500$-digit precision; at 40 digits the arc computation gives garbage) | `kb_gauss.py` |
| sub-cluster constants $30.88$ / $432.87$ / $7946.5$, attained by the stated subsets | $30.8826$, $432.871$, $7946.46$, same subsets | `kb_sub2.py` |
| family $\{\pm a,\pm1/a\}$ degenerate | all 24 orderings have a coinciding root, for every $a$ | `kb_degen.py` |
| second class $y=(-\frac{25}6,\frac23,\frac{18}5,-\frac1{10})$: $1.89\cdot10^7$ | see §2.1 | `kb_class.py` |

Caveats that stand:
* The 2-adic and 3-adic maxima come from finite sampling (mine: denominators up to $2^5,3^3$, residues mod
  $2^{13},3^8$). They agree with the lane's deeper enumeration. A missed higher maximum would only lower $C_\infty$, so
  $145877.66$ is a proved realised upper bound in any case, attained along the explicit progression.
* "Smallest over the family" is evidence only: 28 of 38 classes were computed, and the rest of the surface was not
  searched. The lane says so in §6 but not in the summary sentence "the smallest $C_\infty$ found".

### 2.1 Second class
My result is $\Gamma=527.928$ and $C_\infty=1.8898\cdot10^7$, matching the lane's $1.89\cdot10^7$ and its factor
$\approx130$ over $y_J$. The set $P$ has 27 primes. I skipped the local tests at the three largest
($280673$, $391921$, $6107657$). A gain at any of them would only lower $C_\infty$, and the agreement suggests there
is none.

## 3. Records (§3)

* **Lemma 3.1, re-proved.** Good $\pi$ divides $t_0-\rho$ for at most one of the $2D$ roots. The column of $W$ at
  $\pi$ is $v\cdot\mathbf 1_{\{j:\rho\in z_j\}}$, which is a function-field column up to complement. $|W|\le M$ is
  trivial.
* **Lemma 3.2, re-proved.** For $p=\pi\bar\pi\,\|\,n$ we have $v_\pi(w)+v_\pi(\bar w)=1$, so every cut is a
  transversal of the 4 pairs. The transversal indicators span a space of dimension $\le5$ containing constants, so the
  rank is $\le4$.
* **The 5.754 cluster, recomputed (exact; `kb_records.py`, own factorisation).** $n=5\cdot13\cdot17\cdot29\cdot53\cdot181\cdot457\cdot733\cdot3389$
  is squarefree, $C=5.754210$, there are 8 distinct cuts, and the set of cuts is **equal to the set of all 8 transversal
  cuts** of the conjugate pairs $(a-d,\pm y_d)$. Rank 4. Confirmed. Every 8-point record listed is an axis cluster
  $(a-d,\pm y)$, so conjugation-closed. Only the squarefree one is forced into the cube code; the others have $5^2$,
  $5^3$ layers, as stated.
* **Golden family "level $D/4$": not re-derived beyond the degree count.** That count is
  $|Z|\asymp X^4$, bounded $C$ needs differences $\asymp X^2$, so $\deg(z_w-z_{w'})=6=D/2+2$ at $D=8$, which is
  consistent. The rest rests on `families.md` §2 (refereed earlier).
* **5.674 record and 10-point cluster "sporadic": evidence,** correctly labelled. The 10-point statement only says that
  a source would need $M(D,0)\ge10$, which is trivially true and correct.

## 4. Question B claims (§4)

* **Recomputed exactly (`kb_q5.py`, symbolic in $y$):**
  - all 56 five-subsets have $g=0$, $\mu_0=4$, $\alpha=2(8-g)/\mu_0=4$;
  - code ranks: 32 of rank 4, 24 of rank 3;
  - cone identities: $0$ of $840$.
  The counts 24/32 also follow by pigeonhole: 5 points from 4 conjugate pairs contain one full pair (32 ways) or two
  (24 ways).
* **$\sigma$-invariance reduction of the cubic threefolds to the cones, re-proved.** $\sigma C=C$ because
  $1/\overline{z_j(t)}=z_j(t)/N(t)$ for real $t$, and $\mathcal Q$ is defined over $\mathbb Q$.
* **Swept dimension: upgraded from numerical to PROVED.** I parametrised the surface by
  $(y_0,y_1,y_2)$, $y_3=-\sum$, constrained by $\prod y=1$. At the exact rational point $(y_J,t=1/3)$ I computed the
  exact Jacobian of $(y,t)\mapsto(z_m/z_{S_0})_{m\in S}$ along the 2 tangent directions and $t$. It has **rank 3 for all
  56 subsets**. Hence each 5-subset family sweeps an irreducible 3-fold. For the 24 rank-3 subsets that 3-fold lies in
  the connected 3-dimensional subtorus $u_au_b=u_cu_d$ (the character is primitive), so it is dense there. Proved.
* **Correction (labelling).** The lane's summary states "24 of the 56 are dense in the subtori" without the numerical
  label it carries in §4.2. With the exact Jacobian above the statement is now proved, so the summary is right after
  the fact.
* **$W$, $\iota W$: evidence only** (not recomputed). The lane labels this correctly.
* **Dictionary $\alpha<4\iff g>2c$ (when $\mu_0=D/2-c$), re-proved.** $-K\cdot C=2D'$ needs each nonconstant column to
  be a primitive Losev–Manin ray met transversally, which holds for squarefree $N/A\bar A$.
* **$\phi(M)$ averaging, re-proved.** A column of weight $w$ is constant on $\binom w5+\binom{M-w}5$ five-subsets, and
  this is minimised at $w=M/2$ by convexity. $\phi>0$ iff $M\ge9$, and $\phi\to1/16$.

## 5. Calibration and the (S4) correction (§5)

**Re-derived; the correction is right.**
* $X_r$ is smooth toric, $e\notin\partial X_r$, and $Y=\mathrm{Bl}_eX_r$, so $K_Y=\pi^*K_X+(r-1)E$ and
  $K_Y+\pi^*\partial+E=rE$. The divisor is snc ($E$ is disjoint from $\pi^*\partial$).
* For the cluster curve, $\deg\tilde f^*E$ is the multiplicity at $e$, namely
  $\min_i\mathrm{ord}_{s}(z_i/z_1-1)=D-K=L+1$. $N^{(1)}_{\mathcal D}$ counts at most $|S|$ boundary points and the
  single point $t=\infty$ on $E$.
* Truncated Vojta with $2g-2=-2$ then gives $r\mu\le|S|-1+\varepsilon h$. Without $E$ one gets the coordinator's
  $(r-1)\mu\le|S|-2+\varepsilon h$. Both are instances, and the one with $E$ is stronger. At $c=0$, $\mu=\lceil D/2\rceil$
  (also for odd $D$, where $\mu=(D+1)/2$). So $r=4$ is allowed once $\varepsilon h\ge1$, and $r=5$ needs
  $\varepsilon h\ge D/2+1$ (even $D$), which fails for small $\varepsilon$.
* The threshold "$\varepsilon<1/4$" assumes $h=2D$. It depends on the ample class $A$ (with
  $A=\pi^*H-\delta E$, $h_A=h_H-\delta\mu$), so only "for $\varepsilon$ below an absolute constant" is
  canonical. This is a cosmetic point.
* $\varepsilon=0$ gives exactly Theorem D, $r(L+1)+1\le2D$: checked.
* The remark that the cube at $D=8$ is Vojta-exceptional for $\varepsilon<1/16$ is true at $D=8$ only. Its
  compositions ($D=8q'$, $h\asymp16q'$) are allowed for $\varepsilon\ge1/(16q')$, so composed cubes are not forced into
  $Z_\varepsilon$ for fixed $\varepsilon$. The statement as written is not wrong, but the reader should not infer that
  "the cube family lies in $Z_\varepsilon$".
* **So (S4)'s "r=5 is borderline, (R$_4$) sharper than Vojta" is refuted:** (R$_4$) is exactly the shape the full
  divisor predicts, and its content is the exclusion of $Z_\varepsilon$.

## 6. Lower bounds and §6

* Prop 6.1(a)–(c) re-proved: they follow from Theorem 1.2. Composition evaluates the same polynomials at
  $h(t_0)$, which still has bounded denominators.
* The implications from unbounded Q(i)(t)-clusters to refutations of $A_5(\delta)$ ($\delta<1/30$) and $V_E(5)$ but
  not of $(L_5^\circ)$ match the refereed J4 statement verbatim (`expert_note.md` lines 392–401). The record's
  "$M(C)\ge5$ for every $C>0$" is confirmed in `growth/families.md` l.58.
* **Correction (§1, "Consequence for odd $D$").** "ff Prop 3.3 at $\theta=4/7$ gives $\le5$ at $D=7$" holds **only for
  codes of real rank $\le4$**. A $D=7$, $c=0$ cluster can a priori have rank up to 6. The proved general bound is
  Plotkin $A_2(7,4)=8=D+1$. So "$M(7,0)\le5$" is not proved, and a 6–8 point $D=7$ cluster is not excluded. This is
  exactly the lane's own lead. I did not search for one.

## 7. Literature

* Corvaja–Zannier, "Some cases of Vojta's conjecture on integral points over function fields", J. Algebraic Geom. 17
  (2008), arXiv math/0512074. **Confirmed** by WebSearch (arXiv listing, Udine repository), together with an addendum
  in Asian J. Math. (2010).
* Corvaja–Zannier, "Greatest common divisors of $u-1,v-1$ in positive characteristic…", JEMS 15 (2013) 1927–1942.
  **Confirmed** (EMS Press / EuDML / mathdoc listings).
* **Correction.** Search snippets give the bound $3\sqrt[3]2\,(\deg u\deg v\,\chi)^{1/3}$ as the **characteristic-0**
  statement recalled in the JEMS paper. The positive-characteristic version takes a maximum with a term in
  $\deg u\deg v/p$. The lane calls $3\sqrt[3]2$ "the stated positive-characteristic constant". That is mis-attributed.
  The vacuity conclusion stands: $3\sqrt[3]2\,(2D^3)^{1/3}\approx4.76D>D/2$. Both checks are from search results only;
  no full text was read.
* Brownawell–Masser, Voloch and Wang are not cited by the lane beyond the coordinator's list; I did not check them.

## 8. Corrected statement

> **Lane bridge, round 30 (corrected).**
> * **(Templates; proved.)** A primitive Q(i)(t)-cluster of size $M$ with $D$ even at $c=0$ gives integer clusters
>   with $C(t_0)\to\Gamma/\sqrt{|g(t_0)|}$ along $t_0\to\infty$ with bounded denominators. The optimum
>   $C_\infty=\Gamma/\sqrt{\prod_{\ell\in P}g^*_\ell}$ is attained along a progression, and $M_\infty(C)\ge M$ for
>   $C>C_\infty$. $K<D/2$ (e.g. odd $D$) gives $C\to0$. This is a routine but correct formalisation of (S6).
> * **(Cube; exact.)** At $y=(-\frac94,-\frac23,\frac14,\frac83)$, $C_\infty=145877.658\ldots$, independently
>   recomputed. The 74-digit example is exact. Sub-clusters give $30.88$ (5 points), $432.87$ (6), $7946.5$ (7).
>   Minimality over the family is evidence (28 of 38 low-height classes). No record is improved: the golden family
>   gives $5.8977$.
> * **(Records; exact.)** The squarefree 5.754 cluster carries exactly the cube code (Lemma 3.2). The other records
>   are not specialisations of any known $c=0$ cluster (evidence).
> * **(Question B; exact.)** All 56 cube 5-subsets have $\alpha=4$, so they are not A-points. They lie on no quadric
>   cone and no cubic threefold. They sweep 3-folds (**proved** by an exact Jacobian), and the 24 subsets with two
>   conjugate pairs are dense in subtori $u_au_b=u_cu_d$. $W$, $\iota W$: numerical only.
> * **(Calibration; conditional on truncated function-field Vojta.)** With $\mathcal D=\pi^*\partial+E$:
>   $r\mu\le|S|-1+\varepsilon h$. This forbids rank $\ge5$ at $c=0$ off $Z_\varepsilon$ and allows rank 4, so (R$_4$)
>   has exactly the Vojta shape and (S4) is corrected. The Corvaja–Zannier bound (char-0 constant $3\sqrt[3]2$) is
>   vacuous at this scale.
> * At $D=7$, $c=0$ the proved bound is $M\le8$; $\le5$ holds only for rank $\le4$.
> * No rank-5 cluster and no rank bound. Under the pre-registered stop rule this lane gives no reason to keep the
>   function-field line open.
