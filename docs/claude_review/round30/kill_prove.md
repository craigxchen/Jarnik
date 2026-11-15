# Round 30: adversarial review of lane `prove`

Scripts: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round30/kill_prove_checks/`. All of them are
my own code; I did not run the lane's scripts. Notes: `../kill_prove_notes.md`.

## 0. Verdict

**The lane survives, with one substantive correction and one upgrade.** Its theorems check out:
* Lemmas 1.1–1.3;
* Props 2.1–2.3 and 3.3;
* Theorem 3.1 (cube) and Theorem 3.2 (horizontal line).

Where they rest on computation, I recomputed them independently.

1. **The main correction concerns Prop 2.5 and the claim that "balance does not separate".** The claim is false for the
   saturated-lattice form of balance, which is the lane's own Lemma 1.2.
   * The Theorem-B fakes contain an integer kernel vector of $\ell^1$-weight $D/4<L+1$, and so do their doubles.
     Both therefore violate Lemma 1.2(1). The same holds for (3).
   * Only *pairwise* balance (the coordinator's (S2) for difference vectors) is satisfied by the doubled fakes.
   * The statement "the real structure gives no inequality that kills fakes" therefore has to be withdrawn. It is
     open, not proved negative.
2. **Upgrade from evidence to proof:** there is **no 6-point (rank-5) equidistant cluster at $D=8$, $c=0$.** Every one
   of the 3 classes of 6 pairwise-orthogonal codes has a weight-2 integer vector in the span of its differences. Such
   a vector forces $\rho_k=\rho_l$ or $\rho_k=\bar\rho_l$.
3. **Upgrade:** the 5-point rank-4 equidistant clusters at $D=8$ now have a Krawczyk certificate. One correction: the
   full fibre profile at the certified member is $(6,5,5,4,\dots)$. The 6-point fibre has rank 3. The lane reported
   $(5,5,4,4,\dots)$.
4. Smaller corrections:
   * Prop 2.5 "for any code" needs the hypothesis that no unit vector lies in $\mathrm{span}(C-C)$.
   * The (S4) phrase "forbids $r\ge6$ for $\varepsilon<1/2$" depends on the normalisation of $h$.
   * I could confirm the LST citation only by its arXiv number (my search did not find it).

There is still no rank bound in general and no rank-5 cluster. On the stop rule I agree with the lane's literal
reading: the two locus theorems are genuine. But the recommendation to close the line has to be reconsidered,
because the "balance is faked" argument behind it is wrong (§5).

## 1. (S1)–(S3), Lemmas 1.1–1.3

I re-derived all three. They are correct.

* **Lemma 1.1.** On $S$, $z_j$ carries $\bar B$, so $z_i-z_j=A(B-\bar B)$.
* **Lemma 1.2.** $U_a/\bar U_a=\exp(-2i\sum_m t^{-m}m^{-1}\sum_k a_k\mathrm{Im}\rho_k^m)=1+O(t^{-L-1})$.
  * $U_a$ is not real. With $y_k>0$, its upper roots come from $a_k>0$ and its lower roots from $a_l<0$, and the
    $2D$ numbers are distinct.
  * The IVT count on $\Psi_a$ is right.
  * Lemma 1.2 coincides with the search lane's Lemma 1.1.
* **Lemma 1.3.** $dt/(t-\rho)-dt/(t-\bar\rho)=2iy\,dt/|t-\rho|^2$, and the residue of $W_a/N$ at $\rho_k$ is
  $-ia_k$. Both are right.
* **(S3) correction.** The coordinator's $N_v$ version would need the *real* power sums $\mathrm{Re}\rho^m$ to vanish,
  since $\log N_v$ involves $\rho^m+\bar\rho^m$. So it is false in general, as the lane says. The Blaschke version is
  correct.

## 2. Real structure

* **Prop 2.1: correct.** The roots of $P=-i\kappa$ as $\kappa\to+\infty$ are near $(-i)^{1/w}\mu_w$, which gives
  $\lfloor w/2\rfloor$ upper roots. The count is constant on each half-plane of values.
* **Prop 2.2: correct.** It is the identity $\omega(v_2,v_3)v_1+\omega(v_3,v_1)v_2+\omega(v_1,v_2)v_3=0$ in $\mathbb R[t]^2$.
  * $\sum\kappa=0$ follows by comparing top coefficients of monic $m_l$.
  * $\kappa_{ij}=c_i-c_j$ is just the cocycle condition and carries no content.
* **Prop 2.3: correct.** For $w=D/2$ even, $|p-n|\le1$ forces $p=n$, so $\theta_{ij}\to0$ at both ends, and
  $\sin\theta_{ij}$ never vanishes.

**Prop 2.5: false as used** (`fake_lattice.py`, exact).

* *Setup.* Take the Theorem-B fake at $D=8$, $c=0$: $B=\{e_1,e_2,e_3\}$ and $A=\{110,101,011,111\}$, which gives
  10 points in one fibre.
* *The short vector.* $\mathbf 1_H$ with $H=\{000,111\}$ satisfies $W\mathbf 1_H=0$. It lies in the span of the
  differences, $\mathbf 1_H=\frac14(\mathbf 1+\chi_{110}+\chi_{101}+\chi_{011})$, and has weight $2<L+1=4$.
* *The double.* The doubled fake at $D=16$, $c=0$ ($L=7$) is built exactly as in the lane, with rows $(x,x+e)$:
  * it is numerically one fibre (spread $10^{-15}$);
  * its columns are distinct and non-opposite;
  * its first row is positive;
  * it contains the integer kernel vector $(\mathbf 1_H,-\mathbf 1_H)$ with $\|a\|_1=4<8=L+1$ and
    $a^+=2<L/2=3.5$.

  So the double violates both (1) and (3) of Lemma 1.2.
* *Where the proof goes wrong.* The step "$\sum a=0$ so balance holds" forgets that
  $|\sum a|\le\|a\|_1-L$ also requires $\|a\|_1\ge L$.
* *It is not a one-off.* When $B$ spans $\mathbb F_2^k$, the set $A$ has $\ge2^{k-1}$ elements and is not a
  hyperplane complement. So it contains a triple $a,b,a+b$, and the same weight-$D/4$ vector exists.
  (This uses the classification of maximal sum-free sets in $\mathbb F_2^k$, which I quote from memory.)
* *"Any code" also needs a hypothesis.* $C_2=\{1111,--11,-1-1,1--1\}$ at $D_0=4$ has $e_1\in\mathrm{span}(C_2-C_2)$.
  Its annihilator is $(0,0,0,1)$, so columns $k$ and $D_0+k$ coincide for $k=0,1,2$.

**What survives:** *pairwise* balance plus distance is satisfied by doubled fakes with $D/2+2$ points. Whether the
saturated filter bounds $M$ is **open**. Doubling the cube code passes it, but has only $D/2$ points.

## 3a. Cube locus, Theorem 3.1: confirmed

* **Special members** (`cube_special_mine.py`). For the 35 classes of $d\in\{-2..2\}^4\setminus\mathbb Z\mathbf 1$
  (up to $S_4$, sign and shift), I eliminated by resultants on the plane $\{\sum y=\sum dy=0\}$.
  * Real roots were isolated rigorously with flint `complex_roots`.
  * The degenerate components are $d=(a,a,b,b)$, which gives $y_i+y_j=0$, and $d\sim(-1,0,0,1)$, which gives
    $y_0=y_3$. Both force equal $|y|$.
  * The nondegenerate real solutions form one $|y|$-set, $\{0.123790,0.911178,2.609743,3.397130\}$, and its orbit
    under $S_4\times\{\pm1\}$ has **48** points. This agrees with the lane.
* **Fibres** (`cube_fibres_mine.py`), computed directly from the 8 roots at 50 digits.
  * At 30 random generic members the profile is $(8,4,4,4,\dots)$, the max fibre has rank 4, the within-fibre error
    is $10^{-46}$, and the gap between fibres is $\ge0.51$.
  * At all 48 special members the profile is $(8,4,\dots)$, with gap $\ge4$ against error $10^{-42}$.
  * `cube_comb.py` shows that the combinatorial fibres (coincidences only in $\mathbb Z\mathbf 1$) are
    $8,4,4,\dots$, with a unique 8-fibre of rank 4.
* **Scope.** The theorem holds for the symmetric surface $\{\pm r_k\}$. Other components of the cube-type incidence
  variety, if there are any, are not covered.

## 3b. Horizontal line, Theorem 3.2: proof correct

* The rows triangulate to $x^0,\dots,x^{L-1}$.
* Equal $e_1,\dots,e_{L-1}$ for disjoint sets of size $<L$ forces them to be empty.
* The Bonferroni bound gives $3L\le D$.

No counterexample is possible. The locus is very degenerate, so as a "rank bound on a positive-dimensional locus" it
is weak.

## 3c. Equidistant clusters

* **Prop 3.3: correct.**
* **No 6-point equidistant cluster at $D=8$, $c=0$: PROVED** (`equi_lattice.py`, `equi_witness.py`).
  * *Enumeration.* There are 20160 five-cliques of vectors orthogonal to $\mathbf 1$ (frame $s_1=\mathbf 1$, reached
    by column flips, which preserve $\ell^1$). Up to $S_8\times S_5$ they form 3 classes. This agrees with the lane.
  * *Witnesses.* Each class has a weight-2 integer vector in $\mathrm{span}_{\mathbb Q}$ of the differences. In
    class 0, $e_2-e_3=\frac12(d_2+d_3-d_4-d_5)$. All three witnesses were checked exactly.
  * *Conclusion.* By Lemma 1.2(1) none of the 3 classes is realisable.
* **Five-point classes.** Of the 4 classes, two are killed the same way and two survive (6 sign normalisations each).
  This explains the lane's numerics: only those two classes were realised.
* **5-point rank-4 clusters exist: CERTIFIED** (`equi5_cert.py`, `krawczyk5.py`).
  * *Square system.* I gauge-fixed $\rho_0=i$, $\rho_1=0.752+0.999i$ (rational), which leaves 12 equations in 12
    unknowns.
  * *Certificate.* After Newton the residual is $10^{-61}$. The Krawczyk operator maps the box of radius $10^{-30}$
    into its interior, so the box contains a unique zero.
  * *Nondegeneracy, on the box.* $|\mathrm{Im}\rho|\ge0.281$, and the $\ell^1$ separation of the 16 roots is
    $\ge0.32$. The roots are listed in the script output.
* **Fibre profile** (`equi5_fibres.py`).
  * It is $(6,5,5,4^8,\dots)$, with gap between distinct fibres $\ge0.227$.
  * The 6-fibre (rank 3) and the second 5-fibre have differences inside the real span of the certified cluster's
    differences. Their coincidences are therefore exact.
  * The lane's "$(5,5,4,4,\dots)$" is wrong for this member.
* **Not done:** $D=12$ equidistant 6-codes, where the filter could decide the lane's lead.

## 4. Tools, (S4), literature

* **§4.1–4.3: heuristics, labelled correctly.**
  * The Wronskian blindness argument is right as a dimension argument.
  * In §4.3, "the cube saturates $d=3$" is unverified.
* **(S4).** $K+\partial=(r-1)E$ and $\kappa=0$ are right, and so is $(r-1)m\le2D-2+\varepsilon h$.
  * $r=6$ is excluded only when $\varepsilon h<D/2+2$, so "$\varepsilon<1/2$" assumes $h\le D$.
  * The truncated inequality is a conjecture shape.
* **Literature (WebSearch).**
  * Lehmann–Riedl–Tanimoto, *Non-free curves on Fano varieties*, arXiv 2304.05438, Osaka J. Math. 62 (2025) 19–50:
    **confirmed**. The subject is Mor(B,X) for Fano X, with boundedness for non-dominant families.
  * Lehmann–Tanimoto, arXiv 1702.08508, is about smooth Fano X over $\mathbb C$: **confirmed**.
  * Lehmann–Sengupta–Tanimoto, arXiv 1805.10580: **not found by my search**.

  The inapplicability argument stands: $\mathrm{Bl}_eX$ is not Fano, there is a contact condition, and rigid curves
  lie in closed sets anyway. The gcd-theorem shapes remain memory-quotes.

## 5. Stop rule, novelty, corrected statement

* **Stop rule.**
  * There is no rank-5 cluster for $D\le12$ from this lane. It is now proved absent among equidistant clusters at
    $D=8$.
  * There are rank bounds on two special positive-dimensional loci. The literal reading of the rule, "does not
    fire", is correct.
* **Recommendation.** The case for closing the line rested partly on "balance is faked". That is false for the
  saturated form, which kills Theorem-B fakes and decides $D=8$ equidistant. The one concrete, cheap continuation is
  to run the saturated-lattice filter on all candidate rank-5 codes at $D=8$–$12$ (equidistant first, at $D=12$). If
  that filter kills everything, it is a combinatorial route to (R$_4$)-type statements. If it does not, close the
  line.
* **Novelty.** Lemma 1.2 (lattice form) is shared with the search lane. Props 2.1–2.3 are elementary but new in this
  program. Theorem 3.1 (with the search lane) and Theorem 3.2 are new and modest.
