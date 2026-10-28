# Round 20 — wildcard lens

Status: complete author's write-up, **unrefereed**. Labels: proved / exact / evidence / heuristic / quoted (literature, not re-verified).

## 0. Summary

**Lens chosen:** unlikely intersections in tori (Bombieri–Masser–Zannier bounded height on curves, Habegger's
bounded height theorem, Rémond's cones, Zilber–Pink in anomalous form). It is absent from the record, and it is the
only mature unconditional theory whose constants are uniform in the finitely generated group, so it is immune to
checklist items 4–5 by design and is not a local (fake-satisfied) input.

**Outcome: reduced to a cap. No growth improvement, no uniform bound.** Proved here (unrefereed):

* **Fixed carriers are vacuous** (Lemma 6.1, Theorem 6.2). The equal-norm cone $X_M$, the chord space $L_M$
  ($M\ge4$) and $M_{0,k}$ in its Plücker torus ($k\ge5$) all have empty open-anomalous locus $V^{oa}$: forgetting to a
  sub-configuration is a torus projection onto a non-full-dimensional image, so its fibres are anomalous. Habegger's
  theorem and the anomalous form of Zilber–Pink say nothing on them.
* **The moving carrier encodes the arc in its height.** For $u\in\mathbb P^n(\mathbb Q(i))$ let
  $\ell_u=\mathrm{span}_{\mathbb Q}(\mathrm{Re}\,u,\mathrm{Im}\,u)$ (Plücker coordinates $\mathrm{Im}(\bar u_au_b)$). It is
  never in a coset when coordinates are pairwise non-$\mathbb R$-proportional (Lemma 6.3). A short-arc cluster has
  $h(z)=\log R$, $h(\ell_z)\le\frac32\log R+\log C$, chart sum heights $\ge(M-1)(\frac12\log R-\log C)$ (Prop 6.4):
  ratio $\frac23$, the exact Dirichlet exponent (full circles: exactly $\frac12$; Remark 6.5).
* **Uniform-in-the-line comparisons are false at the needed constants** (Prop 6.6). Rank-one points
  $F_1=(1:\gamma:\dots:\gamma^{n-1}:\gamma^K)$, $F_2=(1:\gamma^{K-n+1}:\dots:\gamma^K)$, $\gamma=A+i$, lie in subgroups of
  codimension $n-1$ and have ratios $\to1$ and $\to\frac n2$, above the cluster's $\frac23$ and $\frac n3$; every
  nonnegative mixture of the two heights is beaten too. $(1:\gamma:\gamma^2)$ has $h(u)=h(\ell_u)=\log(A^2+1)$
  exactly. These are the lens's own "algebraic fakes": what they lack is only the equal-norm constraint, and with that
  constraint the needed inequality is Vojta / $(L_k)$ on $X_M$ or $M_{0,k}$.
* **The ceiling of the lens is $M\le\omega(N')+O(1)$** (Lemma 6.7, Prop 6.8, Cor 6.9). A cluster lies in a proper
  algebraic subgroup iff its exponent rows have rank $\le M-2$; full-rank clusters (like every fake profile) are
  invisible. So even an ideal unlikely-intersection input gives only "rank-deficient clusters are bounded", the same
  ceiling as the Subspace route QS (lit18 Prop B), which is not a growth improvement.

**Evidence** (§7): on actual clusters the Plücker gcd is 1 and $h(z)/h(\ell_z)$ equals the bound of Prop 6.4 to four
digits; the referee's 10-point cluster on $1176852625$ is rank-deficient (rank 5), a 7-point window on $3159797225$ is
full rank.

**Proposed barrier-map row:** "unlikely intersections (BMZ / Habegger / Rémond / Zilber–Pink)": fixed carriers have
$V^{oa}=\varnothing$; moving-line comparisons refuted by rank-one algebraic points; ceiling $M\le\omega(N')+O(1)$.
Residual open statement (S$_{\min}$) (§6.3, §8) is non-standard and not recommended.

## 1. Setup recalled (what the record says)

Notation is that of `overview/overview.md` §1.2. $N=R^2$, $W=\log N=2\log R$, and a *normalised cluster* is a set
of $M$ distinct $z_x\in\mathbb Z[i]$ with $|z_x|^2=N$, Gaussian gcd $1$, on an arc of length $\le C\sqrt R$. Then
$z_x=\varepsilon_x\prod_j\pi_j^{a_{xj}}\bar\pi_j^{e_j-a_{xj}}$, $g_{xy}=\gcd(z_x,z_y)$,
$d_{xy}=\log(N/\mathrm{Norm}\,g_{xy})$, and the pair margin bound $d_{xy}\ge W/2-2\log C$ holds for every pair
(from $|z_x-z_y|^2=\mathrm{Norm}(g_{xy})\,\mathrm{Norm}(c_{xy})\ge\mathrm{Norm}(g_{xy})$ and $|z_x-z_y|\le C\sqrt R$).

What the record says about the mechanism class used here (Zilber–Pink / Bombieri–Masser–Zannier "unlikely
intersections"): a search of `docs/claude_review/` and of the research branch finds no use of it. "Zilber",
"Chabauty", "Skolem", "Frey", "Dobrowolski" and "Hilbert symbol" do not occur; "Bombieri–Masser–Zannier" occurs
twice, both times only to say that their anomalous-subvariety theorems do not exclude certain translates in the
zero estimate ZE of route 5 (`round10/README.md`, `round10/referee_determinant.md`); Dimitrov–Gao–Habegger
occurs once in `lit18/tools.md` ("not applicable (genus 0)"). So the lens is untried in the record.

## 2. The proposed mechanism

### 2.1 Why this lens

The judge's summary (`overview.md` §7.6, end) says that a new mechanism must be *global, multi-place and tied to the
actual arithmetic of $N$*, and items 4–5 of the checklist say that every Diophantine tool in the record either
sees one linear form or has constants and exceptional sets that depend on the prime set $S$. There is exactly one
mature body of unconditional Diophantine theorems whose constants are **uniform in the finitely generated group**,
i.e. do not depend on $S$, on $|S|$ or on the rank at all: the theory of **unlikely intersections in tori**.

* **Bombieri–Masser–Zannier (1999).** If $\mathcal C\subset\mathbb G_m^n$ is an irreducible curve over
  $\overline{\mathbb Q}$ not contained in a translate of a proper subtorus, then the points of $\mathcal C$ lying in
  the union of all algebraic subgroups of codimension $\ge1$ have bounded height, the bound depending only on
  $\mathcal C$. (Points in subgroups of codimension $\ge2$ form a finite set.)
* **Habegger (2009), bounded height conjecture.** For $V\subset\mathbb G_m^n$, the points of the open anomalous
  complement $V^{oa}$ lying in algebraic subgroups of codimension $\ge\dim V$ have bounded height, the bound
  depending only on $V$. Here $Y\subseteq V$ (irreducible, $\dim Y\ge1$) is *anomalous* if $Y\subseteq gH$ for a
  coset of a subtorus $H$ with $\dim Y>\max(0,\dim V+\dim H-n)$, and $V^{oa}=V\setminus\bigcup_{Y\ \rm anomalous}Y$.
* **Rémond's cones.** The same with algebraic subgroups $H$ replaced by cones $\mathcal C_\varepsilon(H)=
  \{\gamma y:\gamma\in H,\ h(y)\le\varepsilon(1+h(\gamma))\}$ for small $\varepsilon$ depending on $V$.

(These statements are quoted from memory of the literature, not re-verified in this round; nothing proved below
depends on them. They are used only to formulate what the lens would need.)

The appeal is precise. A cluster's coordinates live in the group generated by $i,\pi_j,\bar\pi_j$, whose rank
$\omega(N')$ is unbounded; every $S$-dependent tool dies on this (item 5, `lit18` Props A, B; overview §5.4 row 14). A BMZ-type
theorem would not see $S$ at all. And it is not a local input: it concerns algebraic points, which fake circles
are not (§5).

### 2.2 How the lens could apply: carriers

To use the lens one needs a **carrier**: a variety $V$ and a map sending each cluster to a point of $V$ that lies in
(or near, in Rémond's sense) an algebraic subgroup of large codimension, together with a height comparison that the
cluster violates. There are two kinds of carrier.

* **Fixed carriers** (one $V$ for all $N$): the equal-norm variety $X_M\cap\mathbb G_m^{2M}$, the chord linear space
  $L_M=\{(z_a-z_b)_{a<b}\}$, the moduli space $M_{0,k}$ of the Ptolemy data in its Plücker/cross-ratio torus.
  Habegger's theorem would apply with constants depending only on $V$, hence uniform in $N$.
* **Moving carriers** (one curve per cluster): the natural one is the **real line of the cluster**. For a projective
  point $u=(u_0:\dots:u_n)\in\mathbb P^n(\mathbb Q(i))$ with Gaussian-integral coordinates write
  $u_x=X_x+iY_x$, $X,Y\in\mathbb Z^{n+1}$. When $X,Y$ are $\mathbb Q$-independent, put
  $$\ell_u=\mathrm{span}_{\mathbb Q}(X,Y),$$
  a $\mathbb Q$-line in $\mathbb P^n$ containing $u$ and $\bar u$. It depends only on the projective point
  ($\lambda u$, $\lambda\in\mathbb Q(i)^*$, has real and imaginary parts in $\mathrm{span}_{\mathbb Q}(X,Y)$). Its
  Plücker coordinates are $p_{ab}=X_aY_b-X_bY_a=\mathrm{Im}(\bar u_au_b)$, and
  $$h(\ell_u)=\log\max_{a<b}|p_{ab}|-\log\gcd_{a<b}(p_{ab}).$$
  The arc condition is encoded in the height of this line: for a short-arc cluster all
  $\mathrm{Im}(\bar z_az_b)$ are $\le CR^{3/2}$, so $h(\ell_z)\le\frac32\log R+\log C$, while $h(z)=\log R$
  (Proposition 6.4). For a generic point $h(\ell_u)\approx2h(u)$. BMZ applies to $\ell_u$ (Lemma 6.3), but with a
  constant depending on $\ell_u$, so a **uniform-in-the-line** height comparison is needed.

## 3. Exact statement needed

Heights: $h(u)=\log\max_x|u_x|$ for Gaussian-integral coordinates with Gaussian gcd $1$ (the absolute Weil height
on $\mathbb P^n(\mathbb Q(i))$); for $\alpha\in\mathbb Q(i)^*$, $h(\alpha)$ is the absolute Weil height; the
**anchor-$a$ sum height** is $h^{(a)}_\Sigma(u)=\sum_{x\ne a}h(u_x/u_a)$, the height of the point
$(u_x/u_a)_{x\ne a}\in\mathbb G_m^n$ in the chart of anchor $a$.

**The class $\mathcal U_n(\varepsilon)$.** Points $u\in\mathbb P^n(\mathbb Q(i))$ with all $u_x\ne0$, with
$u_a/u_b\notin\mathbb R$ for $a\ne b$ (pairwise non-$\mathbb R$-proportional), lying in the cone
$\mathcal C_\varepsilon(H)$ of some algebraic subgroup $H\subset\mathbb G_m^n$ of codimension $\ge1$ (for
$\varepsilon=0$: the coordinates $u_x/u_0$ are multiplicatively dependent).

> **(UI-V) (fixed carrier).** For some $M_0$ there is a variety $V$ (independent of $N$), a map
> $z\mapsto\Phi(z)\in V(\overline{\mathbb Q})$ defined on short-arc $M_0$-clusters, with $\Phi(z)\in V^{oa}$ and
> $\Phi(z)$ in an algebraic subgroup of codimension $\ge\dim V$ (or in a cone around one), and with
> $h(\Phi(z))\to\infty$ as $R\to\infty$. Then Habegger's theorem bounds $h(\Phi(z))$, hence $R$, for every
> $M_0$-cluster landing there.

> **(UI-$\ell$) (moving carrier, projective height).** There are $\kappa<\frac23$, $\varepsilon\ge0$ and a function
> $c(n)=o(n\log n)$ such that for every $n\ge1$ every $u\in\mathcal U_n(\varepsilon)$ satisfies
> $h(u)\le\kappa\,h(\ell_u)+c(n)$.

> **(UI-$\ell_\Sigma$) (moving carrier, chart height).** The same with $h^{(0)}_\Sigma(u)\le\kappa_\Sigma(n)\,h(\ell_u)+c(n)$
> and $\kappa_\Sigma(n)\le(\frac13-\delta)n$ for some fixed $\delta>0$.

**What they would give** (Proposition 6.8): (UI-$\ell$) or (UI-$\ell_\Sigma$) bounds the size of every
**rank-deficient** cluster (exponent rows of affine rank $\le M-2$) by a constant $M_1(C)$, hence every cluster
satisfies $M\le\omega(N')+M_1(C)$. (UI-V) gives the same with "rank-deficient" replaced by "landing in $V^{oa}$".
The thresholds $\frac23$ and $\frac n3$ are exactly the cluster's own ratios (Proposition 6.4), and they sit where
the record's hardness fact puts them: a ratio $h(u)/h(\ell_u)=1/(1+\alpha)$ corresponds to arcs of length
$R^\alpha$, so $\kappa<\frac23$ is the exact Dirichlet exponent $\alpha=\frac12$ (Remark 6.5).

**Verdict on the statements (proved in §6).** (UI-V) is vacuous for every natural fixed carrier (Theorem 6.2:
$V^{oa}=\varnothing$). (UI-$\ell$) is **false** for every $\kappa<1$ and (UI-$\ell_\Sigma$) is false for every
$\kappa_\Sigma(n)<\frac n2(1-o(1))$, already on points of rank one (Proposition 6.6). And even in an ideal form the
lens cannot see clusters of full rank $M-1$, so its best possible output is $M\le\omega(N')+O_C(1)$, which is the
same ceiling as the Subspace route QS (`lit18` Prop B) and is **not** a growth improvement (Proposition 6.8,
Corollary 6.9).

## 4. Barrier checklist, item by item

The lens is checked twice: *a priori* (is it excluded by an existing item?) and *after §6* (where it stops).

| item | a priori | after §6 |
|---|---|---|
| (1) no fake-satisfied input (pairs, gcds, cut, residues, characters, GAP-Fourier, energy, $L^p$) | **passes.** UI hypotheses are statements about algebraic points lying in algebraic subgroups; an angle model is not an algebraic point and has no line $\ell_z$ | its best conclusion, $M\le\omega(N')+O(1)$, is satisfied by every fake profile (Cor 6.9); and the lens has its own **algebraic fakes**, the rank-one points $F_1,F_2$, which satisfy every hypothesis it can state and beat the cluster's anomaly (Prop 6.6) |
| (2) no robust/averaged analytic statement (round 10 Thm B, Prop C, Thm D) | **passes**: pointwise statements about single points | — |
| (3) no single-place auxiliary polynomial on $X_M$ ($\tau^*<2$), none on $\overline M_{0,k}$ | **passes**: bounded-height theorems are height comparisons on a fixed variety, not Liouville at $\infty$ | on $X_M$ and $M_{0,k}$ they are **vacuous**: $V^{oa}=\varnothing$ (Thm 6.2). This is a new cap of kind V, independent of item 3 |
| (4) no single linear form in logarithms; SLW false; SLW$_w$ Vojta-strength | **passes**: no linear form appears | — |
| (5) no $S$-integer argument with threshold depending on primes only through $\lvert S\rvert$ | **passes by design**: UI constants do not depend on $S$, on $\lvert S\rvert$ or on the rank | lands at the **same ceiling** as lit18 Prop B, $M\le\omega(N')+O(1)$ (Cor 6.9): uniformity in $S$ is bought by needing a multiplicative relation among the points, which is guaranteed when $M\ge\omega(N')+2$ and absent for every full-rank cluster (Lemma 6.7) |
| (6) route 5's $(R_H)$ and CIRCLE$(q)$ | not used | — |
| (7) no $(L_5)/(L_6)$ family mechanism at $k=8$; $(L_8)$ not computationally testable | not used | the equal-norm form of the needed comparison lives on $X_M$ or $M_{0,k}$, i.e. in the record's Vojta / $(L_k)$ territory (§6.3; a reading, not a theorem), which the lens cannot reach (Thm 6.2) |
| (8) abc unit-pencil lane closed | not used | — |
| Dirichlet exponent (`round4/hardness.md`) | the needed constant $\kappa<\frac23$ is exactly the Dirichlet exponent $\alpha=\frac12$ (Remark 6.5) | no UI theorem has sharp constants; and without the equal-norm constraint the constant is false (Prop 6.6) |

**Project constraints.**
* *Uniform exceptional sets from fixed-$S$ Subspace:* none used.
* *Gaussian gcd or denominator costs:* every height is computed with the gcds in place (Gaussian gcd of the
  coordinates for $h(u)$; the gcd of the Plücker coordinates for $h(\ell_u)$; the reduced fraction $\beta/\delta$ for
  $h(z_x/z_a)$). The only inequality using a gcd is $\gcd\ge1$ in an **upper** bound for $h(\ell)$, which is the safe
  direction.
* *Drifting phase:* $\ell_z$ is built from the actual points; it is invariant under $z\mapsto\lambda z$,
  $\lambda\in\mathbb Q(i)^*$, and no phase is set to zero.
* *Template-dependent thresholds:* the BMZ constant for the curve $\ell_z$ depends on $\ell_z$, i.e. on the cluster.
  It is **not** treated as uniform; the uniform version is stated separately as (UI-$\ell$) and is shown false
  (Prop 6.6).

## 5. Fake circles, explicitly

* **The lens is fake-violated in its hypotheses.** A fake circle (Definition 5.2 of the overview) has genuine primes,
  a genuine profile and real angles $\phi_j$ that are not arguments of Gaussian integers. Its "points" are not in
  $\mathbb Q(i)$, so neither "the point lies in an algebraic subgroup" nor the line $\ell_z$ is defined. In this
  sense the lens passes the fake test, like $(L_8)$.
* **But its conclusion is fake-satisfied.** By Lemma 6.7 the lens can only act on rank-deficient clusters
  ($\mathrm{rank}\,B\le M-2$). Every fake of the record has full rank $M-1$ (`overview.md` §3.3) and at least $M-1$
  primes, so the strongest possible output $M\le\omega(N')+O(1)$ holds for the fakes' parameters at $W=O(M\log M)$
  (Corollary 6.9). Concretely: if an actual cluster existed with a fake's profile (rank $M-1$, Theorem 3.12 remark of
  the overview) and genuine angles, it would lie in no proper algebraic subgroup (Lemma 6.7), and no
  unlikely-intersection theorem would say anything about it.
* **The lens has algebraic fakes of its own.** The rank-one points $F_1,F_2$ of Proposition 6.6 are genuine algebraic
  points, lie in subgroups of codimension $n-1$ (far more "unlikely" than any cluster), have non-anomalous carriers
  $\ell_u$ (Lemma 6.3), and have height ratios $\to1$ (resp. $\to\frac n2$) above the cluster's $\frac23$ (resp.
  $\frac n3$). They play for this lens the role that angle models play for local inputs: they satisfy everything the
  lens can state and are not excluded by its conclusions. What they lack is exactly the equal-norm constraint.

## 6. First step: what is proved

All statements in this section are **proved here (unrefereed)**. Tori are split tori over $\overline{\mathbb Q}$;
"subtorus" means a connected algebraic subgroup; "coset" means a translate $gH$ of a subtorus by any point $g$.

### 6.1 Fixed carriers are anomalous

> **Lemma 6.1 (projection anomaly).** Let $V\subseteq\mathbb G_m^n$ be irreducible and let
> $\varphi:\mathbb G_m^n\to\mathbb G_m^m$ be a surjective homomorphism of tori. Let $V'=\overline{\varphi(V)}$. If
> $\dim V'<\dim V$ and $\dim V'<m$, then every point of $V$ lies on an anomalous subvariety of $V$; hence
> $V^{oa}=\varnothing$.

*Proof.* Let $v\in V$ and let $Y$ be an irreducible component of $V\cap\varphi^{-1}(\varphi(v))$ through $v$. By
upper semicontinuity of fibre dimension, $\dim Y\ge\dim V-\dim V'\ge1$. The fibre $\varphi^{-1}(\varphi(v))$ is a
coset of $\ker\varphi$, which is a finite union of cosets of the subtorus $H=(\ker\varphi)^0$, with
$\dim H=n-m$; since $Y$ is irreducible it lies in one coset $gH$. Now
$\max(0,\dim V+\dim H-n)=\max(0,\dim V-m)$, and $\dim Y\ge\dim V-\dim V'$ exceeds both $0$ and $\dim V-m$ by the
two hypotheses. So $Y$ is anomalous. $\square$

> **Theorem 6.2 (natural fixed carriers have empty $V^{oa}$).**
> (a) $V=X_M\cap\mathbb G_m^{2M}$, the torus part of the equal-norm cone (coordinates $u_x=z_x$, $v_x=\bar z_x$),
> projectivised in $\mathbb G_m^{2M}/\mathbb G_m$, has $V^{oa}=\varnothing$ for every $M\ge2$.
> (b) The chord space $V=L_M$, the closure in $\mathbb G_m^{\binom M2}/\mathbb G_m$ of
> $\{(z_a-z_b)_{a<b}:z\in\mathbb C^M\text{ with distinct coordinates}\}$, has $V^{oa}=\varnothing$ for every $M\ge4$.
> (c) $V=M_{0,k}$ in its Plücker torus $T_k=\mathbb G_m^{\binom k2}/\mathbb G_m^k$ (Plücker coordinates of $2\times k$
> matrices with pairwise independent columns, modulo column scaling; $\dim T_k=k(k-3)/2$) has
> $V^{oa}=\varnothing$ for every $k\ge5$.

*Proof.* (a) $V=\{u_xv_x=u_1v_1\ (x\ge2)\}$ is itself a subtorus $H$ of $\mathbb G_m^{2M-1}$, of dimension $M$. Take
$Y=V$: $\dim Y=M>\max(0,2M-(2M-1))=1$ for $M\ge2$. So $V$ is anomalous in itself.
(b) $V$ is irreducible of dimension $M-2$. Let $\varphi$ be the homomorphism to $\mathbb G_m^3/\mathbb G_m\cong\mathbb G_m^2$
keeping the three coordinates $(D_{12}:D_{13}:D_{23})$. Its image lies on the line $D_{12}-D_{13}+D_{23}=0$, so
$\dim V'=1<2=m$, and $1<M-2=\dim V$ for $M\ge4$. Lemma 6.1 applies.
(c) $V$ is the image of the irreducible variety $M_{0,k}$ under an injective morphism (the Plücker vector determines
the columns up to $\mathrm{PGL}_2$), so it is irreducible of dimension $k-3$. Let $\varphi:T_k\to T_4$ be induced by the coordinate projection keeping the six
Plücker coordinates of the first four columns; it is compatible with column scaling, hence well defined and
surjective. $\dim T_4=6-4=2$ (the scaling action of $\mathbb G_m^4$ on $\mathbb G_m^6$ has finite kernel), and the
cross-ratios of the four columns, which are characters of $T_4$, have exponent vectors of rank $2$ (exact computation,
`wild_checks/anomaly_ranks.out`; the full set of cross-ratios has rank $k(k-3)/2=\dim T_k$ for $k=4,\dots,9$). Its image is
$M_{0,4}=\{(\lambda,1-\lambda)\}$, a curve: $\dim V'=1<2$, and $1<k-3$ for $k\ge5$. Lemma 6.1 applies. $\square$

*Reading.* The obstruction is structural. Any carrier whose coordinates contain those of a sub-configuration
(three chords, four Plücker coordinates, the norm pair $u_x,v_x$) is fibred by "fix the sub-configuration", and
these fibres are anomalous. Habegger's bounded height theorem, and the Zilber–Pink conjecture in its anomalous
form, are therefore **vacuous** on every fixed carrier the record uses. This is the precise sense in which the
cluster's structure is "all of it anomalous": the arc condition is archimedean, not a membership in a subgroup.

### 6.2 The moving carrier $\ell_u$ and the cluster's heights

> **Lemma 6.3 ($\ell_u$ is not in a coset).** Let $u\in\mathbb P^n(\mathbb Q(i))$, $n\ge1$, have nonzero
> coordinates with $u_a/u_b\notin\mathbb R$ for $a\ne b$. Then $X,Y$ are $\mathbb Q$-independent, and
> $\ell_u\cap\mathbb G_m^n$ is not contained in any translate of a proper subtorus.

*Proof.* If $X,Y$ were dependent, all $u_x$ would be real multiples of one complex number, so $u_a/u_b\in\mathbb R$.
Parametrise $\ell_u\cong\mathbb P^1$ by $(s:t)\mapsto sX+tY$. The $x$-th coordinate is the linear form
$L_x=sX_x+tY_x\not\equiv0$, with zero $\zeta_x=(-Y_x:X_x)$. $L_a$ and $L_b$ are proportional iff
$(X_a,Y_a)\parallel(X_b,Y_b)$ iff $u_a/u_b\in\mathbb R$; so the $\zeta_x$ are pairwise distinct. A coset of a proper
subtorus is contained in $\{\prod_xL_x^{m_x}=c\}$ for some $m\ne0$ with $\sum m_x=0$ (in the projective torus).
On $\mathbb P^1$ the rational function $\prod L_x^{m_x}$ has divisor $\sum_xm_x[\zeta_x]$, which is nonzero because the
$\zeta_x$ are distinct and $m\ne0$; so it is not constant on $\ell_u$. $\square$

So BMZ's curve theorem applies to each $\ell_u$, but with a constant depending on $\ell_u$.

> **Proposition 6.4 (the cluster's heights).** Let $z=(z_1,\dots,z_M)$, $M\ge2$, be a normalised cluster on an arc
> of length $\le C\sqrt R$ with $R>C^2/\pi^2$, viewed as $z\in\mathbb P^{M-1}(\mathbb Q(i))$. Then:
> (i) $z_a/z_b\notin\mathbb R$ for $a\ne b$;
> (ii) $h(z)=\log R$;
> (iii) $h(\ell_z)\le\frac32\log R+\log C$;
> (iv) for every anchor $a$, $h^{(a)}_\Sigma(z)\ge(M-1)\big(\frac12\log R-\log C\big)$.
> Hence $h(z)/h(\ell_z)\ge\log R/(\frac32\log R+\log C)\to\frac23$ and
> $h^{(a)}_\Sigma(z)/h(\ell_z)\ge(M-1)(\frac12\log R-\log C)/(\frac32\log R+\log C)\to\frac{M-1}3$.

*Proof.* The angular span of the arc is $\le C\sqrt R/R=CR^{-1/2}$. (i) $z_a/z_b\in\mathbb R$ with $|z_a|=|z_b|$ and
$z_a\ne z_b$ forces $z_a=-z_b$, at angular distance $\pi>CR^{-1/2}$. (ii) For Gaussian-integral coordinates with
Gaussian gcd 1, every finite place has $\max_x|z_x|_v=1$, and the unique (complex) archimedean place contributes
$\frac12\cdot2\log\max|z_x|=\log R$. (iii) $|\mathrm{Im}(\bar z_az_b)|=N|\sin(\theta_b-\theta_a)|\le N\cdot CR^{-1/2}=CR^{3/2}$,
and the gcd of the Plücker coordinates is $\ge1$. (iv) $z_x/z_a=\beta/\delta$ with $\beta=z_x/g_{xa}$,
$\delta=z_a/g_{xa}$ coprime, so $h(z_x/z_a)=\log|\beta|=\frac12d_{xa}\ge\frac W4-\log C=\frac12\log R-\log C$ by the
pair margin bound (§1). $\square$

> **Remark 6.5 (where the thresholds sit).** For every $u$, $|\mathrm{Im}(\bar u_au_b)|\le|u_a||u_b|$ gives
> $h(\ell_u)\le2h(u)$, i.e. $h(u)/h(\ell_u)\ge\frac12$ always. The same computation as (iii) on an arc of length
> $CR^\alpha$ gives $h(\ell_z)\le(1+\alpha)\log R+\log C$, ratio $\ge1/(1+\alpha)$: full circles ($\alpha=1$) give
> $\frac12$ (the exact data give exactly $0.5000$ on four full circles, `cluster_heights.out`), short arcs
> ($\alpha=\frac12$) give $\frac23$, and ratio $\to1$ corresponds to "arcs of bounded length". So a height comparison
> with $\kappa<\frac23$ is a statement at exactly the Dirichlet exponent $\alpha=\frac12$ of `round4/hardness.md`.
> On every actual cluster computed (best $k$-windows, $k=4..10$, three norms) the gcd of the Plücker coordinates is
> $1$ and the ratio equals the bound of (iii) to four digits (**evidence**, §7).

### 6.3 Uniform-in-the-line height comparisons are false at the needed constants

> **Proposition 6.6 (rank-one points beat the cluster).** Let $A\ge2$, $n\ge2$, $K\ge n$ be integers and
> $\gamma=A+i$. Put
> $$F_1=(1:\gamma:\gamma^2:\dots:\gamma^{n-1}:\gamma^K),\qquad F_2=(1:\gamma^{K-n+1}:\gamma^{K-n+2}:\dots:\gamma^K)\in\mathbb P^n(\mathbb Q(i)).$$
> (a) $F_1,F_2\in\mathcal U_n(0)$; their coordinate ratios generate a group of rank one, so they lie in algebraic
> subgroups of codimension $n-1$.
> (b) $h(F_1)\ge\frac K{K+n-1}\,h(\ell_{F_1})$.
> (c) $h^{(0)}_\Sigma(F_2)\ge\frac{n(K-n+1)}{2K}\,h(\ell_{F_2})$.
> (d) For $\psi=\alpha h+\beta h^{(0)}_\Sigma$ with $\alpha,\beta\ge0$, $(\alpha,\beta)\ne0$, one of the two families has
> $\liminf_{K\to\infty}\psi/h(\ell)\ge\frac23\alpha+\frac n3\beta$, the cluster's guaranteed ratio of Proposition 6.4.
> Consequently (UI-$\ell$) is false for every $\kappa<1$, (UI-$\ell_\Sigma$) is false whenever
> $\kappa_\Sigma(n)<\frac n2$ for some $n\ge2$, and every inequality $\psi(u)\le\kappa_\psi h(\ell_u)+c$ on
> $\mathcal U_n(\varepsilon)$ (any $\varepsilon\ge0$, any $c$) with $\kappa_\psi<\frac23\alpha+\frac n3\beta$ is false.

*Proof.* (a) $\gamma/\bar\gamma=(A^2-1+2Ai)/(A^2+1)$ has norm one and is not in $\mu_4=\{\pm1,\pm i\}$ for $A\ge2$
(it equals $\pm1$ only if $A=0$, and $\pm i$ only if $A=1$); $\mathbb Q(i)$ contains no other roots of unity, so
$\gamma/\bar\gamma$ has infinite order and $\gamma^m\notin\mathbb R$ for $m\ge1$. Hence all coordinate ratios
$\gamma^{a-b}$ ($a\ne b$) are non-real, the coordinates are nonzero, and the ratios lie in $\langle\gamma\rangle$, of
rank one. (b) The coordinates are Gaussian integers including $1$, so the Gaussian gcd is $1$ and
$h(F_1)=\log|\gamma^K|=K\log|\gamma|$ (Proposition 6.4(ii)'s computation). Every Plücker coordinate satisfies
$|\mathrm{Im}(\bar u_au_b)|\le|u_a||u_b|\le|\gamma|^{K+n-1}$ and their gcd is $\ge1$, so
$h(\ell_{F_1})\le(K+n-1)\log|\gamma|$. (c) For $\alpha\in\mathbb Z[i]\setminus0$, $h(\alpha)=\log|\alpha|$, so
$h^{(0)}_\Sigma(F_2)=\sum_{k=K-n+1}^Kk\log|\gamma|\ge n(K-n+1)\log|\gamma|$, while $h(\ell_{F_2})\le2K\log|\gamma|$.
(d) For $F_1$: $h^{(0)}_\Sigma(F_1)\ge h(\gamma^K)=h(F_1)$, so $\psi(F_1)/h(\ell)\ge(\alpha+\beta)K/(K+n-1)\to\alpha+\beta$.
For $F_2$: $h(F_2)=K\log|\gamma|\ge\frac12h(\ell_{F_2})$, so $\psi(F_2)/h(\ell)\ge\frac\alpha2+\beta\frac{n(K-n+1)}{2K}\to\frac{\alpha+n\beta}2$.
If $\alpha+\beta<\frac{2\alpha+n\beta}3$ then $\alpha<(n-3)\beta$, and then
$\frac{\alpha+n\beta}2-\frac{2\alpha+n\beta}3=\frac{n\beta-\alpha}6>\frac{3\beta}6>0$. The consequences follow
because the inequality would have to hold at $F_1$ or $F_2$ for all large $K$. $\square$

The exact values (`cluster_heights.out`) are $h(F_1)/h(\ell_{F_1})=1.000000$ for $u=(1,\gamma,\gamma^2)$ (indeed
$h(u)=h(\ell_u)=\log(A^2+1)$ exactly, since the Plücker coordinates are $1$, $2A$, $A^2+1$), $0.9906$ at
$(A,n,K)=(10,20,2000)$, and $h^{(0)}_\Sigma(F_2)/h(\ell)=14.93$ at $(2,30,3000)$ against $n/3=10$.

*What Proposition 6.6 says.* The cluster's anomaly (ratio $\frac23$, resp. $\frac n3$) is **smaller** than what
rank-one points achieve without the equal-norm constraint (ratio $\to1$, resp. $\to\frac n2$). So no statement whose
hypotheses are "$u$ lies on $\ell_u$ and in (a cone around) an algebraic subgroup" can separate clusters from
these points; a height comparison that excludes clusters must use the equal-norm constraint
$|u_0|=\dots=|u_n|$. With that constraint the carrier is $X_M$ (Theorem 6.2(a), anomalous) or the Ptolemy point of
$M_{0,k}$ (Theorem 6.2(c), anomalous), and the comparison becomes a height inequality on those varieties of the
Vojta / $(L_k)$ type (a reading of the record, not a theorem), which is not an unlikely-intersection statement.

*Not covered by Proposition 6.6.* The anchor-free variant $\min_ah^{(a)}_\Sigma(u)\le\kappa h(\ell_u)+c$ is not refuted by
$F_1,F_2$ (with an anchor inside the block, $F_2$ has min-anchor sum height $\approx K\log|\gamma|$). It is
not a UI-type statement: its strength comes from requiring **every** pairwise ratio to have large height, i.e. from
the pair margins $d_{xy}\ge W/2-O(1)$, and it constrains unequal-norm tuples as well. Call it **(S$_{\min}$)**: for some
$\delta>0$, every $u\in\mathcal U_n(0)$ has $\min_ah^{(a)}_\Sigma(u)\le(\frac13-\delta)n\,h(\ell_u)+o(n^2\log n)$. It is open
and non-standard, it covers strictly more tuples than the circle problem, and no technique is known for it; it is
recorded in §8 only as the residue of the lens.

### 6.4 Even an ideal unlikely-intersection input stops at $M\le\omega(N')+O(1)$

> **Lemma 6.7 (full rank means no subgroup).** Let $z$ be a normalised cluster with exponent rows $a_x\in\mathbb Z^r$
> ($r=\omega(N')$) and let $B$ be the $(M-1)\times r$ matrix with rows $a_x-a_1$. The point
> $(z_x/z_1)_{x\ge2}\in\mathbb G_m^{M-1}(\mathbb Q(i))$ lies in a proper algebraic subgroup iff $\mathrm{rank}\,B\le M-2$.
> In particular it lies in none when $\mathrm{rank}\,B=M-1$ (which needs $r\ge M-1$), and it always lies in one when
> $M\ge r+2$.

*Proof.* $z_x/z_1=(\varepsilon_x/\varepsilon_1)\prod_j(\pi_j/\bar\pi_j)^{a_{xj}-a_{1j}}$. A proper algebraic subgroup is
contained in $\ker\chi_m$ for some $m\in\mathbb Z^{M-1}\setminus0$, and $\chi_m(z)=\prod_x(z_x/z_1)^{m_x}$ equals a unit
times $\prod_j\pi_j^{c_j}\bar\pi_j^{-c_j}$ with $c=B^{\mathsf T}m$. By unique factorisation in $\mathbb Z[i]$ (the
$\pi_j,\bar\pi_j$ are pairwise non-associate), this is $1$ only if $c=0$. If $\mathrm{rank}\,B=M-1$ then $c\ne0$ for
every $m\ne0$, so $z$ lies in no proper subgroup. If $\mathrm{rank}\,B\le M-2$ pick $m\ne0$ with $B^{\mathsf T}m=0$;
then $\chi_m(z)\in\mu_4$ and $z\in\ker\chi_{4m}$, of codimension 1. $\square$

> **Proposition 6.8 (what the lens can give).** Assume (UI-$\ell$) with $\kappa<\frac23$, $\varepsilon=0$ and
> $c(n)=o(n\log n)$, or (UI-$\ell_\Sigma$) with $\kappa_\Sigma(n)\le(\frac13-\delta)n$ and $c(n)=o(n^2\log n)$. Then for
> each $C$ there are $M_1(C)$ and $R_0(C)$ such that every normalised cluster on an arc of length $\le C\sqrt R$
> with $R>R_0(C)$ and $\mathrm{rank}\,B\le M-2$ has $M\le M_1(C)$. Consequently every such cluster satisfies
> $$M\le\max\big(M_1(C),\ \omega(N')+1\big).$$

*Proof.* By Lemma 6.7 and Proposition 6.4(i), $z\in\mathcal U_{M-1}(0)$. Under (UI-$\ell$), Proposition 6.4(ii)–(iii)
give $\log R\le\kappa(\frac32\log R+\log C)+c(M-1)$, i.e. $(1-\frac32\kappa)\log R\le\kappa\log C+c(M-1)$. The proved
inequality (Theorem 2.1 of the overview, with $W'=W$ for a normalised cluster) gives
$3M\log M\le2\log R+(14+4\log^+C)M$. Combining,
$3M\log M\le\frac{2}{1-3\kappa/2}\big(\kappa\log C+c(M-1)\big)+(14+4\log^+C)M$, which fails for $M\ge M_1(C)$ since
$c(n)=o(n\log n)$. Under (UI-$\ell_\Sigma$), Proposition 6.4(iv) gives
$(M-1)(\frac12\log R-\log C)\le(\frac13-\delta)(M-1)(\frac32\log R+\log C)+c(M-1)$, so
$\frac32\delta\log R\le\frac43\log C+c(M-1)/(M-1)$, and the same comparison with Theorem 2.1 fails for large $M$
because $c(n)/n=o(n\log n)$. Finally, if $M\ge\omega(N')+2$ then $\mathrm{rank}\,B\le\omega(N')\le M-2$. $\square$

> **Corollary 6.9 (the ceiling of the lens).** Every conclusion that the unlikely-intersection lens can deliver
> through Lemma 6.7 is of the form "rank-deficient clusters are bounded", i.e. $M\le\omega(N')+O_C(1)$. This is not a
> growth improvement: $\omega(N')\le(2+o(1))\log R/\log\log R$ is the only general bound, which is weaker than the
> proved $(\frac23+\varepsilon)\log R/\log\log R$, and the profiles of all the record's fakes have
> $\omega\ge M-1$ (full rank $M-1$, `overview.md` §3.3, Theorem 3.12 remark), so the conclusion is satisfied by
> every fake profile at $W=O(M\log M)$. It is the same ceiling as the Subspace route QS (`lit18` Prop B:
> "at best $M\le\omega(N')$").

*Remark (cones).* Rémond's cones enlarge the class $\mathcal U_n(\varepsilon)$, so Proposition 6.6 still refutes the
cone versions. They could in principle reach full-rank clusters by discarding primes of small total weight; but
lowering the rank from $\omega(N')$ to $M-2$ means discarding $\omega(N')-M+2$ primes, and the only proved bound is
$\omega(N')\le(\frac2\varepsilon+o(1))M$ for near-extremal clusters (`overview.md` Thm 3.11), so the discarded weight is
not $o(W)$ in general and a fixed small cone parameter does not suffice. (**Heuristic** for the general case; not
used.)


## 7. Computations

All scripts are in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round20/wild_checks/`, run
under `prlimit --as=3500000000`, exact integer / Gaussian-integer arithmetic (logs taken only at the end).

| script → output | what it checks | label |
|---|---|---|
| `anomaly_ranks.py` → `anomaly_ranks.out` | rank of all cross-ratio exponent vectors in $\mathbb Z^{\binom k2}$ is $k(k-3)/2$ for $k=4..9$; cross-ratios on 4 (resp. 5) points have rank 2 (resp. 5); the anomaly inequality of Lemma 6.1 for forgetful fibres of $M_{0,k}$ ($k=5..9$) and for triangle fibres of $L_M$ ($M=4..11$) | **exact** (used in Thm 6.2(c) for the rank-2 fact) |
| `line_examples.py` → `line_examples.out` | $h(u)$, $h(\ell_u)$ for $u=(1,\gamma,\gamma^b)$, $\gamma=A+i$: ratio exactly $1$ for $b=2$ ($h(u)=h(\ell_u)=\log(A^2+1)$), $\to1$ for $b=3,5$; $(\gamma^0,\dots,\gamma^8)$ gives $4/7$ | **exact** values; illustrate Prop 6.6 |
| `cluster_heights.py` → `cluster_heights.out` | Prop 6.6 families $F_1$, $F_2$ against their rigorous bounds; four full circles ($N=1105,\dots,48612265$, $M=32..256$): ratio exactly $0.5000$, Plücker gcd $1$, rank $=\omega(N)$; best $k$-windows on $N=1176852625$ ($k=6,8,10$; the referee's 10-point cluster has $C=8.2386$ by chord), $N=3159797225$ ($k=5,6,7$), $N=2576450045$ ($k=4,5,6$): Plücker gcd $G=1$ in every case, $h(z)/h(\ell_z)$ equal to the bound $\log R/(\frac32\log R+\log C)$ to 4 digits (0.583–0.633), min-anchor sum ratios $1.03$–$2.76$ against $n/3$ | **evidence** (no proof uses the cluster rows) |

Two data points worth recording (**evidence**): the 10-point cluster on $1176852625=5^3\,13^2\,17\,29\,113$ has
multiplicative rank $5<9=M-1$, so it is rank-deficient and lies in an algebraic subgroup of codimension 4 (Lemma 6.7);
the 7-point window on $3159797225$ has rank $6=M-1$ (full rank), so it lies in no proper algebraic subgroup and is
invisible to the lens.


## 8. Verdict and open problem

**Status: reduced to a cap (new barrier facts; no growth, no uniform bound).**

1. The most promising mechanism not in the record, chosen because it is the only mature unconditional Diophantine
   theory with constants **uniform in the finitely generated group** (hence immune to items 4–5 by design), is the
   Bombieri–Masser–Zannier / Habegger / Rémond theory of unlikely intersections in tori.
2. It needs a carrier. **Fixed carriers are vacuous**: $X_M$, the chord space $L_M$ ($M\ge4$) and $M_{0,k}$ ($k\ge5$) in
   its Plücker torus all have empty open-anomalous locus (Theorem 6.2, via the projection-anomaly Lemma 6.1).
3. The **moving carrier** is the real line $\ell_z=\mathrm{span}(\mathrm{Re}\,z,\mathrm{Im}\,z)$, which encodes the arc
   condition exactly in its height ($h(z)/h(\ell_z)\ge\frac23-o(1)$, the Dirichlet exponent; Prop 6.4, Remark 6.5). A
   uniform-in-the-line height comparison strong enough to exclude clusters is **false**: rank-one points reach ratio
   $1$ (projective) and $\frac n2$ (chart height), above the cluster's $\frac23$ and $\frac n3$, and every nonnegative
   mixture is beaten too (Prop 6.6). The equal-norm constraint is what separates clusters, and with it the carrier is
   $X_M$ or $M_{0,k}$ again, where the needed inequality is Vojta / $(L_k)$, not an unlikely-intersection statement.
4. Even an ideal unlikely-intersection input acts only on rank-deficient clusters (Lemma 6.7), so its ceiling is
   $M\le\omega(N')+O_C(1)$ (Prop 6.8, Cor 6.9), which is no growth improvement, is satisfied by every fake profile,
   and coincides with the ceiling of the Subspace route QS (lit18 Prop B).

**New entries for the barrier map** (proposed; unrefereed):
* *Row "unlikely intersections (BMZ, Habegger's bounded height, Rémond cones, Zilber–Pink in anomalous form)":*
  kind V + "algebraic fakes"; fixed carriers have $V^{oa}=\varnothing$ (Thm 6.2); moving-line comparisons are refuted
  by rank-one points (Prop 6.6); ceiling $M\le\omega(N')+O(1)$ (Cor 6.9).
* *Lesson for item 5:* uniformity in $S$ is not enough. The lit18 recommendation "use equality of norms or a threshold
  tied to $\prod p$" is sharpened: an $S$-uniform tool that does not use the equality of norms is beaten by rank-one
  algebraic points (Prop 6.6), and one that uses only multiplicative dependence stops at $\omega(N')$ (Cor 6.9).

**Residual open statement (non-standard; not recommended as a target).** (S$_{\min}$) of §6.3: for some $\delta>0$,
every $u\in\mathbb P^n(\mathbb Q(i))$ with pairwise non-$\mathbb R$-proportional coordinates that are multiplicatively
dependent satisfies $\min_a\sum_{x\ne a}h(u_x/u_a)\le(\frac13-\delta)\,n\,h(\ell_u)+o(n^2\log n)$. It would bound
rank-deficient clusters (and so give only $M\le\omega(N')+O(1)$ by Cor 6.9); it covers unequal-norm tuples; the
examples of Prop 6.6 do not refute it. A refutation attempt (unequal-norm tuples of $S$-units, all pairwise ratios of
large height, within relative distance $H^{-1/2}$ of a common ray) is the cheap direction; a heuristic count
suggests such tuples are rarer than circle clusters because the equal-norm condition must then be met approximately
(**heuristic**).

**Recommendation.** Do not open an internal round on this lens. Add the row above to the barrier map. The record's
conclusion stands: the needed statement is the equal-norm height inequality at the Dirichlet exponent, i.e. Vojta on
$\mathrm{Bl}_YX_7$ or $(L_8)$; uniformity in $S$ alone does not reach it.

## Appendix: other wildcard lenses screened (short reasons; unrefereed)

| lens | why it was not chosen | closest checklist item |
|---|---|---|
| o-minimal / Pila–Wilkie counting in families | the cluster's exponent vectors are integer points of a semilinear slab $\{t:\langle2t-e,\phi\rangle\in\theta_0+[0,\Delta]\}$; the algebraic part is everything, and the algebraicity of $\phi_j=\arg\pi_j$ is invisible to definable geometry; Galois orbits have degree 2 | 4 (linear forms), 1 |
| Chabauty–Skolem / $p$-adic interpolation of exponent vectors | needs rank $<$ dimension; here the torus is 1-dimensional and the rank is $\omega(N')$; the arc condition is an archimedean inequality, and the available equalities are $S$-unit equations $z_x-z_y=g_{xy}c_{xy}$ | 5 |
| Erdős–Turán / Bilu equidistribution of roots | integer polynomials may have all roots at one point ($(X-1)^D$); for $\prod(\bar z_xX-z_x)$ the Erdős–Turán ratio is $2^M$, trivial | — |
| Elekes–Szabó / sum–product | the circle group is the exceptional group case | 1 |
| quartic / Hilbert reciprocity | $[\pi_j/\varpi]_4=[\varpi/\pi_j]_4$ only re-expresses genuine prime residues; subsumed by P3-G (genuine residues) | 1 (P3 classes) |
| Frey curves / modular method | conductor $\approx\mathrm{rad}$; all relevant relations are squarefree | 15 (abc vacuity) |
| SOS / Positivstellensatz certificates over $\mathbb Z[i]$ | integrality enters only as "$|w|\ge1$ for a chosen nonzero Gaussian integer $w$", i.e. an auxiliary polynomial with forced divisibility | 3, 6 |
| ergodic / measure rigidity, Hecke equidistribution | statements about measures or almost all $N$; a cluster is $M\ll r_2(N)$ points | 2 |
| geometry of numbers in the $C\sqrt R\times O(C^2)$ strip | $M\le2(2/\lambda_2+1)$ with $\lambda_2\gtrsim(C^3\sqrt R)^{-1/2}$, i.e. $M\ll R^{1/4}$; bounded only near directions of small height (Chan-type) | — |
