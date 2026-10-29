# McKinnon-Roth, full-text check against Question B

Read in full: all 55 pages of `f5856c6f-1306.2976v2.pdf` and all 17 pages of `75f151c0-1306.2977v1.pdf`.
Page numbers below are the numbers printed on the pages; they match the PDF page index. Text in quotation marks is
copied from the PDFs (mathematical symbols transcribed). Anything marked **[inference]** is my own deduction and is
not stated in either paper.

## 1. Bibliographic data as printed in the PDFs

* **[MR-I]** David McKinnon and Mike Roth, *Seshadri constants, diophantine approximation, and Roth's theorem for
  arbitrary varieties*. Stamp: "arXiv:1306.2976v2 [math.AG] 27 Apr 2015". 55 pages. "2010 Mathematics Subject
  Classification. Primary 14G05; Secondary 14G40." Both authors are funded by NSERC; affiliations are Waterloo and
  Queen's. **The PDF does not print the journal reference** (Invent. Math. 200 (2015) 513-583), so these files can
  neither confirm nor contradict it. Reference [15] in it is the second paper, listed as "submitted".
* **[MR-L]** David McKinnon and Mike Roth, *An analogue of Liouville's theorem and an application to cubic surfaces*.
  Stamp: "arXiv:1306.2977v1 [math.AG] 12 Jun 2013"; footnote "Date: Revised on January 18, 2013". 17 pages, including
  three appendix tables (pp. 13-16). **The PDF does not print the journal reference** (Eur. J. Math. 2 (2016)). This
  is arXiv **v1**, so theorem numbers in the published version may differ. [MR-I] p. 24 cites "[15, Theorem 3.3]",
  which matches the v1 numbering.

## 2. Definitions

**Distance** ([MR-I] §2 pp. 8-12; [MR-L] p. 2).
* There is one fixed place $v$ of $k$, extended to $\bar k$.
* Archimedean case ([MR-I] (2.1)): embed $X\hookrightarrow\mathbb P^r_k$ and set
  $d_v(x,y)=\big(1-\frac{|\sum x_i\bar y_i|^2}{\sum|x_i|^2\sum|y_j|^2}\big)^{[k_v:\mathbb R]/2}$.
* Non-archimedean case ([MR-I] (2.2)): $d_v(x,y)=\frac{\max_{i<j}\|x_iy_j-x_jy_i\|_v}{\max\|x_i\|_v\max\|y_j\|_v}$.
* Absolute values are normalized relative to $k$.
* By Prop 2.4, different embeddings give equivalent distances.
* By Lemma 2.5 and Lemma 2.6, near $x$ the distance $d_v(x,\cdot)$ is equivalent to
  $\max\|u_i(\cdot)\|_v$, where the $u_i$ generate the maximal ideal of $x$.

**Height** ([MR-I] p. 13). "multiplicative and relative to $k$"; $H(x)=\prod_v\max(\|x_0\|_v,\dots,\|x_r\|_v)$ for
$\mathcal O_{\mathbb P^r}(1)$. Heights are defined only up to equivalence. By Remark (b) on p. 13, $\alpha$ depends only
on the quasi-equivalence class, hence only on the class of $L$ in $NS$.

**Approximation constants** ([MR-I] p. 13, Defs 2.7-2.9).
* "Definition 2.7. Let $X$ be a projective variety, $x\in X(\bar k)$, $L$ a line bundle on $X$. For any sequence
  $\{x_i\}\subseteq X(k)$ of distinct points with $d_v(x,x_i)\to0$ … we set
  $A(\{x_i\},L)=\{\gamma\in\mathbb R \mid d_v(x,x_i)^\gamma H_L(x_i)\text{ is bounded from above}\}$."
* "Definition 2.8. … $\alpha_x(\{x_i\},L)$ [is] the infimum of $A(\{x_i\},L)$ (in particular if $A=\varnothing$ then
  $\alpha_x(\{x_i\},L)=\infty$)."
* "Definition 2.9. The approximation constant $\alpha_{x,X}(L)$ of $x$ with respect to $L$ is defined to be the
  infimum of all approximation constants of sequences of points in $X(k)$ converging to $x$."
* Equivalent form (Def 2.10 and Prop 2.11): if $L$ has the Northcott property, then $\alpha_x(L)=\sup B_x(L)$, where
  $B_x(L)$ is the set of $\gamma\ge0$ such that for every $C>0$ only finitely many $x_i\in X(k)$ satisfy
  $d_v(x,x_i)^\gamma H_L(x_i)<C$.
* [MR-L] Def 2.1 (p. 3) gives the same definition for $L$ ample and $x\in X(\bar k)$: "$\alpha_x$ … the unique
  extended real number $\in(0,\infty]$ such that $d_v(x,y)^{\alpha_x+\delta}<H_L(y)^{-1}$ has only finitely many
  solutions $y\in X(k)$ (respectively infinitely many) for any $\delta<0$ (respectively $\delta>0$)."
* $\alpha_{x,Z}(L|_Z)$, the constant "computed on a subvariety", is the same definition for $Z$, with
  $\alpha_{z,Z}(L|_Z)\ge\alpha_{z,X}(L)$ (Prop 2.14(c)). The phrase "$\alpha_x(L)$ is computed on a proper subvariety"
  means $\alpha_{x,X}(L)=\alpha_{x,Z}(L|_Z)$ (p. 4).
* **No "essential" or Zariski-dense constant is defined in either paper.** Essential-type statements appear only as
  theorems with an exceptional proper subvariety $Z$ or base locus $B$ (see §3).

**Seshadri constant** ([MR-I] Def 3.1, p. 20; [MR-L] Def 2.9, p. 5).
* "Let $X$ be a projective variety over $\mathrm{Spec}(k)$, $x$ a point of $X(\bar k)$, and $L$ a nef line bundle on
  $X$. … $\epsilon_{x,X}(L):=\sup\{\gamma\ge0\mid\pi^*L-\gamma E\text{ is nef}\}$ where
  $\pi:\widetilde X\to X_{\bar k}$ is the blowup … at $x$."
* Prop 3.2 gives the curve form $\epsilon_x(L)=\inf_{C\ni x} (L\cdot C)/\mathrm{mult}_x C$.
* The Seshadri constant is geometric (it depends only on $X_{\bar k}$).
* The cover version is $\hat\epsilon_x(L)=\sup_{\varphi:Y\to X}\min_{y\in\varphi^{-1}(x)}\epsilon_y(\varphi^*L)$, with
  the supremum over unramified covers (Def 8.2, p. 43).

**$\gamma_{\rm eff}$ and $\beta$** ([MR-I] §4, p. 25; Def 4.3, p. 28).
* $L_\gamma:=\pi^*L-\gamma E$, and "$\gamma_{\rm eff}=\gamma_{{\rm eff},x}(L)=\sup\{\gamma\ge0\mid L_\gamma\text{ is effective}\}$".
* $f(\gamma)=\mathrm{Vol}(L_\gamma)/\mathrm{Vol}(L)$ and
  "$\beta_x(L)=\int_0^\infty f(\gamma)\,d\gamma=\int_0^{\gamma_{\rm eff}}f(\gamma)\,d\gamma$".
* For reducible $X_{\bar k}$, $\beta$ is the minimum over the components through $x$ (Def 4.3).
* Interpretation (§9, p. 50): "$\beta_x(L)$ is the expected order of vanishing at $x$ of a section of $L$".
* $f(\gamma_{\rm eff})=0$ by continuity (p. 25).
* For nef $L$: $\mathrm{Vol}(L_\gamma)=\mathrm{Vol}(L)-(\mathrm{mult}_xX)\gamma^n$ on $[0,\epsilon_x(L)]$ (Remark, p. 27).
* $\hat\beta$ is the cover version (Def 8.2).
* The record's $\tau_x(L)$ is the same number as MR's $\gamma_{{\rm eff},x}(L)$ (the big cone is the interior of
  the effective cone, p. 25). The record's expectation constant $\beta$ is MR's $\beta_x(L)$.

**Best approximation and Seshadri curves** ([MR-L] Def 4.1, p. 8; Def 4.3, p. 9).
* "A sequence $\{x_i\}\to x$ whose approximation constant is equal to $\alpha_x(L)$ … is called a sequence of best
  approximation to $x$. A curve $C$ passing through $x$ is … a curve of best approximation … if $C$ contains a sequence
  of best approximation."
* A Seshadri curve is a curve with $\epsilon_{x,X}(L)=(L\cdot C)/\mathrm{mult}_x(C)$.

## 3. Theorems, with exact hypotheses

Conventions in [MR-I] (p. 7): $k$ is a number field. "Variety over Spec(k)" means a possibly reducible, possibly
singular projective variety, i.e. a reduced projective $k$-scheme. A "sequence of $k$-points" means an infinite
sequence of **distinct** points of $X(k)$. Approximating points are always $k$-rational. The target $x$ may be
$k$-rational or algebraic.

### 3.1 Liouville-type theorems (target may be rational)

* **[MR-I] Introduction, p. 3, "Roth theorems".** "If $X$ is a variety over Spec(k), and $x\in X(\bar k)$ with field
  of definition $K$, then for any ample line bundle $L$ on $X$ we have $\alpha_x(L)\ge\frac1d\epsilon_x(L)$, where
  $d=[K:k]$. … we regard this as the general version of Liouville's theorem."
* **[MR-I] p. 24 (Arakelov discussion).** For $x\in X(k)$: "If $0<\gamma<\epsilon_x(L)$ then the line bundle
  $L-\gamma E$ is ample, and therefore satisfies the Northcott property. … This proves the inequality
  $\alpha_x(L)\ge\epsilon_x(L)$". The paper calls this "a stronger inequality than
  $\alpha_x(L)\ge\frac12\epsilon_x(L)$". It also says: "(See also [15, Theorem 3.3] for an extension involving
  $\gamma$ in the big cone and the asymptotic base locus.)"
* **[MR-L] Lemma 3.1 (p. 6).** For $x\in X(k)$, embedded at $[1:0:\dots:0]$:
  $e_w(y)=\min(1,\max_i|u_i(y)|_w)$ satisfies (a) $e_w\le1$, (b) $e_v$ is equivalent to $d_v$, and
  (c) $H_E(y)=(\prod_w e_w(y))^{-1}$ for $y\ne x$.
* **[MR-L] Lemma 3.2 (p. 6).** "Suppose that $x\in X(k)$ and let $\pi:\widetilde X\to X$ be the blow up at $x$ with
  exceptional divisor $E$. Let $L$ be an ample line bundle on $X$ and $\gamma>0$ a rational number such that
  $L_\gamma:=\pi^*L-\gamma E$ is in the effective cone of $\widetilde X$. Let $B'$ be the asymptotic base locus of
  $L_\gamma$ and set $B=\pi(B')$. Then for any sequence $\{x_i\}\to x$ such that all points of $\{x_i\}$ are outside of
  $B$, $\alpha(\{x_i\},L)\ge\gamma$."
  * Proof (pp. 6-7): $H_{L_\gamma}\ge c$ on $(\widetilde X\setminus B')(k)$, and
    $H_{L_\gamma}=H_LH_E^{-\gamma}=H_L\prod_we_w^\gamma\le H_Le_v^\gamma$.
* **[MR-L] Theorem 3.3 (p. 7), the main Liouville theorem.** "Let $X$ be an algebraic variety defined over $k$,
  $x\in X(\bar k)$ any point, and set $d=[K:k]$ where $K$ is the field of definition of $x$. Let
  $\pi:\widetilde X\to X$ be the blowup of $X$ at $x$, with exceptional divisor $E$, $L$ an ample line bundle on $X$,
  and $\gamma>0$ a rational number such that $L_\gamma:=\pi^*L-\gamma E$ is in the effective cone of $\widetilde X$.
  Finally let $B'$ be the asymptotic base locus of $L_\gamma$ and set $B=\pi(B')$. Then
  * (a) For any sequence $\{x_i\}\to x$ of $k$-points approximating $x$ if infinitely many points of $\{x_i\}$ are
    outside $B$ then $\alpha(\{x_i\},L)\ge\gamma/d$.
  * (b) If $\alpha_x(L)<\gamma/d$ then $x\in B$ and $\alpha_x(L)=\alpha_x(L|_B)$.
  * (c) If $x\in B$ and $\alpha_x(L|_B)\ge\gamma/d$ then $\alpha_x(L)\ge\gamma/d$."
  * "Remark. Theorem 3.3 still holds if we replace $B$ by the Zariski closure of $B(k)$."
* **[MR-L] Corollary 3.4 (p. 7).** "For all ample line bundles $L$ on $X$ we have $\alpha_x(L)\ge\epsilon_x(L)/d$."
* **[MR-L] p. 8, non-split quadric example.** This is the only non-split example in either paper. Take a smooth
  quadric $X\subset\mathbb P^3$ with no $k$-lines, $L=\mathcal O(1)|_X$ and $x\in X(k)$. Then $\epsilon_x(L)=1$, but
  $\pi^*L-2E$ is effective, with base locus the two (conjugate, non-$k$) lines through $x$. So $B(k)=\{x\}$, and Thm
  3.3(a) gives $\alpha_x(L)=2$ for all $k$-points. Here Galois non-rationality of the base locus is used to raise
  $\alpha_x$ (the infimum) from $\epsilon=1$ up to $\gamma=2$. The value $\gamma$ is still limited to
  $L_\gamma$ effective.

### 3.2 Roth-type theorems ([MR-I] §§5-6)

* **Theorem 5.1 (p. 29), the central theorem.** Setting: $S$ a finite set of places, each extended to $\bar k$;
  $X$ irreducible of dimension $n$ over $k$; $x_v\in X(\bar k)$ for $v\in S$; $L$ an ample $\mathbb Q$-bundle.
  "Given a collection $\{R_v\}_{v\in S}$ of positive real numbers, if (5.3) $\sum_{v\in S}\beta_{x_v}(L)R_v>1$ then
  (5.1) and (5.2) hold". Here:
  * "(5.1) There is a proper subvariety $Z$ of $X$ such that for all infinite sequences $\{x_i\}$ of distinct points
    of $X(k)\setminus Z(k)$, there is at least one $v\in S$ so that $\alpha_v(\{x_i\},L)\ge 1/R_v$."
  * "(5.2) There is a proper subvariety $Z$ … such that for any collection $\{\delta_v\}$ with each $\delta_v>0$, there
    are only finitely many solutions $y\in X(k)\setminus Z(k)$ to $d_v(x_v,y)\le H_L(y)^{-(R_v+\delta_v)}$ for all
    $v\in S$."
  * Proof: filter $H^0(mL)$ by order of vanishing at $x_v$ and apply Faltings-Wüstholz [7, Thm 9.1] (Theorem 5.4).
    $Z$ is the base locus of the maximal destabilizing subspace $W_\circ$. If $X(k)$ is not dense, take
    $Z=\overline{X(k)}$.
  * Corollary 5.2 allows equality in (5.3) when $\delta_v$ is fixed before $Z$. Corollary 5.3: if $\dim X=1$, equality
    is allowed and $Z=\varnothing$.
* **Theorem 6.1 (p. 33).** "For any ample $\mathbb Q$-bundle $L$ and any $x\in X(\bar k)$ either (a)
  $\alpha_x(L)\ge\beta_x(L)$ or (b) There exists a proper subvariety $Z\subset X$, irreducible over $\bar k$, with
  $x\in Z(\bar k)$ so that $\alpha_{x,X}(L)=\alpha_{x,Z}(L|_Z)$". Here $X$ is irreducible of dimension $n$ (§6
  preamble, p. 32).
* **Corollary 4.4 (p. 28).** "Let $X$ be an irreducible $n$-dimensional variety over Spec(k). Then for any
  $x\in X(\bar k)$ and any ample $L$ we have $\beta_x(L)\ge\frac{n}{n+1}\epsilon_x(L)$. In general, this inequality is
  strict."
* **Corollary 4.2 (p. 27).** $\beta_x(L)\ge\frac{n}{n+1}\sqrt[n]{\mathrm{Vol}(L)/\mathrm{mult}_xX}\ge\frac n{n+1}\epsilon_x(L)$.
* **Theorem 9.1 (pp. 50-51).** For $X$ irreducible of dimension $n$, $x\in X$ and $L$ ample, the following are
  equivalent: (a) $\beta=\frac n{n+1}\sqrt[n]{\mathrm{Vol}/\mathrm{mult}}$; (b) $\sqrt[n]{\mathrm{Vol}/\mathrm{mult}}=\epsilon_x$;
  (c) $\beta=\frac n{n+1}\epsilon_x$; (d) $\epsilon_x(L)=\gamma_{{\rm eff},x}(L)$.
* **Theorem 6.2, "Schmidt type" (p. 4 and p. 33).** "Let $X$ be an irreducible $n$-dimensional variety over
  Spec(k). For any ample $\mathbb Q$-bundle $L$ and any $x\in X(\bar k)$ either (a)
  $\alpha_x(L)\ge\frac n{n+1}\epsilon_x(L)$ or (b) there exists a proper subvariety $Z\subset X$, irreducible over
  $\bar k$, with $x\in Z(\bar k)$ so that $\alpha_{x,X}(L)=\alpha_{x,Z}(L|_Z)$".
  * Alternate statement (p. 4): "Then there is a proper subvariety $Z\subset X$ so that for each $\delta>0$ there are
    only finitely many solutions $y\in X(k)\setminus Z(k)$ to
    $d_v(x,y)<H_L(y)^{-(\frac{n+1}{n\epsilon_x(L)}+\delta)}$."
* **Theorem 6.3, "Roth type" (p. 4; full form p. 33).** "For all varieties $X$ over Spec(k) (possibly reducible),
  all $x\in X(\bar k)$ and all ample line bundles $L$, $\alpha_x(L)\ge\frac12\epsilon_x(L)$". Full form: "with
  equality if and only if both $\alpha$ and $\epsilon$ are computed on a $k$-rational curve $C$ such that (1) $C$ is
  unibranch at $x$, (2) $\kappa(x)\ne k$, (3) $\kappa(x)\subset k_v$, and (4) $\epsilon_{x,C}(L|_C)=\epsilon_{x,X}(L)$."
  * **Corollary 6.4:** for any $\delta>0$ only finitely many $y\in X(k)$ satisfy
    $d_v(x,y)<H_L(y)^{-(2/\epsilon_x(L)+\delta)}$.
  * No exceptional set appears.
  * **Point type.** The theorem holds for all algebraic $x$, rational ones included. By condition (2), equality is
    impossible at a $k$-rational $x$.
* **Corollary 6.5.** If $\alpha_x(L)<\frac m{m+1}\epsilon_x(L)$ for some $m<n$, then $\alpha_x(L)$ is computed on a
  subvariety of dimension $<m$.
* **Corollary 6.6.** "If $x\in X(\bar k)$, and no rational curve passes through $x$ then
  $\alpha_x(L)\ge\frac23\epsilon_x(L)$."
* **Simultaneous approximation (§7).**
  * Prop 7.3 gives equivalent formulations.
  * Thm 7.4: the conditions hold for some proper $Z$ if (a) $R_v>1/\beta_{x_v}(L)$, or (b)
    $R_v>\frac{n+1}{n\epsilon_{x_v}(L)}$.
  * Thm 7.5: they hold with $R_v=2/\epsilon_{x_v}(L)$ and $Z=\varnothing$.
  * Cor 7.6: $\sum_{v\in S}\epsilon_{x_v}(L)/\alpha_v(\{x_i\},L)\le2$.
  * Thm 7.7 characterizes equality.
  * Thm 7.8 and Cor 7.9 give a simultaneous weighted Dirichlet theorem on $\mathbb P^1$.
  * Cor 7.10, 7.11 and 7.12 are further variants.
* **Unramified covers (§8).** By Cor 8.9, "Theorems 5.1, 6.1, 6.2, 6.3, 7.4, 7.5, and 7.7, Corollaries 5.2, 5.3, 6.4,
  6.5, 6.6, 7.6, 7.10, 7.11, and 7.12 hold with $\beta$ and $\epsilon$ replaced by $\hat\beta$ and $\hat\epsilon$."
* **Theorem 10.1 (p. 53), the Vojta special case.** "Let $k$ be a number field, $X$ an irreducible $n$-dimensional
  variety over Spec(k) such that $-K_X$ is ample, $D$ a finite subset of $X(\bar k)$, and $S$ a finite set of places
  of $k$. If $\epsilon_x(-K_X)>\frac{n+1}n$ for every $x\in D$, then Vojta's Main Conjecture is true for $X$ and $D$.
  Specifically, for every $\delta>0$ and any big divisor $A$, there is a closed subset $Z\subset X$ such that for all
  $k$-rational points $P\in X(k)\setminus Z(k)$, we have
  $\sum_{v\in S,x\in D}-\log d_v(x,P)+h_{K_X}(P)<\delta h_A(P)+O(1)$."
  * Remark (p. 54): for $k$-points of $D$, $\epsilon_x(-K_X)\ge1$ suffices, by Liouville $\alpha_x(-K_X)\ge\epsilon_x(-K_X)$.

### 3.3 Basic computations used below

* Lemma 2.13: $\alpha_{x,\mathbb P^n}(\mathcal O(1))=1$ for $x\in\mathbb P^n(k)$.
* Prop 2.14(d): if $x\in X(k)$ and $L$ is very ample, then $\alpha\ge1$.
* Lemma 2.15: on $\mathbb P^1$, $\alpha_x(\mathcal O(d))=d$ if $\kappa(x)=k$, $d/2$ if $k\ne\kappa(x)\subseteq k_v$, and
  $\infty$ if $\kappa(x)\not\subseteq k_v$.
* **Theorem 2.16 (p. 19).** For a $k$-rational curve $C$ with normalization $\varphi$:
  $\alpha_{x,C}(L)=\min_{q\in\varphi^{-1}(x)}d/(r_qm_q)$, with $r_q=1$ if $\kappa(q)=k$. Also stated as [MR-L] Thm 2.8.
* Examples: $\mathbb P^1\times\mathbb P^1$ with $\mathcal O(a,b)$ has $\alpha_x=\epsilon_x=\min(a,b)$ at $k$-points
  (pp. 20, 22), and $\beta=(d_1+d_2)/2$ (p. 27). For abelian varieties $\alpha=\infty$.

### 3.4 Cubic-surface application ([MR-L] §4)

* **Theorem 4.5 (pp. 9-10).** "Let $X$ be a smooth cubic surface in $\mathbb P^3$ defined over $k$, and isomorphic over
  $k$ to the blowup of $\mathbb P^2$ at six $k$-rational points in general position. Let $x\in X(k)$ be any
  $k$-rational point, and let $C_x$ be the curve of intersection of $X$ with the tangent plane to $X$ at $x$." Then:
  * $\epsilon_x(-K)=1$ if $x$ lies on one of the 27 lines, and $3/2$ otherwise.
  * $\alpha_x(-K)=1$ on a line.
  * $\alpha_x(-K)=3/2$ if $x$ is off the lines and $C_x$ is cuspidal, or nodal with tangent slopes in $k_v$ but not in
    $k$.
  * $\alpha_x(-K)=2$ otherwise.
  * Proof mechanism: Thm 3.3(b) with $\gamma=2$ ($\pi^*L-2E$ has asymptotic base locus $\widetilde C_x$), plus
    Thm 2.8 on $C_x$ and conics.
* **Theorem 4.6 (p. 11).** If $x$ is not on a $(-1)$-curve and $D$ is nef, then
  $\epsilon_x(D)=\min\{D.L_i,D.L_{ij},D.B_i,(D.h)/2\}$.
* **Theorem 4.7 (p. 11).** Same hypotheses. If $C_x$ is cuspidal, or nodal with tangent lines defined over $k_v$ but
  not over $k$, then $\alpha_x(D)=\epsilon_x(D)$. Otherwise $\alpha_x(D)=\min\{D.L_i,D.L_{ij},D.B_i\}$.
* p. 12: "Conjecture 4.2 holds for every $x\in X$ not on a $(-1)$-curve."
* **Conjecture 4.2 (p. 9), McKinnon's conjecture as stated here.** "Let $X$ be an algebraic variety defined over
  $k$, and $L$ any ample divisor on $X$. Let $x$ be any $k$-rational point on $X$ and assume that there is a rational
  curve defined over $k$ passing through $x$. Then there exists a curve $C$ (necessarily rational) of best
  approximation to $x$ on $X$ with respect to $L$." The paper attributes it to [9, Conjecture 2.7] (McKinnon, JAG 2007).
  * [MR-L] p. 8 adds: in [9, §4], Vojta implies that a finite $\alpha_x(L)$ "is computed on a subvariety $V\subseteq X$
    of negative Kodaira dimension".
* The cubic surface is **split** throughout: an isomorphism over $k$ with the blow-up at six $k$-points, so all 27
  lines are rational.

### 3.5 Toric varieties and tori

Neither paper contains a theorem about toric varieties, algebraic tori, or anisotropic or non-split tori. The only
toric examples are $\mathbb P^n$, $\mathbb P^1\times\mathbb P^1$ and products of projective spaces. The only non-split
example is the quadric in §3.1. [MR-I] p. 54 notes $\epsilon_x(-K_X)\ge2$ on $G/P$.

## 4. Mechanism

* **Liouville ([MR-L] Lemma 3.2 and Thm 3.3).** If $L_\gamma=\pi^*L-\gamma E$ is effective, then $H_{L_\gamma}$ is
  bounded below off its base locus. Since $H_E=\prod_we_w^{-1}\ge e_v^{-1}$, this gives
  $H_L(y)\,d_v(x,y)^\gamma\ge c$ for $y\notin B$.
  * This is exactly a fixed-section, forced-divisibility argument.
  * The usable $\gamma$ ranges up to $\gamma_{\rm eff}=\tau_x(L)$, and never beyond it, because $L_\gamma$ must be
    effective.
  * Over a non-closed field, the only extra freedom is replacing $B$ by $\overline{B(k)}$ (quadric example). This
    changes the exceptional set, not the exponent bound $\gamma\le\gamma_{\rm eff}$.
* **Roth ([MR-I] Thm 5.1).**
  * Filter $H^0(mL)$ by $\mathrm{ord}_{x}$ at levels $m\gamma_j$.
  * Faltings-Wüstholz then gives finiteness off the base locus of the destabilizing subspace, provided the slope
    $\sum_jc_j(\dim V^j-\dim V^{j+1})/\dim V>1$ with $c_j=R\gamma_j$.
  * Asymptotically the slope is $R\beta_x(L)$.
  * Approximating points must be $k$-rational. The target may be algebraic, and Galois conjugates play no role in the
    single-place bound.
  * The exponent obtained is $1/R<\beta_x(L)$, the *expected* order of vanishing.
  * Because $f\le1$, $f$ is continuous and $f(\gamma_{\rm eff})=0$, we get $\beta_x(L)<\gamma_{\rm eff}$ strictly
    **[inference, immediate from the definitions on p. 25]**.
  * For an algebraic irrational $x$, Liouville loses the factor $1/d$ while Roth does not. For a $k$-rational $x$
    ($d=1$), Liouville with base locus gives $\gamma_{\rm eff}$ and Roth gives $\beta<\gamma_{\rm eff}$.
* **Seshadri** enters only through $\beta\ge\frac n{n+1}\epsilon$ (Cor 4.4) and through the induction on dimension,
  using $\epsilon_{x,Z}(L|_Z)\ge\epsilon_{x,X}(L)$ (Prop 3.4(c)).

## 5. Relevance to Question B

Setting of the record: $X=\overline L{}^{\rm tw}_5$ (weak Fano), $X'$ (Gorenstein toric Fano, $-K_{X'}$ ample),
$k=\mathbb Q$, $v=\infty$, target $e\in T_5(\mathbb Q)$ ($d=1$), $L=-K$, $n=4$, $\mathrm{Vol}=70$,
$\epsilon_e=2$, $\tau_e=\gamma_{\rm eff}=10/3$. A needed bound is $\alpha^{\rm ess}_e(-K)>15/4$.

**Applicability [inference].**
* MR require $L$ ample. Work on $X'$, where $-K_{X'}$ is ample. The point $e$ lies in the open torus, where
  $\varphi:X\to X'$ is an isomorphism; we have $X(\mathbb Q)=X'(\mathbb Q)$, and heights agree.
* $\epsilon,\beta,\gamma_{\rm eff}$ at $e$ can be computed on either model: $\mathrm{Bl}_eX\to\mathrm{Bl}_eX'$ is small
  and crepant for $\pi^*(-K)$, and volumes are birational invariants.
* All MR theorems apply. None needs splitness.

**What each theorem gives at $e$ [inference: substitution of the record's values into the quoted statements].**

| theorem | type of bound | value at $(X',-K,e)$ | compared with the needed $>15/4$ |
|---|---|---|---|
| Liouville, [MR-I] p. 24 / [MR-L] Cor 3.4 | $\alpha_e$ (all sequences) | $\ge\epsilon_e=2$ (sharp: diagonal conics, Thm 2.16 gives $\alpha=2$) | no |
| Roth type, Thm 6.3 / Cor 6.4 | $\alpha_e$ | $\ge1$ (strictly weaker than Liouville at a rational point) | no |
| Schmidt type, Thm 6.2 | $\alpha$ off $Z$ | $\ge\frac45\cdot2=8/5$ | no |
| Thm 5.1 / 6.1 ($\beta$) | $\alpha$ off $Z$, i.e. $\alpha^{\rm ess}\ge\beta_e$ | $\beta_e(-K)\in[2.31,\,2.61]$ (see below) | no |
| Cor 8.9 ($\hat\beta,\hat\epsilon$) | same with covers | no gain: $X$ (a complete smooth toric variety) is simply connected, so there are no nontrivial geometric étale covers | no |
| [MR-L] Thm 3.3(a), $\gamma\uparrow\gamma_{\rm eff}$ | $\alpha$ off $B=\pi(\mathbf B(L_\gamma))$, i.e. $\alpha^{\rm ess}\ge\tau_e$ | $10/3$ (attained, because $L_{10/3}$ is effective via the Vandermonde in $H^0(3(-K))$ of order 10) | no: Cilleruelo-Córdoba |
| Thm 3.3 applied to other toric $L$, restricted to balanced points | via record Prop 1.3 | $\le15/4$ (record, not MR) | no |
| Thm 10.1 (Vojta with $D=\{e\}$ on $X$, not on $\mathrm{Bl}_eX$) | off $Z$ | $d(e,P)\ge H_{-K}(P)^{-(1+\delta)}$, i.e. exponent 1 | no; this is **not** the record's $V_E(5)$ |

**$\beta_e(-K)$ bounds [inference, computed in this check].**
* Lower bound, from Cor 4.2: $\beta\ge\frac45\,70^{1/4}=2.314$.
* Upper bound:
  * $f(\gamma)=1-\gamma^4/70$ on $[0,2]$ (nef range, [MR-I] Remark p. 27).
  * $\mathrm{Vol}^{1/4}$ is concave and $C^1$ on the big cone (Lazarsfeld 11.4.9; Boucksom-Favre-Jonsson), so on
    $[2,10/3]$, $f^{1/4}$ lies below its tangent line at $2$.
  * $f=0$ beyond $10/3$.
  * Hence $\beta\le1.9086+0.6949=2.603$. The crude bound $f\le f(2)$ gives $2.937$.
* So the Roth/Faltings-Wüstholz constant at $e$ lies in $[2.31,2.61]$, well below $\tau_e=10/3$ and $15/4$.

**Answer.**
* **Nothing in either paper gives $\alpha^{\rm ess}_e(-K)>15/4$, or $A_5(\delta)$ for any $\delta<1/30$.**
* The best unconditional output at the rational point $e$ is [MR-L] Thm 3.3 at $\gamma=\gamma_{\rm eff}=10/3$,
  i.e. Cilleruelo-Córdoba ($\delta=1/10$). On balanced points the record's 15/4 is the ceiling of the same mechanism.
* No theorem in either paper bounds approximation to a rational point from below beyond $\gamma_{\rm eff}$ on a
  Zariski-dense set. The structure of the Liouville proof (effective $L_\gamma$ needed) and of the Roth proof
  ($\beta<\gamma_{\rm eff}$) excludes this.
* The only twist-sensitive device in either paper is "replace $B$ by $\overline{B(k)}$" (non-split quadric). It
  acts on the exceptional set, not on the exponent. It could at most enlarge the set of rational points where the
  $10/3$ bound holds.

## 6. Corrections and confirmations for the record

1. **Confirmed: "$\alpha_x(L)\ge\epsilon_x(L)/2$ (Thm 6.3)".**
   * This is [MR-I] Theorem 6.3 (p. 4, full form p. 33). It holds for all $x\in X(\bar k)$, possibly reducible $X$,
     and ample $L$.
   * Equality holds iff there is a $k$-rational curve $C$, unibranch at $x$, with $\kappa(x)\ne k$,
     $\kappa(x)\subset k_v$ and $\epsilon_{x,C}=\epsilon_{x,X}$.
   * round22 README §1.1 lists this citation as "unconfirmed"; it is now confirmed. The $P^n$ constant $n/(n+1)$ is
     Thm 6.2, also confirmed.
2. **Correction to lit.md §3: "No exceptional set appears in the inequality itself (it bounds the infimum over all
   sequences)".** This is true only of Thm 6.3 and Cor 6.4. The central results are essential-type statements with an
   exceptional proper subvariety $Z$, valid for every $x\in X(\bar k)$ including rational $x$:
   * Thm 5.1, where $Z$ is the base locus of the Faltings-Wüstholz destabilizing subspace and depends on $R$;
   * Thm 6.1 and 6.2, where the alternative (b) is that "$\alpha$ is computed on a proper subvariety";
   * Thm 7.4.

   Consequence [inference]: $\alpha^{\rm ess}_x(L)\ge\beta_x(L)$. The record's conclusion is unchanged, because
   $\beta_e(-K)\le2.61<10/3$.
3. **Confirmed, and sharpened: "$\beta\le\tau$".**
   * MR's $\beta_x(L)=\int_0^{\gamma_{\rm eff}}\mathrm{Vol}(L_\gamma)/\mathrm{Vol}(L)\,d\gamma$ ([MR-I] p. 25) is the
     record's expectation constant.
   * MR's $\gamma_{\rm eff}$ is the record's $\tau$.
   * $\beta<\gamma_{\rm eff}$ strictly, since $f\le1$, $f$ is continuous and $f(\gamma_{\rm eff})=0$.
   * At $e$: $\beta_e(-K)\in[2.31,2.61]$, numerically below $10/3$ [inference; bounds above].
4. **Confirmed, with nuance: "Roth-type weaker than Liouville at a rational target".**
   * [MR-L] p. 1: "On $\mathbb P^1$, except for the case that $x\in\mathbb P^1$ is a rational point of the base number
     field, Liouville's theorem is weaker than Roth's."
   * [MR-I] p. 24: $\alpha_x(L)\ge\epsilon_x(L)$ for $x\in X(k)$ is "a stronger inequality than
     $\alpha_x(L)\ge\frac12\epsilon_x(L)$".
   * For the essential versions the same holds: Roth gives $\beta$ and Liouville with base locus gives
     $\gamma_{\rm eff}>\beta$ [inference].
5. **Confirmed: "McKinnon's conjecture concerns $\alpha_x$".**
   * [MR-L] Conjecture 4.2 (p. 9), quoted in §3.4, concerns curves of best approximation.
   * Def 4.1 defines these through sequences with approximation constant equal to $\alpha_x(L)$, the infimum.
   * The conjecture requires $x$ $k$-rational, $L$ ample, and a $k$-rational curve through $x$.
   * Neither paper defines or uses an essential constant.
6. **Correction or precision on the cubic-surface entry (lit.md §2 "cubic surfaces, points off lines (McKinnon-Roth
   2016)" and §4 table).**
   * Thm 4.5 assumes $X$ is "isomorphic over $k$ to the blowup of $\mathbb P^2$ at six $k$-rational points in general
     position". So the result is for **split** cubic surfaces only.
   * Thm 4.5 covers all $k$-points for $-K$, including points on lines.
   * Thms 4.6 and 4.7 cover nef $D$ at points not on a $(-1)$-curve.
   * The values are $\epsilon_x(-K)\in\{1,3/2\}$ and $\alpha_x(-K)\in\{1,3/2,2\}$. So $\alpha_x>\epsilon_x$ can occur at
     a rational point, through Thm 3.3 with $\gamma=2>\epsilon$. This is still $\le\gamma_{\rm eff}$.
   * The record's "mechanism = Liouville + lines" is accurate. More precisely: Liouville with base locus (Thm 3.3(b),
     $\gamma=2$, base locus $\widetilde C_x$) plus Thm 2.8 on the tangent cubic and on conics through a line.
7. **Precision on lit.md §3, "Hypotheses: … a single place (archimedean distance), so $S=\{\infty\}$ is fixed".**
   * [MR-I] allows any single place, archimedean or not. §§5 and 7 give multi-place (simultaneous) versions for
     finite $S$.
   * Neither changes anything for $e$: A-points carry no $p$-adic proximity [inference].
8. **Precision on the exceptional set of the Liouville extension.**
   * [MR-L] uses the "asymptotic base locus $B'$ of $L_\gamma$", with $B=\pi(B')$. Its Remark allows $B$ to be
     replaced by $\overline{B(k)}$.
   * It requires $L$ ample, $\gamma$ rational, and $L_\gamma$ "in the effective cone".
   * For $x$ of degree $d$ the bound is $\gamma/d$.
   * lit.md's "stable base locus of $\pi^*L-\gamma E$ … fixed closed set depending only on $(X,L,x,\gamma)$" is
     consistent with this.
9. **New caution: [MR-I] Theorem 10.1 ("special case of Vojta's main conjecture") is not the record's $V_E(5)$.**
   * It is Vojta for $X$ itself with $D$ a finite set of points, $-K_X$ ample and $\epsilon_x(-K_X)>(n+1)/n$. It has
     no $(n-1)E$ term from the blow-up.
   * At $e$ on $X'$ its conclusion amounts to exponent $1$ for $-K$.
   * It must not be cited as a proved case of $V_E(5)$ [inference from the displayed inequality].
10. **Bibliographic caveat.** The journal data (Invent. Math. 200 (2015) 513-583; Eur. J. Math. 2 (2016)) are not
    printed in these arXiv PDFs. The EJM file is arXiv v1 (June 2013), so its numbering (Thm 3.3, Cor 3.4,
    Conj 4.2, Thms 4.5-4.7) may differ from the published version.
11. **Dictionary check.** MR's normalisation matches the record's $\gamma=4/(1+2\delta)$: $\alpha$ is the infimum of
    $\gamma$ with $d_v^\gamma H_L$ bounded, heights are multiplicative relative to $k=\mathbb Q$, and $d_v$ at the real
    place is the chordal/Fubini-Study distance, $\asymp$ the local coordinate distance (Lemma 2.5). With $H=R''^2$ and
    $d\asymp\mathrm{diam}/R''$, an A-point is $d\le H^{-(1+2\delta)/4}$. **Confirmed.**
12. **Confirmed: the record's "$\epsilon_e=\alpha_e=2$".** It is consistent with MR's Liouville ($\alpha\ge\epsilon$ at
    rational points) and with Thm 2.16 ($\alpha=\deg=2$ on a smooth $\mathbb Q$-rational conic through $e$).
