## 5. Calibration: the function-field analogue of Question B, and (R$_4$)

Work over $\mathbb C$. Over $\mathbb C(t)$ the twist is invisible, since $X_{\mathbb C}=\overline{LM}_5$ is split.
Let $Y=\mathrm{Bl}_eX$. For a non-constant $f:\mathbb P^1\to X$ through $e$, put $h(f)=\deg f^*(-K_X)$,
$\mu(f)=$ the contact order at a point over $e$, and $\alpha(f)=h(f)/\mu(f)$.

> **FF-$A_5(\eta)$.** There is a proper Zariski-closed $Z\subset T_5$ such that every rational curve $f$ through $e$
> with $\alpha(f)\le4-\eta$ has $f(\mathbb P^1)\subset Z$.

By §4.1 this says: *primitive 5-point $\mathbb C(t)$-clusters at negative level
($\max\deg(z_i-z_j)\le D(2-\eta)/(4-\eta)$) have their similarity class in a fixed proper closed set.* It is Question
B with "A-point" replaced by "rational curve with $\alpha<4$". Over $\mathbb Q(i)(t)$, curves with $\alpha<4$
specialise to 5-point integer clusters with $C\to0$ (Theorem 1.2(b)). The record's only 5-point $C\to0$ source is the
cyclotomic family, which is not of this type.

**Proposition 5.1 (proved, given the conjecture named).** Fix a smooth projective toric $X_r\supset\mathbb G_m^r$
whose fan contains all rays in $\{0,\pm1\}^r$. Let $\mathcal D=\pi^*\partial X_r+E$ on $Y_r=\mathrm{Bl}_eX_r$. This is
a simple normal crossings divisor, since $e\notin\partial X_r$, and $K_{Y_r}+\mathcal D=\pi^*(K_{X_r}+\partial X_r)+(r-1)E+E=rE$.
The truncated function-field Vojta conjecture for $(Y_r,\mathcal D)$ reads
$\deg f^*(K+\mathcal D)\le-2+N^{(1)}_{\mathcal D}(f)+\varepsilon h_A(f)$ off $Z_\varepsilon$. For a curve with one
branch through $e$ of multiplicity $\mu$, meeting $\partial X_r$ at $|S|$ points, it becomes
$$r\mu\le|S|-1+\varepsilon h_A(f).\tag{5.1}$$
* $r=4$ (Question B, on $X$): $4\mu\le h-1+\varepsilon h$. So off $Z_\varepsilon$, $\alpha\ge4/(1+\varepsilon)$, which
  is FF-$A_5(\eta)$ for $\varepsilon\asymp\eta$. Dropping the truncation gives the analogue of $V_E(5)$ with the
  same conclusion.
* A rank-$r$ cluster at level $c$ gives such an $f$ on $X_r$: take a basis of difference vectors, so that the
  cocharacters lie in $\{0,\pm1\}^r$. Then $\mu=L+1$ and $|S|\le2D$. At $c=0$, (5.1) allows $r=4$ only when
  $\varepsilon\ge1/(2D)$; the cube at $D=8$ has $16\le15+16\varepsilon$, so it lies in $Z_\varepsilon$ for $\varepsilon<1/16$.
  It **forbids $r\ge5$ off $Z_\varepsilon$ for every $\varepsilon<1/4$**, since $5D/2\le2D-1+2D\varepsilon$ fails.
* With $\varepsilon=0$, (5.1) is exactly Theorem D ($r(L+1)+1\le2D$ for dominating families). Theorem D is the case
  of (5.1) that is a theorem: free curves.

**Correction to (S4).** (S4) used $\mathcal D=\pi^*\partial X$ without $E$. That gives $(r-1)\mu\le|S|-2+\varepsilon h$,
which is weaker and allows $r=5$. Vojta's conjecture applies to every snc divisor, and adding $E$ is legitimate and
standard: it is the divisor of $V_E$, and the one behind Question B. **So (R$_4$) has exactly the Vojta shape; it is
not sharper, and $r=5$ is not borderline.** The cube's rank 4 is the borderline case. What Vojta does *not* give is
(R$_4$) itself: it puts rank-$\ge5$ cluster curves into an exceptional set $Z_\varepsilon\subset\mathbb G_m^5$, and
excluding curves inside $Z_\varepsilon$ is the excess-intersection content of (R$_4$).

*Known theorems are too weak.* The function-field gcd bound for $\gcd(u-1,v-1)$, from Corvaja–Zannier (J. Alg. Geom.
2008, arXiv math/0512074) and its positive-characteristic version (JEMS 2013), has the shape
$\sum_\nu\min(\nu(1-u),\nu(1-v))\le C(\deg u\cdot\deg v\cdot\chi)^{1/3}$. I found this through search results only,
and the constants are not checked. For two cluster ratios $\deg u,\deg v\asymp D$ and $\chi\asymp2D$, so the bound is
$\asymp D$, the same order as the contact $L+1\approx D/2$. It is vacuous here unless the constant is below about
$1/2^{4/3}$; the stated positive-characteristic constant $3\sqrt[3]2$ is not. This confirms the last sentence of (S4)
at the level of orders of magnitude.

**What the 5-point version sees.** A 5-subset of a rank-$r$ cluster lies in a subtorus of rank $\le\min(r,4)$, with
contact $\ge L+1$. Its $\alpha$ drops below 4 only through the gcd ($g>2c$, §4.1). So:
1. FF-$A_5$ never sees the rank of the whole cluster. A rank-5 cluster whose 5-subsets all have $g\le2c$ is invisible
   to it. In the other direction, (R$_4$) says nothing about isolated 5-point clusters at negative level.
2. Large clusters force $\alpha<4$. If a cluster at level $c$ has $M$ points, averaging over 5-subsets gives one with
   $g\ge D\,\phi(M)$, $\phi(M)=2\binom{\lceil M/2\rceil}5/\binom M5$ (Jensen on the column densities), so
   $\alpha\le4(1-\phi(M))\cdot D/(D-2c)$. Here $\phi(M)>0$ iff $M\ge10$, and $\phi\to1/16$, giving
   $\alpha\to15/4$. This is the function-field form of `gauss.md` Lemma 1.2. So $\sup_DM(D,c)\ge10$ produces
   $\alpha<4$ curves; the cube (8 points) does not.
3. Dimension counts. Rank-4 5-point clusters at level $c'$ have expected dimension $2D-2-4L$ modulo affine maps. At
   $K=(D-1)/2$ (odd $D$) this is $0$ (rigid; $\alpha=4D/(D+1)$), and for $K<(D-1)/2$ it is negative. FF-$A_5(\eta)$
   needs only $D\le4/\eta$, so it is consistent with any finite number of rigid solutions per $D$. It fails only if
   positive-dimensional families at negative level exist in infinitely many directions. Neither is known.

**Effect on confidence in Question B.**
* *If the function-field truth is FF-$A_5$ (and (R$_4$)):* weak positive evidence. It is the statement Vojta
  predicts, and its arithmetic analogue is Question B. But the function field is easier, because it has
  $\asymp D/2$ independent linear conditions (ff Lemma 1.1; residue-field obstruction), so a function-field proof
  transfers nothing (README §2.1 item 7).
* *If FF-$A_5$ fails:* say a positive-dimensional family of $\mathbb C(t)$-curves with $\alpha\le4-\eta$, dense in
  $T_5$, or rank-5 clusters for infinitely many $D$ outside every fixed $Z$. Then the function-field Vojta
  conjecture fails for $(Y,\pi^*\partial+E)$. That would be a counterexample to a widely believed conjecture, and it
  would **lower confidence in Question B substantially**, since Question B rests only on Vojta. If the
  family is defined over $\mathbb Q$, Theorem 1.2 turns it into Zariski-dense A-points, i.e. an actual refutation of
  $A_5(\delta)$ for all small $\delta$ (§6).
* The asymmetry is the point. A function-field proof is weak evidence; a function-field counterexample over
  $\mathbb Q(i)$ is arithmetic.
