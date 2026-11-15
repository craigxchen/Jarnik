# Round 29, lane ff: the Cilleruelo–Granville problem over function fields

## 0. Verdict

**The function-field Cilleruelo–Granville problem is not decided here. The function-field model reproduces the
arithmetic exceptional-set barrier; it does not bypass it.** Status: *obstruction*, with new exact results.
All statements below are proved unless labelled; "exact" means a finite exact computation.

1. **Linearisation** (Lemma 1.1). Over $\mathbb C(t)$ the arc condition is $L=D-1-\lfloor D/2+c\rfloor$ *exact real
   linear* equations on exponent vectors. Their columns $v(\rho_k)=(\mathrm{Im}\rho_k^m)_{m\le L}$ lie on a
   2-dimensional "moment surface". A cluster through $z_1$ is $V'\cap\{0,1\}^D$ for a linear subspace $V'$. The pair
   distance is $\ge L+1\approx D/2-c$ (Lemma 1.2, char-free). This corrects the coordinator's $D/4$ and the
   coordinator's count in (i), which should use the rank of the code, not $M-1$.
2. **Pair bound** (Thm A): $M\le2^{2c+2}(D+2)$, linear in $D\asymp\log R$. **Fakes** (Thm B): arbitrary vectors in
   place of the moment vectors give $M\ge D+2c+2$. This is (F1) in the function field: linear structure plus all pair
   inequalities cannot beat linear growth.
3. **Eight points for every $c$** (Thm C, exact): an explicit 2-parameter family over $\mathbb Q$ (and over
   $\mathbb F_q$, $q=19,23,31,43$ checked) of 8-point clusters at $D=8$, $c=0$. Composition extends it to all
   $D\equiv0\ (8)$. So $\sup_DM(D,c)\ge8$ over $\mathbb C(t)$, $\mathbb Q(i)(t)$ and $\mathbb F_q(t)$. By-product: the first
   polynomial template in $\mathbb Z[i][u]$ with 8 points at bounded $C$ (huge $C$; a mechanism, not a record).
4. **Dominating families** (Thm D, uses the derivation through generic smoothness): a dominating family of rank-$r$
   clusters has $r(L+1)+1\le2D$. At $c=0$ this gives $r\le3$ and $M\le6$. This is the function-field "4 false / 5
   true" for *covering* families. The cube family is rank 4 and non-dominating, i.e. an exceptional-set family.
5. **What remains.** $\sup_DM<\infty$ iff the code rank is bounded. The dimension count predicts rank $\le4$ and
   $\limsup M=8$, which Prop. 3.3 makes sharp. Proving it is the excess-intersection statement (R$_4$), the
   function-field Lang–Vojta exceptional-set problem for code subtori. Mason–Stothers, Wronskians and the truncated
   SMT are either vacuous on squarefree $N$ or exactly critical. Deformation theory sees only dominating families.
6. **Arithmetic transfer: negative.** The ultrametric step needs $\asymp\log R$ independent exact linear conditions,
   which the arithmetic $p$-adic twin (one congruence) does not provide. The derivation step needs a family, and an
   integer circle is a single point; its arithmetic counterpart is Question B itself. Fermat quotients (Buium
   $\delta$-characters) are residues mod $p^2$, already faked. No growth improvement follows. The growth rate
   $(\frac23+\varepsilon)\log R/\log\log R$ is unchanged.

## 1. Setting, normalisations, corrections to the coordinator's observations

**Setting (twisted form).** Let $k=\mathbb C$, with complex conjugation acting on coefficients of $\mathbb C[t]$
(the same text works over $\mathbb Q(i)$, and over $\mathbb F_{q^2}$, $q\equiv3\ (4)$, with Frobenius as conjugation,
wherever only coefficients of index $<p$ enter). A *cluster* is a set of distinct $z_1,\dots,z_M\in\mathbb C[t]$ of
degree $D$ with $z_j\bar z_j=N$ for all $j$ and $\deg(z_i-z_j)\le D/2+c$. Put
$$K:=\lfloor D/2+c\rfloor,\qquad L:=D-1-K\quad(=D/2-c-1\text{ for even }D,\ c\in\mathbb Z).$$
$M(D,c)$ is the largest cluster size. The question is whether $\sup_D M(D,c)<\infty$.

**Normalisations (proved, elementary).** (N1) If $D>2c+1$, all $z_j$ have the same leading coefficient (else
$\deg(z_i-z_j)=D$), so WLOG all are monic. (N2) A common factor $g$ of degree $d$ can be removed; the condition becomes
$\deg\le (D-d)/2+c-d/2$, which is stronger. So WLOG $\gcd(z_1,\dots,z_M)=1$; then $N$ has no real root, and
$N=\prod_k\big((t-\rho_k)(t-\bar\rho_k)\big)^{e_k}$ with $\rho_k$ in the upper half-plane, $\sum e_k=D$. Points are
$z_a=\prod_k(t-\rho_k)^{a_k}(t-\bar\rho_k)^{e_k-a_k}$, $a\in\prod_k[0,e_k]$. "Squarefree" means all $e_k=1$; then
points are $a\in\{0,1\}^D$.

**Lemma 1.1 (linearisation; char 0).** $\deg(z_a-z_b)\le K$ iff $\sum_k(a_k-b_k)\,\mathrm{Im}(\rho_k^m)=0$ for
$m=1,\dots,L$.
*Proof.* Two monic polynomials of degree $D$ differ in degree $\le K$ iff their elementary symmetric functions
$e_1,\dots,e_L$ agree, iff (Newton, char 0) the power sums $p_1,\dots,p_L$ of their root multisets agree. Here
$p_m(z_a)-p_m(z_b)=\sum_k(a_k-b_k)(\rho_k^m-\bar\rho_k^m)=2i\sum_k(a_k-b_k)\mathrm{Im}\rho_k^m$. ∎

Equivalently (the coordinator's (ii), which is correct): $u=z_a/z_b$ lies in the *group*
$U_{L+1}=1+t^{-L-1}\mathbb C[[t^{-1}]]$, and $\log:U_1\to t^{-1}\mathbb C[[t^{-1}]]$ is a group isomorphism with
$\log u=\sum_m\frac{t^{-m}}m\sum_k(a_k-b_k)(\bar\rho_k^m-\rho_k^m)$. The arc condition is therefore **$L$ exact
real-linear equations on the exponent vector**, where the archimedean arc condition is one inequality of precision
$R^{-1/2}$. Put $v(\rho):=(\mathrm{Im}\rho^m)_{m\le L}\in\mathbb R^L$ (the *moment surface* $\Sigma_L=v(\mathbb
C\setminus\mathbb R)$, real dimension 2, $v(\bar\rho)=-v(\rho)$).

**Base-point form.** Fix $Q=z_1=\prod_{k\le D}(t-r_k)$ (squarefree case). The cluster through $Q$ is
$\{z_F=Q\prod_{k\in F}\frac{t-\bar r_k}{t-r_k}\}$ over the *zero-sum sets* $F\subseteq[D]$,
$\sum_{k\in F}v(r_k)=0$; equivalently $\deg\mathrm{Im}P_F\le|F|-L-1$ for $P_F=\prod_{k\in F}(t-r_k)$. So the
cluster is $V'\cap\{0,1\}^D$ with $V'=\ker[v(r_1)|\cdots|v(r_D)]$ a **linear** subspace, $\dim V'\ge D-L$.

**Lemma 1.2 (pair bound, char-free).** If $z_a\ne z_b$ are in a cluster then $\|a-b\|_1\ge L+1$.
*Proof.* $z_a=gA$, $z_b=gB$ with $g=\gcd$, $A,B$ monic coprime, $\deg A=\deg B=w:=\|a-b\|_1$. Then
$\deg(z_a-z_b)=D-w+\deg(A-B)\le K$, so $\deg(A-B)\le w-L-1$. If $w\le L$ then $A=B$, hence $A=B=1$, $w=0$. ∎

**Corrections to the coordinator's observations.**
* (i) The equations run over $m=1,\dots,L=D-1-K$ (for even $D$: $m\le D/2-c-1$). The parameter count is right in form
  but is misapplied: a cluster whose difference vectors span a space of real **rank $r$** costs $rL$ equations
  (Lemma 1.1 is linear), not $(M-1)L$, and a rank-$r$ set can hold up to $2^r$ points. The count
  $rL\le 2D-2$ gives $r\le4$ for $D\ge7$ ($c=0$) and for $D>10c+6$; with Prop. 3.3 below this predicts
  $\limsup_DM(D,c)=8$, not $4$ or $5$. The "4 false / 5 conjectured" pattern does appear in function fields, but as the
  dominating-family theorem (Theorem D), not as a parameter count.
* (ii) Correct; Lemma 1.1 is its quantitative form.
* (iii) Corrected: the minimum distance is $L+1=D-K\approx D/2-c$ (Lemma 1.2), not $\approx D/4$. The Vandermonde bound
  (support $s$ with $2s\ge L+2$) is weaker. After a base point the cluster is a linear section $V'\cap\{0,1\}^D$, not
  only a coset.
* (iv) Verified, with a proof in our setting (Theorem D); it uses generic smoothness, i.e. characteristic 0.

## 2. What the record already has

* `round18/theory.md` §2.3: in the function-field setting no covering family of genus-$\le1$ curves on
  $\overline M_{0,8}$ has balanced class of small imbalance; for genus $g$ through a general point the deformation
  argument gives $K\cdot\gamma\le2g-2$. That is a deformation (free-curve) bound for **covering** families, like
  Theorem D below, on a different variety.
* `threshold/L8.md` §7: the function-field analogue of $(L_8)$ asks for the absence of balanced curves of genus $\le1$;
  their expected dimension is negative but this is "not proved". The truncated second main theorem (braid
  arrangement in $\mathbb P^5$) gives no contradiction ($465\le621$).
* `growth/families.md` §4 (Lemma P): polynomial templates $z_j\in\mathbb Z[i][u]$ are exactly function-field
  clusters over $\mathbb Q(i)$ evaluated at $u\in\mathbb Z$. The exhaustive search with small Gaussian-integer roots
  found at most 6 points at bounded $C$ ($D=4$), 5 ($D=6$), 4 ($D=8$).
* Not in the record: the linearisation, the reduction to codes of bounded rank, explicit 8-point function-field
  clusters, the function-field fakes, and the dominating-family theorem.

## 3. Small cases

**Proposition 3.1 ($L\le1$; proved).** If $L\le0$, $M(D,c)=2^D$ (generic squarefree $N$; every transversal is in the
cluster). If $L=1$, $M(D,c)=\binom{D}{\lfloor D/2\rfloor}$.
*Proof for $L=1$.* By Lemma 1.1 the cluster is a fibre $\{a:\sum_ka_ky_k=\text{const}\}$ of $a\in\prod[0,e_k]$ with
$y_k=\mathrm{Im}\rho_k>0$. A fibre is an antichain in the product of chains $\prod[0,e_k]$. Its largest antichain is
its largest rank level (de Bruijn–Tengbergen–Kruyswijk; memory-quote). The block-summing map
$\{0,1\}^D\to\prod[0,e_k]$ is onto and rank-preserving, so that level has at most $\binom D{\lfloor D/2\rfloor}$ elements
(Sperner). Equality: $D$ distinct roots with equal imaginary parts. ∎

**Proposition 3.2 ($L=2$ is a free-column problem; proved).** $v(x+iy)=(y,2xy)$ maps the upper half-plane bijectively
onto $\{(a,b):a>0\}$. So for squarefree $N$ with $L=2$, $M$ is the largest $|(s_1+U)\cap\{\pm1\}^D|$ over subspaces
$U$ with $U\cap\mathbb R^D_{\ge0}=0$ (Gordan: some row of the $2\times D$ matrix is positive) and
$e_k-e_l\notin U$ for $k\ne l$ (distinct columns). `ff_checks/fakeexact.py` enumerates all such closures
(exhaustive breadth-first search; the LP feasibility tests are in floating point). Results: $M(5,0)=3$, $M(6,0)=5$
(after examining 245,105 closures). As a validation, with the separation test switched off, the same code reproduces the Sperner values $3,6,10$ at $L=1$ (`fakeexact_val.py`). The value $5$ agrees with `growth/families.md` §4, which found 5 at $D=6$ with
Gaussian-integer roots. My integer-box searches with entries $\le3$ found only $4$, so small boxes miss the optimum.

**Proposition 3.3 (codes of bounded real rank; computation).** Let $C=V\cap\{0,1\}^D$, where $V$ is a linear
subspace of dimension $\le r$, $0\in C$, and all distances are $\ge\theta D$. Project onto an information set: $C$
corresponds to $S\subseteq\{0,1\}^r$ spanning $\mathbb R^r$, and every coordinate is a linear functional with
$\ell(S)\subseteq\{0,1\}$. Maximising the minimum distance over the weights of the admissible functionals is an LP.
Exhaustive over all $S$ (`ff_checks/rank4codes.py`, `rankcodes_thr.py`):
* $\theta=\frac12$: $|C|\le6$ for $r=3$ and $|C|\le8$ for $r=4$. Both are attained, the second by
  $\{0,\mathbf 1,\mathbf 1_{H_i},\mathbf 1_{H_i^c}\}$ with $H_1,H_2,H_3$ the coordinate half-cubes of $\mathbb F_2^3$.
* $\theta=\frac47$ ($D=7$, $c=0$): $|C|\le4$ for $r=3$ and $|C|\le5$ for $r=4$.

**Table (squarefree $N$).**

| $D$ | $c=0$: $L$ | $M(D,0)$ | $c=1$: $L$ | $M(D,1)$ |
|---|---|---|---|---|
| 2 | 0 | 4 | 0 | 4 |
| 3 | 1 | 3 | 0 | 8 |
| 4 | 1 | 6 | 0 | 16 |
| 5 | 2 | 3 (exact search) | 1 | 10 |
| 6 | 2 | 5 (exact search) | 1 | 20 |
| 7 | 3 | $\ge4$ (numerical), $\le A(7,4)=8$ | 2 | $\ge6$ (exact search, see §6) |
| 8 | 3 | $8\le M\le A(8,4)=16$ (8 exact, Thm C) | 2 | $\ge8$ |

At $D=7$, $c=0$, rank-3 codes are realised numerically (`d7b.py`). No rank-4 code was realised from 20 random starts
per code class; this is evidence that $M(7,0)=4$, not a proof. At $D=8$, 175 rank-5 spans inside the extended Hamming code were
generated and 40 tested with 40 starts each; none was realised (`rank5.py`; evidence). The cyclotomic $N=t^{2D}+1$
gives only $M=4$ for even $D$ and $2$ for odd $D$, for $D\le20$ (`brute.py`).

## 4. Main results

**Theorem A (pair bound; proved, any field, twisted or split).** $M(D,c)\le A_2(D,L+1)\le2^{2c+2}(D+2)$, where
$A_2(n,d)$ is the binary code bound.
*Proof.* Encode $a\in\prod[0,e_k]$ in unary (block $k$: $a_k$ ones, then $e_k-a_k$ zeros); $\ell^1$ distance becomes
Hamming distance, so Lemma 1.2 makes the cluster a binary code of length $D$ and distance $d=L+1=D-K$. If $n=D\ge2d$,
shortening gives $A_2(n,d)\le2^{n-2d}A_2(2d,d)$, with $n-2d=2K-D\le2c$, and Plotkin gives $A_2(2d,d)\le4(d+1)$.
If $D<2d$, Plotkin gives $A_2\le2d$. ∎

This is the function-field Cilleruelo–Córdoba bound: linear in $D\asymp\log R$. Over $\mathbb C$ there are no small
primes, so the pigeonhole step behind $\log R/\log\log R$ has no analogue. Over $\mathbb F_q$ it does (degree-1 places).

**Theorem B (function-field fakes; proved).** Let $D=2^k$, $0\le c\le D/2-k-1$, and $L=D/2-c-1$. There is an
$L\times D$ real matrix $W$ with distinct columns, all with positive first entry after a fixed sign normalisation, such
that $\{s\in\{\pm1\}^D:Ws=\text{const}\}$ has at least $D+2c+2$ elements, with pairwise distances $\ge D/2\ge L+1$.
*Proof.* Index columns by $u\in\mathbb F_2^k$. Choose $B\subseteq\mathbb F_2^k\setminus0$ spanning $\mathbb F_2^k$ with
$|B|=L$, and put $A=\mathbb F_2^k\setminus(B\cup0)$, so $|A|=D/2+c$. Let $W=(\chi_b(u))_{b\in B,u}$. By orthogonality
$\ker W\supseteq\mathrm{span}\{\mathbf 1,\chi_a:a\in A\}$, which contains the $0/1$ vectors
$0,\mathbf 1,(\mathbf 1\pm\chi_a)/2$: $2|A|+2$ of them, at pairwise distance $D/2$ or $D$. The columns are distinct
because $B$ spans. A generic $g\in\mathrm{rowspace}(W)$ has no zero entry ($\chi_b(u)=\pm1$); flipping the sign
of the columns with $g_u<0$ (choosing $\bar\rho$ instead of $\rho$) and using $g$ as the first row gives the
normalisation. ∎

So in the function field a fake is a set of arbitrary vectors in $\mathbb R^L$ in place of the moment vectors
$v(\rho_k)\in\Sigma_L$. Fakes satisfy Lemma 1.1's exact linear and group structure and every pair inequality
(Lemma 1.2), and they have linear growth. Everything that beats them must use the fact that the columns lie on the
2-dimensional surface $\Sigma_L$. This is the exact function-field analogue of (F1).

**Theorem C (8-point clusters for every $c$; proved, exact).** Let $F$ be $\mathbb Q$, or $\mathbb F_q$ with
$q\equiv3\ (4)$, $p\ge5$. Take $y_0,\dots,y_3\in F^*$ with $\sum y_k=0$ and $\prod y_k=1$. Put
$\varepsilon=(1,-1,-1,1)$, $r_k=\varepsilon_k/y_k+iy_k$, and $Q=\prod_{k}(t^2-r_k^2)$ (degree 8). Assume the 16 numbers
$\pm r_k,\pm\bar r_k$ are distinct. Then the eight sets
$$\emptyset,\ \text{all},\ \{\pm r_0,\pm r_1\},\{\pm r_2,\pm r_3\},\{\pm r_0,\pm r_2\},\{\pm r_1,\pm r_3\},\
\{r_0,r_1,r_2,r_3\},\{-r_0,\dots,-r_3\}$$
are zero-sum for $L=3$. Hence $\{z_F\}$ is an 8-point cluster with $D=8$, $c=0$ and squarefree $N=Q\bar Q$ over
$F(i)[t]$. Composing with any $h\in F[t]$ of degree $q'$ gives 8-point clusters with $D=8q'$ and $c=0$, hence for every
$c\ge0$. So $\sup_DM(D,c)\ge8$ for all $c$, over $\mathbb C(t)$, $\mathbb Q(i)(t)$ and $\mathbb F_{q^2}/\mathbb F_q(t)$.
*Proof.* Here $\mathrm{Im}\,r_k^2=2\varepsilon_k$ and $\mathrm{Im}\,r_k^3=3/y_k-y_k^3$ (using $x_k^2y_k^2=1$). Sets of the form
$\{\pm r_a,\pm r_b\}$ have vanishing odd moments, and $m=2$ gives $2(\varepsilon_a+\varepsilon_b)=0$ for the four
listed pairs. For $\{r_0,\dots,r_3\}$: $m=1$ gives $\sum y=0$, $m=2$ gives $\sum2\varepsilon=0$, and $m=3$ gives
$\sum(3/y_k-y_k^3)=0$. The last follows from the identity (sympy, `cube_rational.py`) that for $d=-(a+b+c)$,
$abcd\sum(3/y-y^3)=-3(a+b)(a+c)(b+c)(1-abcd)$. Its negative is zero-sum because $p_m(-S)=(-1)^mp_m(S)$, and the full set
is zero-sum as a union. Over $\mathbb F_q$, Newton's identities are used only for $m\le3<p$. Composition multiplies
all degrees by $q'$. ∎

*Exact checks* (`cube_exact.py`, `cube_Fq.py`). For $y=(-\frac94,-\frac23,\frac14,\frac83)$ and two further rational
points, all 256 transversals were enumerated over $\mathbb Q(i)$: the cluster through $Q$ has **exactly 8** points,
all differences have degree $\le4$, $z\bar z=N$ holds, and $N$ is squarefree of degree 16. Over $\mathbb F_q$,
$q=19,23,31,43$, every tested nondegenerate $y$ gives exactly 8. ($q=7,11$ have no nondegenerate $y$ with roots in
$\mathbb F_{q^2}$.) Structure: $Q=\tilde Q(t^2)$, where $\tilde Q$ has 4 roots with $|\mathrm{Im}|$ equal, i.e. the
$D=4$ Sperner optimum (6 points), plus one extra non-composed face pair (2 points). The symmetric slice already has
dimension 2. Numerically, the unsymmetrised cube system has a 2-dimensional solution set modulo real affine maps
(`zerosum.py`), as the count $16-2-12$ predicts.

*Corollary C′ (integer polynomial template; exact).* Scaling the roots by 72 gives Gaussian-integer roots and an
8-point template in $\mathbb Z[i][u]$ at bounded $C$, the first polynomial template with more than 6 points
(`families.md` §4). Its constant is huge: chord$/\sqrt R\to3.18\cdot10^9$ before gcd reduction (`cube_template.py`).
It adds a mechanism, not a record: $M_\infty(C)\ge8$ for $C>5.898$ is already known.

**Theorem D (dominating families; proved, characteristic 0).** Fix a combinatorial type: the exponent rows
$\delta_2,\dots,\delta_M\in\{0,\pm1\}^D$ relative to $z_1$, spanning a lattice $\Lambda$ of rank $r$. Let $B$ be an
irreducible variety parametrising root configurations realising this type with given $(D,c)$. Each $b\in B$ gives
$f_b:\mathbb P^1\to T_\Lambda\cong\mathbb G_m^r$, $t\mapsto(g_1,\dots,g_r)$ for a basis of the ratio group, with
$f_b(\infty)=e$. If $\bigcup_bf_b(\mathbb P^1)$ is Zariski dense in $T_\Lambda$, then
$$r(L+1)+1\le2D.$$
In particular, for $c=0$ and even $D$, dominating families have $r\le3$, so $M\le6$ (Prop. 3.3). For $c\ge1$ and
$D\ge10c-1$ they have $r\le4$.
*Proof.* Take a smooth projective toric $X\supset T_\Lambda$ whose fan contains the rays $\pm\nu_k$; here
$\nu_k\in\mathbb Z^r$ is the cocharacter of $f_b$ at $\bar\rho_k$, primitive because the rows have entries in
$\{0,\pm1\}$ (resolve a complete fan containing these rays; the rays survive). Then $f_b$ extends to $\mathbb P^1$ and
meets the boundary transversally, only at $\rho_k,\bar\rho_k$ with $\nu_k\ne0$. So
$-K_X\cdot C=\sum_\rho D_\rho\cdot C\le2D$. With $s=1/t$, Lemma 1.1 says $f_b\equiv e\bmod s^{m}$, $m=L+1$.
Consider the space $H$ of morphisms with $f|_{\mathrm{Spec}\,\mathcal O/s^m}\equiv e$. Its tangent space at $f$ is
$H^0(F)$ with $F=f^*T_X(-m\infty)$. If $H\times\mathbb P^1\to X$ is dominant, generic smoothness (char 0) makes the
differential surjective at a general point. So the subsheaf $G\subseteq F$ generated by $H^0(F)$ has rank
$\ge r-1$, and is $\bigoplus\mathcal O(c_j)$ with $c_j\ge0$. Further, $df$ vanishes to order $\ge m-1$ at $\infty$,
so it is a nonzero section of $F(-1)$ ($=f^*T_X\otimes K_{\mathbb P^1}(-(m-1)\infty)$). If $\mathrm{rk}\,G=r-1$, this
section is generically outside $G$, and $\deg F/G\ge1$. If $\mathrm{rk}\,G=r$, some $c_j\ge1$. Either way
$\deg F\ge1$, i.e. $-K_X\cdot C\ge rm+1$. ∎

*How Theorem D relates to the record.* It is the function-field form of "covering families are excluded"
(`gauss.md` Thm 4.3(ii)) and of the coordinator's (iv). At $c=0$ it allows dominating rank-3 families (4 points:
the analogue of "the 4-point version is false") and forbids dominating rank-4 families (5 or more points: the
analogue of $A_5$ for **dominating** families). The cube family of Theorem C has rank 4 and is a 2-parameter family
of curves in the 4-dimensional $T_\Lambda$. It is not dominating, so it is exactly an **exceptional-set** family.
Composition with $h$ only reparametrises the image curves, so all its members give the same 2-parameter family of
image curves.

**What would decide the function-field problem.** By Theorem A and Odlyzko's bound ($\le2^r$ points of
$\{0,1\}^D$ in an $r$-dimensional subspace), $\sup_DM(D,c)<\infty$ **iff** the rank of the ratio lattice is bounded.
By Prop. 3.3 the following implies $M\le8$ for $D\ge D_1(c)$:

> **(R$_4$)** For $D\ge D_1(c)$ there is no nondegenerate configuration $\rho\in(\mathbb C\setminus\mathbb R)^D$ whose
> moment kernel $V'(\rho)$ contains a rational subspace spanned by $\ge5$ independent $0/1$ vectors of weight
> $\ge L+1$.

(R$_4$) is an **excess-intersection** statement. Its incidence varieties have expected dimension
$2D-2-rL<0$ for $r\ge5$, and Theorem D only shows that they do not dominate. It is the function-field Lang–Vojta
exceptional-set statement for the code subtori $T_\Lambda$: rigid families may exist on special loci, as the rank-4
cube family exists in rank 4. No function-field tool known to me proves it. Mason–Stothers and Wronskians work on
single relations: for squarefree $N$ they give only Lemma 1.2-type bounds, and the $k$-term log-Wronskian gives
$m\le\binom{k+1}2(|S|-1)/k$, which gets weaker as $k$ grows. The truncated second main theorem has an error term
$\asymp|S|=2D$ of the same order as the contact $m\approx D/2$. So its constants decide the outcome: exactly
critical, as in (F2). Deformation theory sees only dominating families. **Verdict: open in the function field too.**

## 5. The steps that use the ultrametric place and the derivation; arithmetic substitutes; the filters

**Which step uses what.**
1. *Ultrametric place plus characteristic 0: Lemma 1.1.* $U_{L+1}$ is a group and $\log$ is an isomorphism onto a
   $\mathbb C$-vector space. Together they turn the arc condition into **$L\approx D/2$ exact linear equations** on the
   exponent vector, with $\mathbb C$-valued coefficients $v(\rho_k)$. This is the step that makes the problem a
   *dimension count* ($rL$ equations against $2D-2$ parameters). Lemma 1.2 and Theorem A use only the ultrametric
   inequality and work in every characteristic.
2. *Derivation: Theorem D only* (generic smoothness, hence tangent vectors to the family). Mason–Stothers, also a
   derivation argument, gives for one pair only $\deg c_{ij}\ge d_{ij}-2n_{ij}+1$. That is vacuous for squarefree
   $N$: the function-field twin of the record's "abc on single relations is vacuous on squarefree profiles"
   (`expert_note.md` §5 item 6).

**What an arithmetic substitute would have to do.**
* *For step 1.* The exact non-archimedean twin of the arc is a congruence $z_i\equiv z_j\bmod\pi^k$ at a split
  prime. There $\log_p$ linearises too, but $U_1/U_k$ of $\mathbb Z_p$ is cyclic, so the condition is **one**
  congruence modulo $p^{k-1}$. It carries the same information as the archimedean arc ($\log p^k\approx\frac12\log R$)
  and gives no dimension count. The function-field strength comes from the residue field $\mathbb C$ being infinite:
  $U_1/U_{L+1}\cong\mathbb C^L$. Over $\mathbb F_q(t)$ the same group is finite, of order $q^L$ on the norm-one part,
  and the problem becomes counting again, with the same heuristic status as over $\mathbb Z$. There is no arithmetic
  object giving $\asymp\log R$ **independent exact** linear conditions per pair.
* *For step 2.* One would need an arithmetic *generic freeness*: a Zariski-dense set of A-points carries
  "tangent directions" that force $-K\cdot C\ge rm+1$. Candidates:
  * Buium $p$-derivations. On $\mathbb G_m$ the $\delta$-characters are Fermat quotients
    $q_p(x)=(x^{p-1}-1)/p$ (memory-quote). Applied to cluster ratios they give, for each auxiliary $p$, one
    $\mathbb F_p$-linear functional $\delta\mapsto\sum_q\delta_q\,q_p(\pi_q/\bar\pi_q)$ on exponent vectors. That is a
    point-level residue modulo $p^2$, inside the residue class that residue-design fakes satisfy for $p^2\le M^A$.
    It carries no dominance information.
  * Vojta's arithmetic jets / Arakelov deformation theory. No unconditional analogue of generic smoothness for
    integral points is known to me.

  Moreover Theorem D concerns **families**, and an integer circle is a single point. Its arithmetic counterpart is a
  statement about Zariski-dense *sets* of A-points, which is precisely $A_5$ / Question B (exceptional set), not a
  pointwise inequality.

**(F1) fakes.** Theorem B gives function-field fakes. They satisfy Lemmas 1.1–1.2 (exact linear and group
structure, all pairs) and have $M\ge D+2c+2$. So in the function field, too, nothing local or linear beats linear
growth. Theorem D does not kill fakes: fakes are not images of algebraic curves, and Theorem D is not an inequality
on a single configuration. The only single-configuration statement that would separate the truth from fakes is the
dimension count, i.e. (R$_4$), which is unproved. **Outcome: the function-field argument does not kill fakes; it
relocates the gap to (R$_4$).**

**(F2) criticality.** The single-relation tools (Mason, Wronskian, truncated SMT) either have error terms
$\asymp|S|=2D$ of the order of the contact $m\approx D/2$, or are vacuous on squarefree $N$. Theorem D is not
critical, since it beats the cube's rank 4 for dominating families, but it applies only to dominating families.
The 4-point test: Theorem D allows dominating rank-3 families, consistent with the record's dense 4-point
families. So it does not "prove the 4-point version", and it passes (F2).

## 6. Files and status of computations

All files are under `scratchpad/round29/ff_checks/`. Run them with `python3`; sympy, numpy and scipy are required.

| file | content | type |
|---|---|---|
| `brute.py` | cyclotomic $t^{2D}+1$, $D\le20$: $M=4$ / $2$ | numerical (rounded hashing) |
| `L2search.py`, `L2multi.py`, `L2rand.py` | $L=2$ integer-box searches (found 4 at $D=6$; superseded by `fakeexact.py`) | search |
| `fakeexact.py` | exhaustive closure enumeration for $L=2$: $M(5,0)=3$, $M(6,0)=5$; $D=7$, $c=1$: $\ge6$ (run `fe_7_2.txt`) | exact enumeration, floating LP |
| `fakeexact_val.py` | validation at $L=1$ (no separation): reproduces Sperner $3,6,10$ | check |
| `rank4codes.py`, `rankcodes_thr.py` | Prop. 3.3 (rank $\le3$: 6, rank $\le4$: 8 at $\theta=\frac12$; 4, 5 at $\theta=\frac47$) | exhaustive + floating LP |
| `zerosum.py` | least-squares realisation of zero-sum families; cube at $D=8$ (2-dim family) | numerical |
| `cube_rational.py` | factorisation $-3(a+b)(a+c)(b+c)(abc(a+b+c)+1)$ | exact (sympy) |
| `cube_exact.py` | Theorem C over $\mathbb Q(i)$: cluster exactly 8, $\deg\le4$, $N$ squarefree | exact |
| `cube_Fq.py` | Theorem C over $\mathbb F_q$, $q=19,23,31,43$: exactly 8 | exact |
| `cube_template.py` | integer template, scale 72, chord$/\sqrt R\to3.18\cdot10^9$ | exact polys, numerical $C$ |
| `gauss8.py` | Gaussian-integer roots in $[-7,7]^2$: at most 2 splits (6 points) | exhaustive in box |
| `d7.py`, `d7b.py` | $D=7$, $c=0$: rank-3 codes realised (4 points), rank 4 not found | numerical (evidence) |
| `rank5.py` | $D=8$: rank-5 spans in the extended Hamming code, none realised | numerical (evidence) |

Notes file: `scratchpad/round29/ff_notes.md`. Literature used (memory-quotes, not re-checked this round):
de Bruijn–Tengbergen–Kruyswijk; Sperner; Plotkin; Odlyzko's $2^r$ lemma; Kollár, *Rational curves on algebraic
varieties* II.3 (free curves, generic smoothness); Mason–Stothers; Buium's $\delta$-characters of $\mathbb G_m$.
