# Round 30, lane search: small cases and the hunt for rank >= 5 in the function-field model

## 0. Verdict and table

**No rank-5 cluster at $c=0$ was found** ($D\le12$; structured loci up to $D=12$). The search found several new
certified objects and proved several combinatorial restrictions. It decided no value of $M(D,0)$ for
$D\ge7$. Status: *no new theorem on (R$_4$)*; new exact and certified calibration. Labels: **proved**,
**certified** (Krawczyk interval certificate), **exhaustive** (complete finite computation), **evidence**.

1. **Balance, strengthened (proved, Lemma 1.1).** The coordinator's (S2) is correct. It holds for every integer
   vector of the real span of the differences, not just for the differences themselves. This kills types by linear
   algebra, e.g. one of the two possible 5-point types at $D=7$.
2. **8-point clusters of rank $\le4$ need $4\mid D$** (Prop. 2.1 proved for the half-cube type; exhaustive MILP over all rank-4
   codes for $D\le16$). At $D\equiv2\pmod4$, $c=0$, codes of rank $\le4$ have at most **6** points ($D=10,14$). So the
   round-29 prediction "$\limsup M=8$ under (R$_4$)" holds only along $4\mid D$, and along $D\equiv2\pmod 4$
   (R$_4$) would give $M\le6$.
3. **A new 8-point rank-4 cluster at $D=12$, not a composition** (Thm 5.1, certified). It is the half-cube code with
   multiplicities $(2,1)$, and the cluster through its base point has exactly 8 points. This answers item 2. At
   $D=10,14$ there are no 8-point clusters of rank $\le4$.
4. **The cube locus has rank exactly 4 everywhere** (Thm 3.1, proved). Its excess intersection comes from
   $\mathrm{Im}\,r^2\equiv2\varepsilon$, and every candidate extra set degenerates.
5. **Chebyshev-system loci (proved, §3.2).** On the horizontal-line (PTE) locus $M\le2$ for $D\ge7$. On circles
   centred on $\mathbb R$, $M\le6$ for $D\le10$.
6. **Rank 5 at $c=1$ beyond the covering range (certified, Thm 6.1)** at $D=10,12$. Random rank-5 types fail at
   $D=14,16$, where the expected dimension is $1,0$.
7. **A real obstruction beyond dimension counting (evidence).** Types of non-negative expected dimension are often
   never realised: the 141 six-point rank-4 codes at $D=10$ (dimension 2), 16 of the 46 five-point types at $D=8$, and rank 5 at $c=1$, $D=14$.

| $D$ | $c$ | largest $M$ found | largest rank found | status of the lower bound | proved upper bound | evidence |
|---|---|---|---|---|---|---|
| 7 | 0 | 4 | 3 | exact (record) | 7 (balanced clique) | $M=4$: the only 5-point type left, $\mathcal T_7$, not realised (3000 starts) |
| 8 | 0 | 8 | 4 | exact (cube) | 14 (balanced clique) | $M=8$: no rank-5 cluster in the hunt over all 46 five-point types (§4) |
| 9 | 0 | 4 | 3 | certified | 6 (balanced clique) | $M=4$: all 6 five-point types fail (600 starts each) |
| 10 | 0 | 5 | 4 | certified | 20 (Plotkin); $\le6$ if rank $\le4$ (exhaustive) | $M=5$: all 141 six-point rank-4 codes fail |
| 12 | 0 | 8 | 4 | certified (Thm 5.1) | 24 (Plotkin) | no rank-5 extension of Thm 5.1 (918 sets) |
| 5, 6 | 1 | 10, 20 | — | exact ($L=1$) | 10, 20 | — |
| 7 | 1 | 8 | — | exact (record) | 14 (balanced clique) | — |
| 8 | 1 | 14 | — | exact (record) | 20 (Plotkin, $A(8,3)$) | — |
| 10, 12 | 1 | $\ge6$ | **5** | certified | — | — |
| 14, 16 | 1 | — | — | — | — | random rank-5 types fail |

The arithmetic growth rate $(\frac23+\varepsilon)\log R/\log\log R$ is unchanged. Nothing here bears on it.

## 1. Verified and corrected coordinator observations (S1, S2, S3, S5, S6)

Notation: $K=\lfloor D/2+c\rfloor$, $L=D-1-K$, $h:=\lceil L/2\rceil$. Roots normalised by $\mathrm{Im}\,\rho_k>0$; a point
is a sign vector $\varepsilon\in\{\pm1\}^D$ ($+1$: factor $t-\rho_k$). "Base form": fix a point $Q=\prod(t-q_k)$; the
cluster through $Q$ is $\{z_F:F\subseteq[D]\text{ zero-sum}\}$, $F$ zero-sum iff $\sum_{k\in F}\mathrm{Im}\,q_k^m=0$
($m\le L$) (Lemma 1.1 of round 29).

**(S1) correct.** If $z_i,z_j$ differ exactly on $S$, $|S|=w$, then $z_i=AB$, $z_j=A\bar B$ with $A=\gcd$, so
$z_i-z_j=2iA\,\mathrm{Im}B$ and the pair condition is $\deg\mathrm{Im}B\le w-L-1$ ($=w-D/2+c$ for even $D$). For
$w=L+1$: $B=P+i\kappa$, $P$ real monic, $\kappa\in\mathbb R^*$, i.e. the $w$ roots of $B$ form a fibre $P^{-1}(-i\kappa)$ of
a real polynomial over a non-real value. (The sign $-2i$ in the coordinator's text is immaterial.)

**(S2) correct, and it extends from differences to the whole saturated lattice (proved).**

> **Lemma 1.1 (balance).** Let $V\subseteq\mathbb R^D$ be the real span of the difference vectors
> $\delta=(\varepsilon_i-\varepsilon_j)/2$ of a cluster (normalised coordinates). Every nonzero $a\in V\cap\mathbb Z^D$ satisfies
> $$\|a\|_1\ge L+1,\qquad a^+:=\textstyle\sum_{a_k>0}a_k\ \ge L/2,\qquad a^-:=\sum_{a_k<0}|a_k|\ \ge L/2 .$$

*Proof.* The moment equations $\sum_ka_k\mathrm{Im}\rho_k^m=0$ ($m\le L$) are real-linear, so they hold on all of $V$. For
$a\in V\cap\mathbb Z^D\setminus0$ put $H=\prod_{a_k>0}(t-\rho_k)^{a_k}\prod_{a_k<0}(t-\bar\rho_k)^{|a_k|}$, monic of degree
$w=\|a\|_1$. Newton's identities (char 0) turn the equations into $e_j(H)=e_j(\bar H)$ for $j\le L$, i.e.
$\deg\mathrm{Im}H\le w-L-1$. $\mathrm{Im}H\not\equiv0$: otherwise the root multiset of $H$ is closed under conjugation,
impossible since it uses at most one of $\rho_k,\bar\rho_k$ for each $k$ and the $2D$ numbers $\rho_k,\bar\rho_k$ are
distinct. Hence $w\ge L+1$. On $\mathbb R$, $H/\bar H=e^{i\Theta(t)}$ with
$\Theta=2\sum_ka_k\arg(t-\rho_k)$, $\arg(t-\rho_k)\in(-\pi,0)$ increasing from $-\pi$ to $0$. So $\Theta$ is continuous from
$-2\pi(a^+-a^-)$ to $0$, and by the intermediate value theorem it meets $2\pi\mathbb Z$ at least $|a^+-a^-|-1$ times in
the open line; each meeting is a real zero of $\mathrm{Im}H$. So $|a^+-a^-|-1\le w-L-1=a^++a^--L-1$, i.e.
$\min(a^+,a^-)\ge L/2$. ∎

The coordinator's statement is the case $a=\delta$. Two consequences used below: (i) the condition is not
pairwise but applies to every integer vector of $V$ (in particular to sums of differences); (ii) the one-signed
case ($a^-=0$) forces $L\le0$, the "Blaschke" remark. The flip symmetry $\varepsilon_k\mapsto-\varepsilon_k$ of the
linear equations is **broken** by the normalisation, so types must be counted up to coordinate permutations and the
global flip $\varepsilon\mapsto-\varepsilon$ only (the latter is $z\mapsto\bar z$).

**(S3) correct** but equivalent to Lemma 1.1 of round 29: $\eta_a=d\log F_a$ with
$F_a=\prod((t-\rho_k)/(t-\bar\rho_k))^{a_k}$; $W_a=2\sum_ka_ky_k\prod_{l\ne k}N_l$ is real, residue $a_k$ at $\rho_k$,
and $F_a=1+O(t^{-L-1})$ iff $\deg W_a\le2D-L-2$. Not used further.

**(S5) correct.** For $\rho_k=x_k+i$, $\mathrm{Im}\rho^m=\sum_{j\text{ odd}}\binom mj(-1)^{(j-1)/2}x^{m-j}$ is triangular with
leading term $mx^{m-1}$, so the equations are $\sum\delta_kx_k^j=0$, $0\le j\le L-1$.

**(S6) correct with a caveat.** It needs roots in $\mathbb Q(i)$ (after scaling, in $\mathbb Z[i]$). The new cluster of §5
is certified over $\mathbb R$ (real algebraic numbers), not over $\mathbb Q(i)$, so it is **not** shown to be a template.

**(S4)** is a heuristic about a conjectural inequality; I did not use it and do not assess it.

## 2. Combinatorial filters: balanced codes, codes of rank $\le4$, 5-point types

**2.1 Balanced code bound (proved; exhaustive clique computation).** By Lemma 1.1 applied to differences, a cluster
is a set of sign vectors with, for every pair, $\#\{k:\varepsilon_{ik}=+,\varepsilon_{jk}=-\}\ge h$, the same with $i,j$
swapped, and Hamming distance $\ge L+1$. The exact maximum of such sets (`clique.py`, igraph exact clique number on
the $2^D$ vertices) is:

| $D$ ($c=0$) | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|
| $L$, $h$ | 2, 1 | 2, 1 | 3, 2 | 3, 2 | 4, 2 |
| balanced bound | 4 | 7 | 7 | **14** | **6** |
| without balance ($A(D,L+1)$) | 4 | 8 | 8 | 16 | — |
| true $M(D,0)$ | 3 | 5 | 4 (evidence, §4) | $\ge8$ | $\ge4$ (§4) |

So balance improves Plotkin to $M(8,0)\le14$ and gives $M(9,0)\le6$ (proved). These bounds use pairs only.

**2.2 Codes of rank $\le4$ (exhaustive; MILP).** As in round-29 Prop. 3.3, a base-form cluster of real rank $r$ is
$C\cong S\subseteq\{0,1\}^r$ ($0\in S$, $S$ spanning) with each coordinate a linear functional $\ell$ with
$\ell(S)\subseteq\{0,1\}$. I enumerated **all** labelled $S$ and solved, for each, an integer program for the column
multiplicities (each functional split by the sign of $\mathrm{Im}\,q_k$), with constraints: length $D$, distances
$\ge L+1$, and balance $\ge h$ in both directions for every pair (`codes4.py`; HiGHS MILP, exact integer data).
Maximal $|C|$ at $c=0$:

| $D$ | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 16 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| rank $\le3$ | 4 | 4 | 4 | 6 | 4 | 4 | 4 | 6 | 4 | 4 | — |
| rank $\le4$ | — | 6 | 5 | 8 | 5 | **6** | — | 8 | — | **6** | 8 |

(At $D=7,9,10,14$ the rank-4 values are the same without balance; at 8 points ($D=8,12,16$) exactly 29
labelled $S$ are feasible, and for each of them the MILP solution is antipodal, with all non-antipodal distances equal to
$D/2$, i.e. of half-cube form. I did not enumerate all multiplicity solutions per $S$, so "every 8-point rank-4 code
is half-cube" is verified for the solutions found, not proved.)

> **Proposition 2.1 (half-cube codes need $4\mid D$; proved).** Let $D$ be even, $c=0$ (so the distance is $L+1=D/2$).
> The half-cube code (coordinates labelled by $u\in\mathbb F_2^3$ with multiplicities $m_u$, words
> $\emptyset,[D],H_i=\{u_i=0\},H_i^c$) has all distances $\ge D/2$ iff $m_u=(D+\lambda\chi(u))/8$ with
> $\chi(u)=(-1)^{u_1+u_2+u_3}$ and an integer $\lambda$; integrality forces $D+\lambda\equiv D-\lambda\equiv0\pmod 8$,
> hence $4\mid D$. For $D=8$: $m\equiv1$ (the cube). For $D=12$: $m_u\in\{2,1\}$.

*Proof.* $|H_i|$, $|H_i^c|$, $|H_i\triangle H_j|$ and $|H_i\triangle H_j^c|$ are all $\ge D/2$ and pairwise complementary
in $[D]$, so all equal $D/2$. Thus $\sum_um_u(-1)^{a\cdot u}=0$ for the six $a\ne0,(1,1,1)$; Fourier inversion
gives $m_u=(D+\lambda\chi(u))/8$ with $\lambda=\sum_um_u\chi(u)$. ∎

**Consequence (conditional, cleanly stated).** If (R$_4$) held at $D\equiv2\pmod4$, $c=0$, then $M(D,0)\le6$ there (for
$D=10,14$ by the table). So the round-29 prediction "$\limsup M=8$" can only come from $D\equiv0\pmod4$: the LP
value 8 at $\theta=\frac12$ is not attained at $D\equiv2\pmod4$, because the extremal code needs integral
multiplicities.

**2.3 Five-point types (exhaustive).** `fftypes.py` enumerates all sets of $M$ sign vectors in $\{\pm1\}^D$ with
pairwise balance $\ge h$ and distance $\ge L+1$, up to $S_D$, $S_M$ and the global flip. Then I apply Gordan (the $m=1$
equation needs $y>0$ in the kernel) and Lemma 1.1 on the **saturation** $V\cap\mathbb Z^D$ (LLL basis via python-flint,
all combinations with coefficients in $[-2,2]$; a violation is a proof that the type is not realisable). Counts:

| $D$, $M$ | types | after Gordan | after saturated balance |
|---|---|---|---|
| 7, 5 | 2 | 2 | **1** |
| 8, 5 | 65 | 65 | 46 |
| 9, 4 | 29 | 28 | (21 by lattice filter) |
| 9, 5 | 6 | 6 | 6 |

The killed $D=7$ type is instructive. Its four nonzero base-form sets are $\pi_2\cup\pi_3,\pi_1\cup\pi_3,\pi_1\cup\pi_2$
and one more, with $\pi_1,\pi_2,\pi_3$ disjoint pairs. Linearity then makes each pair $\pi_i$ zero-sum. That is a
saturation vector of weight $2<L+1$, which is impossible. Any 6-point cluster contains a 5-point cluster, so:

> **Corollary 2.2.** $M(7,0)\ge5$ iff the single type $\mathcal T_7$ below is realisable; $M(9,0)\ge5$ iff one of six
> explicit types is realisable. (Proved reductions; the types are in `types_D7_M5_c0.pkl`, `types_D9_M5_c0_sat.pkl`.)
> $$\mathcal T_7:\ \begin{smallmatrix}-&-&-&-&+&+&+\\-&-&+&+&-&-&+\\-&+&-&+&-&+&-\\+&-&-&+&+&-&-\\+&+&+&-&-&-&-\end{smallmatrix}$$

## 3. Exact result on a positive-dimensional locus: the cube locus has rank exactly 4

Cube locus (round 29, Thm C): $y=(a,b,c,d)\in(\mathbb R^*)^4$ with $d=-(a+b+c)$ and $abcd=1$,
$\varepsilon=(1,-1,-1,1)$, $r_k=\varepsilon_k/y_k+iy_k$, $Q=\prod_k(t^2-r_k^2)$, $D=8$, $L=3$. Call a point
**nondegenerate** if the 16 numbers $\pm r_k,\pm\bar r_k$ are distinct (i.e. $N$ is squarefree).

> **Theorem 3.1 (proved).** At every nondegenerate real point of the cube locus, the cluster through $Q$ is exactly
> the 8 points of Thm C. So the rank is exactly 4 on the whole 2-dimensional locus, and no point of the locus extends
> to a rank-5 or 9-point cluster.

*Proof.* A base-form set $F$ is given by $(\alpha_k,\beta_k)\in\{0,1\}^2$ ($r_k\in F$, $-r_k\in F$). Put
$s_k=\alpha_k-\beta_k$ and $t_k=\alpha_k+\beta_k$. Then $\mathrm{Im}\,r_k=y_k$, $\mathrm{Im}\,r_k^2=2x_ky_k=2\varepsilon_k$ and
$\mathrm{Im}\,r_k^3=3x_k^2y_k-y_k^3=3/y_k-y_k^3$, with odd moments changing sign under $r\mapsto-r$. So $F$ is zero-sum iff
$$\textstyle(m{=}1)\ \sum s_ky_k=0,\qquad(m{=}2)\ \sum t_k\varepsilon_k=0,\qquad(m{=}3)\ \sum s_k(3/y_k-y_k^3)=0 .$$
The $m=2$ condition is **combinatorial** on this locus. That is the excess-intersection mechanism of the cube.
Case $s=0$: then $t\in\{0,2\}^4$, and $m=2$ gives the 6 symmetric sets. Case $s=\pm(1,1,1,1)$: $t=\mathbf 1$; these are
the 2 transversals (Thm C). Every other $s\ne0$ ($80-2$ vectors; `cube_locus.py` lists those passing $m=2$, and in
each the support of $s$ has size 2 or 4, since support 3 makes $\sum t\varepsilon$ odd) falls in one of three cases.
* $|\mathrm{supp}\,s|=4$ with three equal signs: then $m=1$ reads $2y_j=0$ for some $j$, which is impossible.
* $|\mathrm{supp}\,s|=4$ with two signs each: then $m=1$ and $\sum y=0$ give $y_i+y_j=0$ for a pair.
* $|\mathrm{supp}\,s|=2$: then $m=1$ gives $y_j=\pm y_i$.

In both of the last two cases $y_j=\pm y_i$. Check the four sign cases with $x=\varepsilon/y$:
* $y_j=y_i$, $\varepsilon_j=\varepsilon_i$ gives $r_j=r_i$;
* $y_j=y_i$, $\varepsilon_j=-\varepsilon_i$ gives $-r_j=\bar r_i$;
* $y_j=-y_i$, $\varepsilon_j=\varepsilon_i$ gives $-r_j=\bar r_i$;
* $y_j=-y_i$, $\varepsilon_j=-\varepsilon_i$ gives $r_j=\bar r_i$.

Each is a degeneracy. ∎

(The sympy factorisations in `cube_locus.py` confirm independently that, for every listed $s$, the $m=3$ numerator
is divisible by the $m=1$ form or the $m=1$ form is a single $y_j$.)

**3.2 Loci where the moment rows form a Chebyshev system (proved).** Suppose that on a locus
$\mathrm{Im}\,\rho^m=w(\rho)\,p_{m-1}(\varphi(\rho))$ for $m\le L$, with $w>0$, $\varphi$ real and injective on the
roots, and $p_{m-1}$ of exact degree $m-1$. Then for every lattice vector $a$, the signed measure
$\sum_ka_kw(\rho_k)\delta_{\varphi(\rho_k)}$ annihilates all polynomials of degree $<L$. So $a$ has at least $L$ sign changes
when the roots are ordered by $\varphi$: otherwise a polynomial of degree $\le L-1$ with the same sign pattern would have
positive integral.
* *Horizontal line* $\mathrm{Im}\rho=1$ (S5): $w=1$, $\varphi=x$. Here the $j=0$ row also forces $a^+=a^-$. Then $a^+\ge L$:
  if $a^+\le L-1$, Newton makes the two root sets equal. So every difference has weight $\ge2L=D-2$ at $c=0$ (even $D$).
  Plotkin ($2d>n$) gives $M\le2\lfloor (D-2)/(D-4)\rfloor=2$ for $D\ge7$. **On the horizontal-line (PTE) locus,
  $M\le2$ for $D\ge7$, $c=0$** (Plotkin gives $\le4$ at $D=6$; for odd $D$ the weight is $\ge2L=D-1$ and again $M\le2$). That is a decisive negative for (S5) as a source of clusters: ideal
  PTE solutions have size $\approx$ degree, while $c=0$ needs degree $L-1\approx D/2$ inside a $D$-set.
* *Circles centred on $\mathbb R$* (affine images of $|\rho|=1$): $\rho=e^{i\theta}$, $w=\sin\theta$, $\varphi=\cos\theta$,
  $p_{m-1}=U_{m-1}$. Every difference needs $\ge L$ sign changes in the angular order, plus balance and distance.
  Exact clique numbers (`circle_clique.py`): $M\le6,4,6,4,6$ for $D=6,7,8,9,10$, $c=0$. So **the circle locus carries no
  cluster with more than 6 points for $D\le10$**, and in particular not the cube.
* *Vertical lines* $\rho=x_0+iy$: only $\lceil L/2\rceil$ independent rows, so $r\lceil L/2\rceil\le D-1$ and $r\le3$ at
  $c=0$. This agrees with kill_ff's "no excess for $L\ge3$".

**Remark (scope of 3.1).** Theorem 3.1 is a rank bound on a positive-dimensional locus, but only on the cube locus. It is
not shown that every realisation of the half-cube type at $D=8$ lies on that locus. Evidence that it does
(`halfcube8_sym.py`): all 8 nondegenerate realisations found from 400 starts are $t\mapsto-t$ symmetric after
centring, to $10^{-10}$. Combined with §2.2 (every 8-point rank-$\le4$ code at $D=8$ is half-cube) the conditional
statement is: **if** the half-cube locus is the cube locus, then $D=8$, $c=0$ has no 9-point cluster of rank $\le5$
containing an 8-point rank-4 subcluster.

## 4. Numerical type search, $c=0$ (task items 1 and 3)

**Method.** Base form, unknowns $q_k=x_k+iy_k$ (any half-plane); residuals $\sum_{k\in F}\mathrm{Im}\,q_k^m$ for the
type's sets; gauges $\sum x_k=0$, $\sum|q_k|^2=D$; Levenberg–Marquardt (trf when underdetermined), starts with
$\mathrm{sign}\,y_k$ prescribed by the type. A solution counts if the residual is $<10^{-10}$ and the minimal distance
among the $2D$ numbers $q_k,\bar q_k$ is $>10^{-4}$ (`basesolve.py`, `realise.py`). Every existence claim was then
Krawczyk-certified (`certify.py`: 60-digit Newton, then a box of radius $10^{-30}$, with surplus coordinates fixed to
rationals). *Calibration:* a first solver with $y_k=e^{u_k}$ (`ffsolve.py`) realised 0 of the 21 four-point types at
$D=9$. The base-form solver realised 20 of 21, and it recovers the cube family. I discarded all results of the
first solver.

**Results.**

| target | types | starts/type | outcome |
|---|---|---|---|
| $D=7$, 5 points (rank 4, exp. dim 0) | 1 ($\mathcal T_7$; the other is killed by saturation) | 3000 | **none**; the best runs converge to a degeneration (separation $3\cdot10^{-5}$) |
| $D=7$, 4 points | 5 | 30 | 3 realised (the record's exact $M(7,0)\ge4$) |
| $D=9$, 4 points (rank 3, exp. dim 4) | 21 | 60 | 20 realised; one **certified** (`cert_D9_M4_c0_t0.npy`; $\min|\mathrm{Im}|=0.144$, sep. $0.123$) |
| $D=9$, 5 points (rank 4, exp. dim 0) | 6 (all) | 600 | **none**; all best runs degenerate (sep. $\le6\cdot10^{-9}$) |
| $D=10$, 6 points, rank 4 (exp. dim 2) | 141 labelled codes (all MILP-feasible $S$, several multiplicity solutions each) | 40 | **none** |
| $D=10$, 5 points (rank 4, exp. dim 2) | 30 random types | 30 | 25 realised; one **certified** (`rt_D10_c0_M5_1`; cluster through $Q$ exactly 5, $\min|\mathrm{Im}|=0.091$, sep. $0.088$) |
| $D=10$, 4 points | 10 random types | 20 | 10 realised |
| $D=12$ half-cube + one more set (rank 5) | 918 sets $F$, $\lvert F\rvert=6$ | 2 seeds on the family | **none** |
| $D=8$: every 5-point type (46) + one more set (rank 5) | 46 types; 30 realised | $\le3$ family seeds + 2–4 random | **no genuine rank 5** (all hits spurious, see below) |

The $D=8$ line is exhaustive in type space. Any rank-5 cluster contains 5 points with independent differences (a
rank-4, 5-point subcluster, which is one of the 46 types) plus one further point, which is one extra set $F$ in base
form. Its weakness is numerical: each type's 2-dimensional family is sampled at a few points. Types 0–9 were run with
all $F$. For the rest, $F$ was first pruned by the necessary conditions: the new point is balanced against all 5
points, and Lemma 1.1 holds on the saturated span with coefficients in $\{-1,0,1\}$ (`ext_types2.py`). Of the 46
types, 16 were never realised in 30 starts, although their expected dimension is 2.

*Hits and why none is genuine* (`check_hits.py`, `limit_check.py`).
* 7 hits extend a **rank-3** five-point type. Polished, they are cubes: exactly 8 zero-sum sets, rank 4.
* The other 8 hits are rank 5, but each has separation $pprox9\cdot10^{-5}$ and residual $pprox10^{-12}$.
  * Forcing separation $\ge s$ raises the minimal residual to $5\cdot10^{-9},10^{-6},1.5\cdot10^{-4},9\cdot10^{-2}$ for
    $s=10^{-3},10^{-2},5\cdot10^{-2},10^{-1}$.
  * Snapping the roots to the four limit points $e^{\pm i\pi/4},e^{\pm3i\pi/4}$ makes all six points **coincide**, so the
    pair differences vanish identically (exact sympy).

  So these are degenerations of a single point, not clusters. Under the stated method, **no rank-5 cluster exists at
  $D=8$, $c=0$** (evidence).

## 5. Rank-4 8-point clusters at $D=10,12,14$ (task item 2)

**$D=10$ and $D=14$: none (proved).** By §2.2 (exhaustive over $S$, MILP), a rank-$\le4$ cluster at $D=10$ or $14$, $c=0$,
has at most 6 points, with or without balance. The half-cube code is excluded by hand (Prop. 2.1). The MILP is
floating-point HiGHS on 0/1 data; I label this "exhaustive computation" rather than a pen-and-paper proof.

**$D=12$: a new 8-point rank-4 cluster that is not a composition (certified).**

> **Theorem 5.1 (Krawczyk-certified).** There is a squarefree real $N$ of degree 24 and an 8-point cluster at $D=12$,
> $c=0$, of rank 4. Its type is the half-cube code with multiplicities $m_u=2$ for $u_1+u_2+u_3$ even and $m_u=1$ for
> $u_1+u_2+u_3$ odd. The cluster through its base point $Q$ has exactly 8 points. It is not a composition $z\circ h$ with
> $\deg h\ge2$.

*Certificate* (`halfcube12.py`, `cert_hc12.py`, `fullcluster_mp.py`).
* Unknowns: the 24 real coordinates of $q_1,\dots,q_{12}$.
* Equations: the 20 equations of the 4 independent sets $[12],H_1,H_2,H_3$; the gauges $\sum\mathrm{Re}\,q_k=0$ and
  $\sum|q_k|^2=12$; and $\mathrm{Re}\,q_1=-1.835$, $\mathrm{Re}\,q_2=-0.705$ (exact rationals).
* Krawczyk test: the operator maps the box of radius $10^{-30}$ about the 60-digit Newton refinement into its
  interior. Hence there is a unique exact solution in the box.
* Nondegeneracy: on it $\min|\mathrm{Im}\,q_k|=0.0261$ and the minimal distance among the 24 numbers $q_k,\bar q_k$ is
  $0.050$, so $N$ is squarefree.
* The 4 remaining sets $H_i^c=[12]\setminus H_i$ follow by linearity.
* Exactly 8: all $4096$ subsets were evaluated at 60 digits. The non-members have residual norm $\ge0.118\gg10^{-30}$,
  so the cluster through $Q$ is exactly these 8.
* Not a composition: composing a $(D',0)$-cluster with $\deg h=12/D'$ preserves the size, and
  $M(6,0)=5$, $M(4,0)=6$, $M(3,0)=3$, $M(2,0)=4$ are all $<8$.

The numerical solution set has the expected dimension $24-2-20=2$: 4 nondegenerate solutions out of 300 starts,
all distinct modulo the gauge. A second, independent certified solution (`cert_hc12.py`) has
$\min|\mathrm{Im}|=0.157$ and separation $0.031$.

*Answer to item (2).* At $D=8$ every 8-point rank-$\le4$ code found by the exhaustive-in-$S$ search is half-cube (§2.2); at $D=8$ the half-cube multiplicities are forced ($m\equiv1$), so this is the cube type. The cube
locus is a 2-parameter real component of it (dimension 2 modulo affine maps, matching $16-2-12$). That the half-cube
locus is **only** the cube locus is evidence (§3 remark), not proved. At $D=10,14$ there are no 8-point clusters of
rank $\le4$. At $D=12$ there is the certified non-composed family of Theorem 5.1. In general: **8-point rank-4
clusters at $c=0$ exist only for $4\mid D$** ($D\le16$ by computation, and for every $D$ for the half-cube type by
Prop. 2.1). They do exist at $D=8,12$ (and at $16$ by composition).

*Rank-5 extension of Theorem 5.1 (evidence).* `extend.py hc12` adds each of the extra sets $F\subseteq[12]$, independent
of the span and with $L+1\le|F|\le D-L-1$, to the half-cube system, and solves from two seeds on the family. See §4
for the outcome.

## 6. $c=1$ (task item 4)

**Small values.** $D=5,6$: $L=1$, so $M=\binom D{\lfloor D/2\rfloor}=10,20$ (round 29, Prop. 3.1). For $L=1$ the
balanced clique number is also 10 and 20, consistent with this. $D=7$ ($L=2$): the record's exact lower bound is
$M(7,1)\ge8$, and the balanced clique bound gives $M(7,1)\le14$ (proved; Plotkin gives 16). $D=8$ ($L=2$): the
record's exact lower bound is $\ge14$; the clique computation did not finish. $L=2$ is the free-column problem, so
exact values need the round-29 closure enumeration, which I did not repeat.

**Rank $\ge5$ beyond the covering range (certified).** Theorem D (covering families) allows rank 5 at $c=1$ only for
$D\le8$. The dimension count $2D-2-5L$ is $8-D/2$ for even $D$, so non-covering rank-5 families should exist up to
$D=16$. Random 6-point types (pairwise balanced, Gordan, saturated balance with coefficients $\le1$; 20 starts each;
`randtypes.py`):

| $D$ ($c=1$) | $L$ | exp. dim (rank 5) | types tried | realised | certified |
|---|---|---|---|---|---|
| 10 | 3 | 3 | 15 | 10 | 2 (`rt_D10_c1_M6_{1,10}`): Krawczyk, cluster through $Q$ exactly 6, rank 5 |
| 12 | 4 | 2 | 15 | 7 | 1 (`rt_D12_c1_M6_1`): Krawczyk, exactly 6, rank 5; $\min|\mathrm{Im}|=0.026$, sep. $0.0085$ |
| 14 | 5 | 1 | 15 | 0 | — |
| 16 | 6 | 0 | (see §7) | | |

> **Theorem 6.1 (certified).** There are rank-5 clusters (6 points) at $c=1$ for $D=10$ and $D=12$. Sporadic
> (non-covering) rank-5 families therefore exist beyond the covering range $D\le8$ of Theorem D.

So at $c=1$ the expected-dimension count is a good guide while the dimension is $\ge2$. At expected dimension 1
($D=14$) random types already fail, so the real (semi-algebraic) constraints bite before the complex count does.
The same happens at $c=0$: the rank-4 six-point codes at $D=10$ have expected dimension 2 and were **never**
realised (141 codes), while the rank-4 five-point types at the same $D$ realise 25 times out of 30.

## 7. Stop rule, what was not done, files

**Pre-registered stop rule.** *No rank-5 cluster at $c=0$ for $D\le12$:* none was found (§4).
*No rank bound on any positive-dimensional locus:* this half is **formally not met**. Theorem 3.1 (cube locus: rank
exactly 4) and §3.2 (horizontal-line locus: $M\le2$; circle locus: $M\le6$ for $D\le10$) are proved bounds on
positive-dimensional loci. In substance they are weak. The cube result says only that the cube's excess is used up
by the 8 points, and the line/circle results use a Chebyshev-system structure that general configurations lack.
None of them touches a locus on which rank 5 was expected. My recommendation: treat the rule as **triggered in
spirit**. Close the function-field line for internal work unless someone proposes an argument that uses the
real-algebraic obstruction seen in §4–§6: types with non-negative expected dimension that are nevertheless never
realised (the $D=10$ six-point codes, most $D=8$ five-point types, $c=1$ at $D=14,16$). That obstruction is the only
new phenomenon here, and it is not understood.

**Not done.**
* The balanced clique bound at $D=10$ (stopped after 34 min).
* Exact closure enumeration at $c=1$, $L=2$.
* A Gröbner or homotopy certificate that $\mathcal T_7$ and the six $D=9$ types have only degenerate solutions:
  no Gröbner or homotopy software was available, and `pip install phcpy` failed.
* Composition with $t^n$ and dihedral orbit loci beyond the quadratic case already in the record.
* Belyi or Chebyshev root loci.

**Files** (all under `scratchpad/round30/search_checks/`).

| file | content | type |
|---|---|---|
| `codes4.py` | rank $\le r$ code maxima with balance (MILP over all $S$) | exhaustive (floating MILP) |
| `clique.py`, `circle_clique.py` | balanced code bounds; circle-locus bounds | exact clique numbers |
| `fftypes.py`, `satlat.py` | type enumeration; saturated-lattice balance filter (python-flint LLL) | exhaustive / exact |
| `basesolve.py`, `realise.py`, `randtypes.py` | LM realisation of types | numerical |
| `certify.py`, `cert_generic.py`, `cert_rt.py`, `cert_hc12.py`, `fullcluster_mp.py` | Krawczyk certificates, 60-digit cluster counts | certified |
| `halfcube12.py`, `hc12_lam4.npy` | $D=12$ 8-point cluster (Thm 5.1) | certified |
| `cube_locus.py` | Thm 3.1 equations | exact (sympy) |
| `halfcube8_sym.py` | half-cube realisations at $D=8$ are symmetric | numerical |
| `extend.py`, `ext_types.py` | rank-5 extension hunts ($D=12$ half-cube; all $D=8$ five-point types) | numerical |
| `*.txt` | run logs | |

Notes: `scratchpad/round30/notes.md`. Literature: none newly cited. Plotkin's bound and the Chebyshev-system
sign-change lemma are standard memory-quotes, the latter with the one-line proof given in §3.2.
