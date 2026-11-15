# Round 30, lane prove: rank bounds in the function-field model

Notation as in round 29 (`ff.md`). We have $c=0$ unless stated otherwise. The roots $\rho_k=x_k+iy_k$ are
normalised to $y_k>0$, and $\phi_k(t)=\arg(t-\rho_k)$ for real $t$. A point is $s\in\{\pm1\}^D$, with $s_k=+1$ meaning
that the point uses $\rho_k$. Put $L=D-1-\lfloor D/2+c\rfloor$, $v(\rho)=(\mathrm{Im}\rho^m)_{m\le L}$, and let
$V'=\ker[v(\rho_1)|\dots|v(\rho_D)]\subseteq\mathbb R^D$ be the real moment kernel.
Labels: **proved**, **exact** (finite exact computation), **certified** (exact identity plus numerics with a
rigorous margin), **evidence**, **sketch**, **memory-quote**.

## 0. Verdict

**Status: obstruction plus three locus theorems. No general rank bound and no rank-5 cluster.**

1. **(S1)–(S3) verified, with corrections** (§1).
   * (S1) holds, up to a sign.
   * (S2) holds and extends from difference vectors to **every** nonzero integer vector $a$ of the moment kernel:
     $\|a\|_1\ge L+1$ and $\min(a^+,a^-)\ge L/2$ (Prop 1.2).
   * (S3) holds for the differentials. Its "factorisation by types" is **false as stated**: the near-identities hold for
     the Blaschke factors $n_v/\bar n_v$, not for the norms $N_v$. An exact counterexample is on the cube.
   * (S4) is the truncated Vojta inequality for $(\mathrm{Bl}_eX,\partial)$, and its arithmetic is right. But
     $K+\partial=(r-1)E$ has Iitaka dimension 0, so the conjecture puts no a-priori limit on the size of the
     exceptional set. The known unconditional function-field inequalities have error $\asymp|S|$ (§4.4).
2. **Real structure (§2).**
   * Phases are monotone and the Cauchy index is bounded, which gives balance.
   * Minimal pairs are fibres $P^{-1}(-i\kappa)$ of real polynomials.
   * Triples obey a symplectic (Plücker) identity among block polynomials. Mutually minimal triples give a real
     pencil relation.
   * Equidistant clusters have totally ordered phase functions.

   All of this is proved, and it is real non-fake structure. But **balance does not separate**: doubled fakes have
   $\ge D/2$ points with every difference exactly balanced (Prop 2.5, proved). I found no single-configuration
   inequality that fakes violate *and* that bounds $M$.
3. **Positive-dimensional loci.**
   * (a) **Cube locus (certified):** every nondegenerate real member has largest cluster **exactly 8**, over all
     256 base points, and rank exactly 4. There are 48 special real members (one $|y|$-set up to symmetry), and they
     are no larger.
   * (b) **Horizontal-line locus (proved):** $M\le2$ (rank $\le1$) for $D>6c+6$. The PTE structure forces the pair
     distance to be $\ge 2L$, not $L+1$.
   * (c) **Minimal-distance (equidistant) clusters:**
     * $M\le D$ (proved), with the structure of §2.
     * 5-point rank-4 equidistant clusters exist at $D=8$ in 2-parameter families (evidence, residual $10^{-51}$).
     * No nondegenerate 6-point (rank-5) equidistant cluster at $D=8$ (evidence).
     * Classification beyond this: open.
4. **General tools (§4), each labelled.**
   * The Wronskian / $W$-space is **provably blind**: integrality enters only through $V'\cap\mathbb Z^D$.
   * Simultaneous Mason is vacuous.
   * The deformation / a-invariant stratification gives at most $\dim Y\le3$ for swept-out loci (sketch). It is
     blind to rigid configurations, which are exactly the expected rank-$\ge5$ ones.
   * Lehmann–Riedl–Tanimoto's closedness theorem is for Fano $X$ without contact conditions (abstract-level check).
     It does not apply.
   * Gcd theorems: $h\le C|S|$, vacuous.
5. **Partial statements (§5):**
   * rank $\le1$ on horizontal lines;
   * $M=8$ on the cube;
   * equidistant $M\le D$, sharpened by structure only.

   **Nothing beats Plotkin in general.**
6. **Stop rule (§6).** This lane found no rank-5 cluster (evidence only, $D=8$ equidistant). It proved rank bounds only
   on two positive-dimensional loci, (a) and (b), and these are *special* loci: the cube surface and the
   horizontal-line $D$-fold. I read the pre-registered rule literally. "Proves no rank bound on any positive-dimensional
   locus" is not met, so the rule does not formally fire. My recommendation is nonetheless to **close the line for
   internal work**. The two locus theorems close by elimination and by PTE combinatorics, and neither mechanism
   extends to general position (§6).

## 1. Verification of (S1)–(S4)

**Lemma 1.1 (S1; proved).** Let $z_i\ne z_j$ differ on $S$, $|S|=w$. Let $A=\gcd(z_i,z_j)$ and let $B=\prod_{k\in S}$ (the
factor of $z_i$ at $k$). Then $z_j=A\bar B$ and
$$z_i-z_j=A(B-\bar B)=2i\,A\,\mathrm{Im}B,$$
so $\deg\mathrm{Im}B\le w-L-1$ ($=w-D/2+c$ for even $D$). If $w=L+1$, then $B=P+i\kappa$ with $P$ real monic and
$\kappa\in\mathbb R^*$. *Proof.* On $S$, $z_j$ uses the conjugate root, so its factor product is $\bar B$. Degrees add:
$\deg(z_i-z_j)=D-w+\deg\mathrm{Im}B\le K$. $\mathrm{Im}B\ne0$, since $B$ cannot be conjugation-invariant: it contains
exactly one root of each pair. ∎ The coordinator's $-2i$ is a sign convention.

**Lemma 1.2 (S2 and its lattice form; proved).** For real $t$, $\phi_k'(t)=y_k/|t-\rho_k|^2>0$, and $\phi_k$ increases
from $-\pi$ to $0$. Let $a\in\mathbb Z^D\cap V'$, $a\ne0$. Put
$U_a=\prod_{a_k>0}(t-\rho_k)^{a_k}\prod_{a_k<0}(t-\bar\rho_k)^{|a_k|}$, of degree $\|a\|_1$. Then:
1. $\deg\mathrm{Im}U_a\le\|a\|_1-L-1$, hence $\|a\|_1\ge L+1$;
2. $\arg U_a=\Psi_a:=\sum a_k\phi_k$ runs from $-\pi\sum a_k$ to $0$;
3. $|\sum_ka_k|\le\|a\|_1-L$, i.e. $\min(a^+,a^-)\ge L/2$.

*Proof.*
1. $U_a/\bar U_a=\prod((t-\rho_k)/(t-\bar\rho_k))^{a_k}=\exp(\sum_ka_k\log\frac{1-\rho_k/t}{1-\bar\rho_k/t})$, and by
   Lemma 1.1 of `ff.md` (Newton) this is $1+O(t^{-L-1})$. So $U_a-\bar U_a=O(t^{\|a\|_1-L-1})$. $U_a$ is not real,
   because its roots $\rho_k$ ($a_k>0$) and $\bar\rho_l$ ($a_l<0$) are never conjugate to each other.
2. A root in the upper half-plane contributes $+\pi$ to the change of argument and a root in the lower half-plane
   contributes $-\pi$.
3. $\mathrm{Im}U_a=|U_a|\sin\Psi_a$ vanishes wherever $\Psi_a\in\pi\mathbb Z$. By the intermediate value theorem
   this happens at $\ge|\sum a_k|-1$ points. So $|\sum a|-1\le\|a\|_1-L-1$.

∎ For $a=\delta\in\{0,\pm1\}^D$ this is exactly (S2): $|p-n|\le w-L$. A one-signed $\delta$ is impossible for $L\ge1$
(for $L=0$, $p-n=w$ is allowed). The coordinator's "contact exactly 1" is the special case $\deg\mathrm{Im}B=w-1$
(Hermite–Biehler: all roots of $B$ lie in one half-plane, so $\mathrm{Re}B$ and $\mathrm{Im}B$ interlace and
$\deg\mathrm{Im}B=w-1$). That means contact $L+1=1$.

**Lemma 1.3 (S3; proved, with a correction).** For $a\in\mathbb R^D$:
* $\eta_a=\sum a_k(\frac{dt}{t-\rho_k}-\frac{dt}{t-\bar\rho_k})=2i\Psi_a'(t)dt=i\,W_a\,dt/N$, where
  $W_a=2N\sum_ka_ky_k/|t-\rho_k|^2$ is real of degree $\le2D-2$;
* $a\in V'$ iff $\deg W_a\le2D-2-L$;
* $\{W_a:a\in V'\}$ is exactly the space of real $W$ with $\deg W\le2D-2-L$ and $\mathrm{Re}\,\mathrm{res}_{\rho_k}(W/N)=0$
  for all $k$.

*Proof.* $\Psi_a(t)=-\sum_mt^{-m}m^{-1}\sum_ka_k\mathrm{Im}\rho_k^m$ for $t\to+\infty$. Vanishing of the first $L$
moments is equivalent to $\Psi_a'=O(t^{-L-2})$, i.e. to $\deg W_a\le 2D-2-L$. Conversely, a real $W$ with
$\deg W\le2D-2$ has $W/N=\sum_k(c_k/(t-\rho_k)+\bar c_k/(t-\bar\rho_k))$, and the residue $c_k$ is $a_k/(2i)\cdot$(real)
iff $\mathrm{Re}\,c_k=0$. ∎

*Correction to the grouping claim.*
* Take independent differences $a^{(1)},\dots,a^{(r)}$. Group roots by the column type $v_k=(a^{(i)}_k)_i$, and put
  $n_v=\prod_{v_k=v}(t-\rho_k)$.
* The correct near-identities are $\prod_v(n_v/\bar n_v)^{v_i}=1+O(t^{-L-1})$ (from Lemma 1.2(1)).
* With $N_v=n_v\bar n_v$ in place of $n_v/\bar n_v$ the statement is **false**. On the cube cluster with
  $y=(-\frac94,-\frac23,\frac14,\frac83)$ and the difference $G=\{r_0,\dots,r_3\}$, one gets
  $N_+/N_-=1+\frac{337}{36}t^{-1}+\dots$ (`prove_checks/s3_check.py`, exact). The Blaschke version gives
  $1-\frac{8645 i}{324}t^{-4}$ for $\{\pm r_0,\pm r_1\}$, as it should.

**(S4) (assessment; the inequality itself is a conjecture shape).**
* Setup. Let $X$ be a smooth projective toric compactification of $T_\Lambda\cong\mathbb G_m^r$ whose fan contains the
  cocharacters $\pm\nu_k$. Let $\tilde X=\mathrm{Bl}_eX$, with exceptional divisor $E$ and boundary $\partial$. Then
  $K_{\tilde X}+\partial=\pi^*(K_X+\partial_X)+(r-1)E=(r-1)E$.
* For the cluster curve $f$ ($f\equiv e\bmod s^m$, $m=L+1$) the strict transform has $\tilde f\cdot E=m$: the order of
  the maximal ideal of $e$ along $f$ is $\min_i\mathrm{ord}_\infty(G_i-1)=m$. Also
  $N^{(1)}_\partial(f)\le|S|\le2D$.
* So the truncated function-field Vojta inequality
  $h_{K+\partial}(f)\le N^{(1)}_\partial(f)+(2g-2)+\varepsilon h(f)+O(1)$ reads $(r-1)m\le|S|-2+\varepsilon h$, as (S4)
  says.
* At $c=0$, with $m=\lceil D/2\rceil$ and $h\asymp D$, this allows $r=5$ once $\varepsilon\gtrsim 2/D$ and forbids
  $r\ge6$ for $\varepsilon<1/2$. **Arithmetic confirmed.**
* Two caveats.
  * $K+\partial=(r-1)E$ is effective but has Iitaka dimension 0. Vojta's conjecture still asserts the inequality
    outside a proper closed $Z_\varepsilon$, but there is no structural reason for the cluster curves to avoid
    $Z_\varepsilon$. The rank-4 cube curves satisfy $(r-1)m=3\cdot 4=12\le14$, so they need no exception.
  * No unconditional function-field theorem of this shape is known to me for this pair (memory-quote; see §4.4).

## 2. The real structure, pushed

**Prop 2.1 (minimal pairs are real fibres; proved).** If $w=L+1$, then $B=P+i\kappa$, and the $w$ roots of $z_i$ on $S$
are exactly $P^{-1}(-i\kappa)$. If $\kappa>0$, then $\lfloor w/2\rfloor$ of them lie in the upper half-plane; if
$\kappa<0$, then $\lceil w/2\rceil$. Conversely, every real monic $P$ of degree $w$ and every $\kappa\ne0$ such that
$P+i\kappa$ is squarefree give a minimal pair. The parameter count is $w+1$ modulo nothing, i.e. $w-1$ modulo real
affine maps, which matches $2w-L$.
*Proof.* Lemma 1.1. The number of roots of $P-c$ in the upper half-plane is constant for $c$ in a fixed open
half-plane, since roots cannot cross $\mathbb R$ where $P$ is real. For $c\to-i\infty$ the roots approach
$c^{1/w}\cdot\mu_w$, which gives the stated count. ∎

**Prop 2.2 (triple identity; proved).** Take three points of a cluster.
* *Blocks.* Up to global sign each coordinate is of type $A$ (all equal) or $O_l$ (point $l$ is the odd one out).
  Let $\beta_l=|O_l|$, and let $m_l$ be the product of the majority's factors over $O_l$.
* *Pairs.* Points $i,j$ differ exactly on $O_i\cup O_j$, with $B_{ij}=\bar m_im_j$. So
  $\omega_{ij}:=\mathrm{Im}(\bar m_im_j)=\mathrm{Re}\,m_i\,\mathrm{Im}\,m_j-\mathrm{Im}\,m_i\,\mathrm{Re}\,m_j$ is a real
  polynomial of degree $\le\beta_i+\beta_j-L-1$.
* *Plücker relation:* $\omega_{23}m_1+\omega_{31}m_2+\omega_{12}m_3=0$. This is the identity for three vectors in the
  rank-2 module $\mathbb R[t]^2$ with its symplectic form $\omega$.
* *Mutually minimal triple* ($c=0$, $4\mid D$, $\beta_l=D/4$). Each $\omega_{ij}=\kappa_{ij}$ is constant, so
  $$\kappa_{23}m_1+\kappa_{31}m_2+\kappa_{12}m_3=0,\qquad \kappa_{12}+\kappa_{23}+\kappa_{31}=0:$$
  the three majority polynomials, monic of degree $D/4$, lie in a **real pencil**. Moreover $\kappa_{ij}=c_i-c_j$
  for reals $c_i$, and $z_i-z_j=2i\kappa_{ij}\gcd(z_i,z_j)$.
* *Mason–Stothers* applied to the pencil relation is vacuous: the radical has degree $3D/4$ and the degree is $D/4$.

**Prop 2.3 (phase ordering for equidistant clusters; proved).** Let $c=0$, $4\mid D$, and let all pairwise distances be
$D/2$. Put $\Theta_i=\sum_ks_{ik}\phi_k$. Then for each pair $\Theta_i(t)-\Theta_j(t)$ has the constant sign
$\mathrm{sgn}\,\kappa_{ij}$ on all of $\mathbb R$, and $|\Theta_i-\Theta_j|<2\pi$. Hence the $M$ phase functions are
totally ordered, never cross, and all tend to 0 at $+\infty$ (at rate $|t|^{-D/2}$ for the gaps).
*Proof.* Here $\theta_{ij}=(\Theta_i-\Theta_j)/2=\arg B_{ij}$ and $\mathrm{Im}B_{ij}=|B_{ij}|\sin\theta_{ij}=\kappa_{ij}\neq0$.
Since $w=D/2$ is even, Lemma 1.2 gives $p=n$, so $\theta_{ij}\to0$ at both ends and never meets $\pi\mathbb Z$. ∎

**What the real structure does not do.** Every statement above is an exact consequence of the moment equations plus
$y_k>0$. The only one that is a priori an inequality on a single configuration is balance (Lemma 1.2). It does not
separate.

**Prop 2.5 (balanced fakes; proved).** Let $C\subseteq\{\pm1\}^{D_0}$ be any code with pairwise distance $\ge D_0/2$ and
difference span of dimension $r<D_0$. Put $D=2D_0$ and fix $c\ge0$ with $L=D_0-c-1\ge1$. Then there is a real
$L\times D$ matrix with distinct, pairwise non-opposite columns and positive first row, such that:
* $\{(s,-s):s\in C\}$ lies in one fibre;
* the pairwise distances are $\ge D/2\ge L+1$;
* every difference vector, and every integer vector of the span of the differences, is of the form $(a,-a)$, so
  $\sum=0$: balance holds exactly.

With a Theorem-B code at $D_0$ this gives $\ge D_0+2=D/2+2$ points.
*Proof.*
* Rows: $(x,x+e)$ with $x\in\mathbb R^{D_0}$ arbitrary and $e\perp\mathrm{span}(C-C)$. These annihilate
  $(\delta,-\delta)$, because $x\cdot\delta-(x+e)\cdot\delta=-e\cdot\delta=0$.
* There are $D_0-r\ge1$ independent choices of $e$. A generic choice of $x$ (first coordinate positive) and of small
  $e$ gives distinct, non-opposite columns.
* The rank of the resulting row space is $\ge L$, since the $x$ are free; take $L$ generic rows.
* Distances double: $2w\ge D_0$. ∎

So the real-structure inequalities that are available (balance, the zero counts of $\mathrm{Im}B$, Cauchy index) are
**satisfied by linear-size fakes**. The genuinely non-fake content, namely the pencil relations, the fibre structure and
the phase ordering, is algebraic and **equivalent to the moment equations**. It yields no inequality on $M$ that I
could prove.

## 3. Positive-dimensional loci

### 3a. The cube locus

**Theorem 3.1 (certified).** Let $y\in\mathbb R^4$ with $\sum y_k=0$, $\prod y_k=1$, and $|y_k|$ pairwise distinct
(nondegeneracy: then the 16 roots $\pm\varepsilon_k/y_k\pm iy_k$ are distinct). Let $N=Q\bar Q$ be the cube polynomial.
Then:
* the largest cluster of $N$, over **all 256 base points**, has exactly 8 points;
* every cluster of size 8 has rank 4;
* the cluster through $Q$ is exactly the 8 sets of ff Theorem C.

*Proof.*
* *Moments.* On the locus, $\mathrm{Im}(\pm r_k)^1=\pm y_k$, $\mathrm{Im}(\pm r_k)^2=2\varepsilon_k$ and
  $\mathrm{Im}(\pm r_k)^3=\pm g_k$ with $g_k=3/y_k-y_k^3$.
* *Coordinates.* A point is given by $(\alpha,\beta)$ with $\alpha_k=(s_k^+-s_k^-)/2$ and
  $\beta_k=(s_k^++s_k^-)/2$, exactly one of them nonzero. Its moment vector is
  $(\sum\alpha_ky_k,\ 2\sum\beta_k\varepsilon_k,\ \sum\alpha_kg_k)$.
* *Coincidences.* Two points in one fibre have $d=\alpha-\alpha'\in\{-2,\dots,2\}^4$ with
  $\sum d_ky_k=\sum d_kg_k=0$.
* *Elimination* (`cube_special3.py`, `cube_fibres3.py`; sympy). For each $d$ not in $\mathbb Z(1,1,1,1)$, a lex
  Gröbner basis of $\{\sum dy,\ \mathrm{numer}\sum dg,\ \prod y-1,\ u\prod_{i<j}(y_i^2-y_j^2)-1\}$ is computed.
  * Either it is $\{1\}$, or the solution set is finite.
  * Taken together, the real solutions form **48 points**, all with the same set
    $\{|y_k|\}=\{0.12379\ldots,0.91118\ldots,2.60974\ldots,3.39713\ldots\}$.
* *Generic members.* Off these points the only coincidences are $d\in\mathbb Z(1,1,1,1)$, so the fibres are
  combinatorial.
  * The fibre of $\alpha=0$ consists of the 6 choices $\beta\in\{\pm1\}^4$ with $\sum\beta\varepsilon=0$, together with
    $\alpha=\pm(1,1,1,1)$: 8 points, the Theorem C cluster.
  * Every other fibre has $\le4$ points.
* *The 48 special points.* The fibre sizes were computed at 50 digits after Newton refinement of the algebraic points;
  the profile is $(8,4,4,4,\dots)$ at every one of them. Values separated by more than $10^{-30}$ are certainly
  distinct, since the evaluation error is $<10^{-40}$, so the upper bound 8 is rigorous. The lower bound 8 is exact
  (Theorem C).
* *Rank.* An 8-point fibre is the cube code, of rank 4 (Prop 3.3 of `ff.md`). ∎

So on the cube surface the answer to (2a) is **yes**: $M=8$ and rank 4 for every nondegenerate member, with no
excess even at special members. The rigidity comes from the $m=2$ row being combinatorial ($2\varepsilon_k$). The
extra coincidences that would add rank (the $d$ with $\{0,\pm1\}$ entries and a valid $\beta$) have **no**
nondegenerate complex solution at all (`cube_special2.py`).

**Compositions** ($Q\circ h$) and **first-order deformations**: not proved.
* kill_ff's numerical test (three $h$ of degree 2, $D=16$: exactly 8 zero-sum sets) stands as evidence.
* *First-order deformations (evidence).* The numerical cube system has local dimension 2 (`ff.md` `zerosum.py`).
  This suggests that the locus is a component of the rank-4 incidence variety of its type, with first-order
  deformations only along the locus; I have not certified this.
* *Proved:* no rank-5 point lies on the locus itself.

### 3b. The horizontal-line locus (S5)

**Theorem 3.2 (proved).** Suppose all $\mathrm{Im}\rho_k=1$ (equivalently, by a real affine map, all
$\mathrm{Im}\rho_k$ are equal). Then:
* every cluster is a family of subsets $T_i\subseteq X=\{x_k\}$, all of the same size $\tau$, with equal power sums
  $p_1,\dots,p_{L-1}$;
* distinct cluster points have $|T_i\setminus T_j|\ge L$, i.e. Hamming distance $\ge2L$;
* $M\le2$ whenever $D>6c+6$. For $c=0$ and even $D\ge8$ this means rank $\le1$.

*Proof.*
* *Moments.* $\mathrm{Im}(x+i)^j=jx^{j-1}-\binom j3x^{j-3}+\dots$ is a polynomial in $x$ of degree $j-1$ with
  leading coefficient $j$. So the $L$ moment rows span $\{x^0,\dots,x^{L-1}\}$, and the condition is
  $\sum_k\delta_kx_k^j=0$ for $j<L$.
* *Pair distance.* $T_i\setminus T_j$ and $T_j\setminus T_i$ are disjoint sets of reals, of equal size, with equal
  $p_1,\dots,p_{L-1}$. By Newton their $e_1,\dots,e_{L-1}$ agree, so if both had size $\le L-1$ they would be equal,
  hence empty. Thus $|T_i\setminus T_j|\ge L$.
* *Range of $\tau$.* $L\le\tau\le D-L$.
* *Three sets.* Suppose $T_1,T_2,T_3$ are cluster points. Then $|T_i\cap T_j|\le\tau-L$, and Bonferroni gives
  $D\ge|T_1\cup T_2\cup T_3|\ge3\tau-3(\tau-L)=3L=3(D/2-c-1)$ for even $D$. This is impossible for $D>6c+6$. (For odd
  $D$, $L=(D-1)/2-c$, and the same argument works for $D>6c+3$.) ∎

The PTE chains of (S5) therefore have length $\le2$. Ideal PTE solutions from Chebyshev fibres give pairs, never
triples. The horizontal line is the *most degenerate* locus, not a rich one: it multiplies the minimal distance by 2.

### 3c. Clusters with all pairwise distances minimal (equidistant)

**Prop 3.3 (proved).** Let $c=0$ and $D$ even, and suppose all pairwise distances equal $L+1=D/2$. Then
$4\mid D$ (if $M\ge3$), $M\le D$, and $\mathrm{rank}=M-1$.
* $M\le D$, rank $=M-1$: in $\pm1$ form, distance $D/2$ means orthogonality, and orthogonal vectors are independent.
* $4\mid D$: three pairwise orthogonal $\pm1$ vectors need it.
* Structure: Props 2.1–2.3. Every pair is a real fibre, every triple a real pencil, and the phases are totally ordered.

Equidistant clusters are therefore the cleanest test of (R$_4$): a 6-point equidistant cluster has rank 5.

**Computations at $D=8$** (`equi_codes.py`, `equi_num.py`, `equi_inspect.py`; evidence).
* All orthogonal codes up to the permutations used: 4 classes for $M=5$, 3 classes for $M=6$.
* $M=5$ (rank 4):
  * Two classes are realised nondegenerately. The minimal root separation and $|\mathrm{Im}|$ are $0.034$–$0.061$
    at scale $\sum|\rho|^2=8$.
  * Newton in mpmath gives residual $4\cdot10^{-51}$.
  * The Jacobian rank gives local dimension 2, the expected value $2D-2-4L=2$.
  * Their full fibre profile is $(5,5,4,4,\dots)$. These are **new rank-4 families, different from the cube**: rank
    4 is not special to the cube.
  * Not certified: the solution set is positive-dimensional, so a certificate would need two slice equations and a
    Krawczyk test (not done).
* $M=6$ (rank 5): 150 starts per class, 450 in all, found only degenerate solutions (separation $<10^{-12}$).

Classification for general $D$: **open**. I could not turn the pencil relations of Prop 2.2 into a bound. The
obstruction is that the relations are linear in the $m_l$ with real coefficients, while the degeneracy constraints
(distinct, non-conjugate roots) are open conditions. Every argument I tried reduced to counting conditions, i.e. to
the dimension count.

## 4. The general statement: tools tried

**4.1 Wronskian of the $W$-space: fails, provably.**
* By Lemma 1.3 the real space $\{W_a:a\in V'\}$ is the *full* space of real $W$ with $\deg W\le2D-2-L$ and
  imaginary residues. It has dimension $\ge D-L$.
* Every subspace of it is realised by real vectors $a$. So any argument that uses only degrees, orders at $\infty$ and
  Wronskian vanishing gives at best $r\le\dim V'$, which is $\ge D/2+c+1$.
* Integrality of the residues is equivalent to $a\in\mathbb Z^D$. That is exactly the rationality of the subspace in
  (R$_4$), and a Wronskian does not see it: the Wronskian of $r$ such $W$'s has degree $\le r(2D-2-L)-r(r-1)$ and
  order $\ge rL$ at $\infty$. This gives only $r\le 2D-1-L$.

**4.2 Simultaneous Mason–Stothers on the near-identities: fails.**
* For the correct (Blaschke) identities $\prod_v(n_v/\bar n_v)^{v_i}=1+O(s^{L+1})$, each single one is
  $U-\bar U=2i\,\mathrm{Im}U$ with radical degree $\ge2w$. Mason is vacuous on it.
* Several of them together form an $S$-unit near-identity system with $|S|=2D$. The Brownawell–Masser / Voloch
  inequalities (memory-quote: height of a nondegenerate solution of a unit equation $\le C(n)(|S|+2g-2)$) give
  $h\le C|S|$. Since $h\asymp D\asymp|S|$, this is vacuous, as (S4) says.

**4.3 Deformation theory on the swept-out locus (sketch).**
* Let $B$ be an irreducible family of rank-$r$ configurations of a fixed type, and $Y$ the closure of
  $\bigcup f_b(\mathbb P^1)\subseteq T_\Lambda$, with $\dim Y=d$.
* Take a log resolution $(\tilde Y,\partial_Y)$. The free-curve argument of Theorem D, run on $\tilde Y$, gives
  $dm+1\le-K_{\tilde Y}\cdot\tilde C=\partial_Y\cdot\tilde C-(K_{\tilde Y}+\partial_Y)\cdot\tilde C$.
* For a closed subvariety of a torus, $\bar\kappa(Y^\circ)\ge0$ (Iitaka–Kawamata; memory-quote). So
  $K_{\tilde Y}+\partial_Y$ is pseudo-effective, and it is $\ge0$ on the movable class $\tilde C$.
* If moreover $\partial_Y\cdot\tilde C\le2D$, we get $d(L+1)+1\le2D$. That is $d\le3$ at $c=0$, an extension of
  Theorem D (the case $Y=T_\Lambda$). The cube saturates it: $d=3$.
* *Gap:* $\partial_Y\cdot\tilde C\le2D$ needs transversality of $\tilde C$ to $\partial_Y$ over the $2D$ points. It
  holds in $X$ (primitive cocharacters), but I have not proved it on a resolution.
* *Why it cannot bound $r$:* the bound is on $\dim Y$, not on $r$. A configuration with 0-dimensional moduli has
  $Y=C$, and the inequality is empty. Expected rank-$\ge5$ configurations have negative expected dimension, so they
  are exactly the rigid ones. **The method is blind to the case that matters.**

**4.4 Literature (WebSearch; arxiv.org could not be fetched, so abstracts and search snippets only).**
* *Lehmann–Riedl–Tanimoto, "Non-free curves on Fano varieties"* (arXiv 2304.05438; Osaka J. Math. 62 (2025) 19–50).
  For Fano $X$ (in char 0) there is a proper closed $V\subset X$ containing all curves of the non-dominant components
  of $\mathrm{Mor}(B,X)$.
* *Lehmann–Tanimoto, "Geometric Manin's conjecture and rational curves"* (arXiv 1702.08508), and
  *Lehmann–Sengupta–Tanimoto, "Geometric consistency of Manin's conjecture"* (arXiv 1805.10580). The a-invariant
  exceptional set is closed in characteristic 0. According to the survey snippets this uses boundedness of singular
  Fano varieties (Birkar).
* I could not confirm the Hacon–Jiang statement beyond its being cited in this program.
* **Applicability.** Our $\tilde X=\mathrm{Bl}_eX$ is not Fano: $-K_{\tilde X}=\partial-(r-1)E$, which is big but not
  nef for $r\ge2$. The relevant objects are morphisms with an $m$-th order contact condition at $e$ (a log/jet
  condition). These are not covered by the abstracts I saw. Even if a closed exceptional set were available, the
  conclusion would be that the cluster curves of large rank lie in a closed $V\subset\tilde X$. Rigid curves always
  lie in some closed set, so this gives no rank bound without a description of $V$.
  **Fails-because:** closedness of the exceptional set is the wrong kind of statement for single rigid curves.
* Function-field gcd theorems (Corvaja–Zannier, Brownawell–Masser, Voloch, Wang; memory-quote, not re-checked) have the
  shape $h\le C|S|$. They are vacuous here.

## 5. Partial statements (what is proved)

1. Lattice balance (Lemma 1.2) for every integer kernel vector. True, but satisfied by doubled fakes (Prop 2.5).
2. Horizontal-line locus: $M\le2$ (Thm 3.2).
3. Cube locus: $M=8$ exactly, rank 4, all base points, including the 48 special real members (Thm 3.1).
4. Equidistant: $M\le D$ and rank $=M-1$, with the pencil and phase-ordering structure (Props 2.2, 2.3, 3.3).
5. **Not obtained:** any bound $M(D,0)\le C D^{1-\eta}$, any bound on the number of minimal-distance neighbours of a
   point, or rank $\le r_0$ under a stated genericity hypothesis.
   * The generic statement "for $\rho$ outside a proper closed set, rank $=0$" is trivial: a generic configuration
     has no nonzero integer kernel vector. It is useless, because clusters live on the closed set.

## 6. Stop rule and recommendation

* Rank-5 at $c=0$, $D\le12$: none found by this lane. I searched only $D=8$ equidistant. Broader searches belong to
  the search lane.
* A rank bound was proved on two positive-dimensional loci: the cube surface (Thm 3.1) and the horizontal-line
  $D$-fold (Thm 3.2). Taken literally, the stop rule's second clause therefore does not fire.
* But neither proof generalises:
  * the cube proof is an elimination on a fixed 2-dimensional surface whose second moment row is combinatorial;
  * the horizontal-line proof uses a collapse of the moment rows to $x^j$, which forces distance $2L$.
* Every general tool fails for a stated reason:
  * $W$-space: blind to integrality;
  * Mason / gcd: error $\asymp|S|$;
  * deformations: blind to rigid configurations;
  * Manin-type closedness: wrong type of statement.
* The real structure reduces to the moment equations, apart from balance, and balance is faked.
* **Recommendation: close the function-field line for internal work.** The remaining content is exactly the
  excess-intersection question (R$_4$), and no tool in reach addresses rigid excess loci.

## Files
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round30/prove_checks/`:

| file | content |
|---|---|
| `cube_special.py`, `cube_special2.py` | base-$Q$ extra zero-sum sets on the cube locus |
| `cube_special3.py` | all base points: coincidences, Gröbner |
| `cube_fibres2.py`, `cube_fibres3.py` | max fibres at generic and at all 48 real special members |
| `cube_fibres3.out` | output of the special-member run |
| `s3_check.py` | the (S3) correction (exact) |
| `equi_codes.py`, `equi_num.py`, `equi_inspect.py` | equidistant codes and numerics at $D=8$ |

`cube_fibres.py` is a superseded first run with faulty hashing. Notes: `../prove_notes.md`.
