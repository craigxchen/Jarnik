## 4. Question B: the 5-subsets of a cube cluster

**4.1 Dictionary (proved).** Let $X=\overline L{}^{\rm tw}_5\supset T_5$, $e$ the identity, $h=h_{-K_X}$. A rational
curve $f:\mathbb P^1\to X$ over $\mathbb Q$ meeting $e$ at a $\mathbb Q$-rational point (moved to $t=\infty$) is the
same thing as a 5-point Q(i)(t)-cluster $f(t)=[z_1(t):\dots:z_5(t)]$, $z_j\in\mathbb Q(i)[t]$ with common leading
coefficient and $z_j\bar z_j=N$ (the real structure $u\mapsto1/\bar u$ of $T_5$ fixes $f(t)$ for real $t$, because
$1/\bar z_j=z_j/N$). Remove the polynomial gcd $A$, $\deg A=g$, and put $D'=D-g$. Then:
* $-K_X\cdot C=2D'$ when $N/A\bar A$ is squarefree: each root of $N/A\bar A$ meets one toric divisor $E_B$
  transversally, and conjugate roots meet $E_B$ and $E_{B^c}$ (the parity lemma);
* the contact order with $e$ is $\mu_0=D-\max_{i<j}\deg(z_i-z_j)$, unchanged by removing $A$;
* along integer specialisations ($t_0\to\infty$ in $\mathbb Q$; a rational branch, McKinnon–Roth $r_q=1$) the
  approximation exponent is $\alpha_C=(-K\cdot C)/\mu_0=2D'/\mu_0$, which equals $4/(1+2\delta)$ for the A-point
  exponent $\delta$ (diameter $\asymp R''^{1/2-\delta}$, $R''$ the reduced radius).

For a 5-subset of a cluster at level $c$ ($D$ even) $\mu_0\ge D/2-c$, with equality in the generic case, so
$$\alpha_C=\frac{2(D-g)}{\mu_0}\le\frac{2(D-g)}{D/2-c},\qquad\text{and if }\mu_0=D/2-c:\quad\alpha_C<4\iff g>2c .$$
So **$\alpha<4$ (the A-point regime) means exactly that the gcd-reduced 5-point cluster has negative level
$c'=c-g/2<0$**. By Theorem 1.2(b) such curves specialise to 5-point integer clusters with $C(t_0)\to0$.

**4.2 The cube (exact, `cube_basic.py`, $y=(-\frac94,-\frac23,\frac14,\frac83)$; the combinatorial statements hold for
every non-degenerate $y$).**
* Every column of the cube code has weight 4 out of 8, so no 5-subset is constant on a column. Hence **$g=0$ for all
  56 five-subsets** (proved: a constant column would need 5 points inside a class of size 4).
* All pair degrees are 4, except $\deg(z_P-z_{P^c})=3$, so every 5-subset has $\mu_0=4$.
* Hence **$-K\cdot C=16$, $\mu_0=4$, $\alpha_C=4$ and $\delta=0$ for all 56 curves.** Along integer specialisations
  $C(t_0)\to C^{(S)}_\infty>1$ (§2), so $\delta(t_0)=-\log C(t_0)/\log R''\to0^-$: **no specialisation is an A-point
  for any $\delta>0$** once $R$ is large. This makes kill_ff's remark ("exponent exactly 1/2, no A-points") precise,
  now with the gcd reduction included.
* *Profiles.* The 5-subsets are very unbalanced. The 32 subsets of rank 4 meet 8 of the 15 cuts, with weight 2 each
  in degree units. The 24 subsets of rank 3 meet 4 cuts, with weight 4 each. So none of them is in the balanced scope
  of $(L_5^\circ)$.
* *Exceptional-set components (round 27).*
  - Quadric cones $(u_b-u_a)(u_c-u_a)=(u_d-u_a)(u_e-u_a)$: for all 56 subsets and all 15 cones the polynomial
    identity fails (exact; 0 hits out of 840).
  - Cubic threefolds $\iota\mathcal Q$: the curves are $\sigma$-invariant, where $\sigma=\overline{\phantom u}\circ\iota$.
    So $C\subset\iota\mathcal Q$ iff $C=\sigma C\subset\sigma\iota\mathcal Q=\bar{\mathcal Q}=\mathcal Q$. None.
  - $W$ and $\iota W$: the same argument reduces $\iota W$ to $W$, since $W$ is defined over $\mathbb Q$. A
    least-squares membership test against the parametrisation of round 23 Thm 2 was validated: residual $10^{-15}$ on
    points of $W$, and $0.2$–$0.4$ on random points. It puts the cube curves off $W$ in every $S_5$-labelling
    (`wtest.py`; **numerical evidence**, numbers in §7).
* *What they do lie on.* The 24 rank-3 subsets contain two complementary pairs $\{z,\bar z\},\{w,\bar w\}$, so they
  lie on the codimension-1 subtorus $S_{ab|cd}=\{u_au_b=u_cu_d\}$. That subtorus is defined over $\mathbb Q$, because
  Galois acts on characters by $-1$. Every 5-tuple containing two conjugate pairs lies on such a subtorus (proved:
  $\zeta\bar\zeta=w\bar w=R^2$). The swept dimension is 3 for every subset: the map $(y_0,y_1,t)\mapsto T_5$ has
  Jacobian rank 3 at random points, for all 56 types (`sweep_rank.py`; numerical). So the 32 rank-4 types sweep
  hypersurfaces $H_S$, and the 24 rank-3 types are **Zariski dense in the subtori $S_{ab|cd}$**.

**4.3 Verdict on part (2).** The cube 5-subsets are **not A-points**. They sit at exactly $\alpha=4=\dim X$, the value
Vojta predicts for $\alpha^{\rm ess}_e$, which is the counting-critical $\delta=0$; they do not go beyond the
threshold, since $A_5(\delta<1/30)\iff\alpha^{\rm ess}>15/4$ is a statement about $\alpha<4$. Consequences:
* nothing forces them into an exceptional set;
* they lie in none of the known split components: $\mathcal Q$ and $\iota\mathcal Q$ exactly, $W$ and $\iota W$
  numerically;
* their union is a finite union of hypersurfaces of $T_5$ (the subtori $S_{ab|cd}$ and the $H_S$), each swept by
  $\alpha=4$ curves.

**So nothing changes for $(L_5^\circ)$, for $A_5$, or for the reduction $A_5\Rightarrow$ uniform bound.** One piece of
structure is new: the locus swept by rational curves through $e$ with $\alpha=4$ contains whole subtori. Any proof of
$A_5$ by counting curves must therefore allow a Zariski-dense set of rational curves at exactly $\alpha=4$. This is
the function-field shadow of "$\delta=0$ is counting-critical".
