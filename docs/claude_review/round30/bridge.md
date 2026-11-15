# Round 30, lane bridge: function-field clusters and the integer problem

Lane bridge, round 30. Write-up: `round30/bridge.md`; scripts: `round30/bridge_checks/`. Labels: **proved**, **exact**, **numerical/evidence**, **conditional**, **open**.

## 0. Summary of results

**Status: obstruction plus calibration. No record is improved and no rank bound is proved. The growth rate
$(\frac23+\varepsilon)\log R/\log\log R$ is unchanged.**

1. **Template theorem (proved; Theorem 1.2).** A primitive Q(i)(t)-cluster with $D$ even at $c=0$ specialises at
   rational $t_0\to\infty$ (bounded denominators) to primitive integer clusters with
   $C(t_0)\to\Gamma/\sqrt{|g(t_0)|}$. Here $\Gamma=\max|\text{coeff}_{t^{D/2}}(z_i-z_j)|$ and $g$ is the fractional
   Gaussian gcd of the values. The best constant $C_\infty=\Gamma/\sqrt{\prod_\ell g_\ell^*}$ is a product of exactly
   computable local maxima. It is an invariant under affine changes and rescaling, it is never improved by
   composition, and it is attained along a progression. Clusters with $K<D/2$ (e.g. odd $D$) give $C\to0$; $K>D/2$
   gives nothing.
2. **Cube family as templates (exact).** The smallest constant found is $C_{\rm cube}=145877.658\ldots$, at
   $y=(-\frac94,-\frac23,\frac14,\frac83)$ with $t_0\in\frac12\mathbb Z$. This is $2.2\cdot10^4$ times better than
   the round-29 figure $3.18\cdot10^9$, and still $2.5\cdot10^4$ times worse than the golden $5.8977$. The other 27
   computed parameter classes (of 38 with $a,b$ of height $\le25$) are all worse, by a factor $\ge130$. Sub-clusters: 5, 6, 7 points at $30.9$, $433$,
   $7946$. An explicit 74-digit integer example is exactly checked. **No record improved.**
3. **Records versus function fields (exact codes).**
   * Conjugation symmetry forces the cube code (Lemma 3.2). The record 8-point cluster with $C=5.754$ has
     **exactly the cube code**, with unequal weights.
   * The golden family is a specialisation, but along units, of a level-$D/4$ function-field cluster. Its asymptotic
     code is a sub-code of the cube.
   * The 5.674 record (non-squarefree, rank 7) and the 10-point cluster (rank 7) are not specialisations of any
     known $c=0$ function-field cluster. The 10-point cluster would need a $\ge10$-point Q(i)(t)-cluster at $c=0$,
     which would itself be new.
4. **Question B (exact plus numerical).** All 56 five-subsets of a cube cluster have gcd 1, contact 4,
   $-K\cdot C=16$, so **$\alpha=4$ and $\delta=0$**. They are not A-points and impose nothing on $Z$.
   * They lie on none of the 15 quadric cones and none of the 15 cubic threefolds (exact, via $\sigma$-invariance),
     and off $W$ and $\iota W$ (numerical).
   * 24 of them are Zariski dense in the subtori $\{u_au_b=u_cu_d\}$; the other 32 sweep hypersurfaces.
   * $(L_5^\circ)$, $A_5$ and the reduction $A_5\Rightarrow$ uniform bound are untouched.
   * General dictionary: $\alpha<4$ iff the gcd-reduced 5-point cluster has negative level.
5. **Calibration (proved conditionally on the named conjecture).** The function-field analogue FF-$A_5$ follows from
   the truncated function-field Vojta conjecture for $(\mathrm{Bl}_eX,\pi^*\partial+E)$.
   * **(S4) is corrected:** with $E$ in the divisor the inequality is $r\mu\le|S|-1+\varepsilon h$. It forbids rank
     $\ge5$ at $c=0$ and allows rank 4; Theorem D is its $\varepsilon=0$ free-curve case. So (R$_4$) has exactly the
     Vojta shape.
   * FF-$A_5$ sees only 5-subsets with large gcd, which a cluster has once $M\ge9$ (then $\alpha\to15/4$). It does
     not see the rank.
   * A function-field proof would be weak evidence for Question B. A function-field counterexample over
     $\mathbb Q(i)$ would refute $A_5$ and $V_E(5)$ arithmetically.
6. **Lower bounds.**
   * The cube gives $M_\infty(C)\ge8$ for $C>145877.66$, against the known $C>5.89771$.
   * Unbounded Q(i)(t)-clusters would give $\sup_CM_\infty(C)=\infty$, which refutes $A_5(\delta)$ for all
     $\delta<1/30$ and $V_E(5)$, but not $(L_5^\circ)$. It refutes Cilleruelo–Granville only if the constants stay
     bounded, and nothing in the model controls them.

**Stop rule.** This lane found no rank-5 cluster and proved no rank bound on a positive-dimensional locus. By the
pre-registered rule this lane's contribution does not keep the function-field line open.

## 1. Templates: the specialisation theorem (S6) and the constant C

Notation. A *Q(i)(t)-cluster* is a set of distinct monic $z_1,\dots,z_M\in\mathbb Q(i)[t]$ of degree $D$ with
$z_j\bar z_j=N$ for all $j$ and $\deg(z_i-z_j)\le K=\lfloor D/2+c\rfloor$; it is *split* if all roots lie in
$\mathbb Q(i)$, *primitive* if $\gcd(z_1,\dots,z_M)=1$ in $\mathbb Q(i)[t]$. For real $t_0$ write $|z(t_0)|=N(t_0)^{1/2}$.
For rational $t_0$ let $g(t_0)\in\mathbb Q(i)^*$ generate the fractional ideal $\sum_j z_j(t_0)\mathbb Z[i]$ (defined up to
units; $\mathbb Z[i]$ is a PID). For a rational prime $\ell$ put $|g|_\ell=\prod_{\pi\mid\ell}N\pi^{v_\pi(g)/2}$, so
$|g|=\prod_\ell|g|_\ell$. For $i\ne j$ let $\gamma_{ij}$ be the coefficient of $t^{D/2}$ in $z_i-z_j$ ($D$ even).

**Lemma 1.1 (proved).** For a cluster with $D$ even and $c=0$, every $\gamma_{ij}$ is purely imaginary, and
$\gamma_{ij}\neq0$ iff $\deg(z_i-z_j)=D/2$.
*Proof.* Put $\delta=z_j-z_i$. From $z_j\bar z_j=z_i\bar z_i$: $z_i\bar\delta+\bar z_i\delta+\delta\bar\delta=0$. If
$\deg\delta=D/2$ the coefficient of $t^{3D/2}$ is $\bar\gamma+\gamma$ ($z_i$ monic, $\deg\delta\bar\delta=D$). ∎

**Theorem 1.2 (templates; proved).** Let $\{z_j\}$ be a primitive Q(i)(t)-cluster, $D$ even, $c=0$, and put
$\Gamma=\max_{i<j}|\gamma_{ij}|$.
1. For all rational $t_0$ outside a finite set, $w_j(t_0)=z_j(t_0)/g(t_0)$ are $M$ distinct Gaussian integers with
   $\gcd=1$ on the circle $|w|=R(t_0):=|z(t_0)|/|g(t_0)|$, inside an arc of length $\ell(t_0)$ with
   $$C(t_0):=\ell(t_0)/\sqrt{R(t_0)}=\bigl(\Gamma+O(|t_0|^{-1})\bigr)\,|g(t_0)|^{-1/2}$$
   along any sequence with bounded denominators and $|t_0|\to\infty$.
2. There is a finite set $P$ of primes ($\ell$ dividing a coefficient denominator of some $z_j$, or the norm of a
   difference of two roots; for non-split clusters, replace the latter by the resultants $\mathrm{Res}(z_i,\bar z_j)$
   and $\mathrm{Res}(z_i, z_j)$) such that $|g(t_0)|_\ell\le1$ for $\ell\notin P$, with equality when $t_0$ is
   $\ell$-integral. For $\ell\in P$ the function $t\mapsto|g(t)|_\ell$ on $\mathbb Q_\ell$ is bounded and attains its
   maximum $g_\ell^*$ on an open set.
3. Hence $C_\infty:=\Gamma/\sqrt{\prod_{\ell\in P}g^*_\ell}$ is the infimum of $\liminf C(t_0)$ over all
   sequences of rationals with bounded denominators tending to $\pm\infty$, and it is attained along an arithmetic
   progression $t_0=r+mn$, $n\to\infty$. In particular $M_\infty(C)\ge M$ for every $C>C_\infty$.
4. Sequences with unbounded denominators are worse: if $t_0=p/q$, $\gcd(p,q)=1$, then $|g(t_0)|\le B\,q^{-D}$ with $B$
   depending only on the cluster. The value $C(t_0)$ is unchanged under rational affine reparametrisation
   $t=as+b$, so $C_\infty$ is an invariant of the cluster modulo $\mathrm{Aff}(\mathbb Q)$ (and of the scaling
   $z_j\mapsto\lambda z_j$).

*Proof.* (1) $z_i\ne z_j$ as polynomials, so their values differ off a finite set; $|z_j(t_0)|^2=N(t_0)$ for real
$t_0$; dividing by a generator of the ideal makes the values Gaussian integers generating the unit ideal. For
$|t_0|$ large the points lie within an angle $O(|t_0|^{-D/2})$ of each other, so the shortest arc containing them
has length $2R\arcsin(\max|w_i-w_j|/2R)$, and $|z_i(t_0)-z_j(t_0)|=|\gamma_{ij}||t_0|^{D/2}(1+O(|t_0|^{-1}))$ when
$\gamma_{ij}\ne0$ (smaller otherwise), $|z(t_0)|=|t_0|^D(1+O(|t_0|^{-1}))$. Then
$\ell/\sqrt R=\max|z_i-z_j|/(|g|\sqrt{|z|/|g|})(1+o(1))$.
(2) Split case. If $\pi$ lies over $\ell\notin P$ and $t_0$ is $\pi$-integral, at most one root $\rho$ satisfies
$\pi\mid t_0-\rho$ (differences of roots are $\pi$-units); by primitivity some $z_j$ omits $\rho$, so
$v_\pi(z_j(t_0))=0$. If $v_\pi(t_0)<0$ then $v_\pi(z_j(t_0))=Dv_\pi(t_0)<0$ for all $j$. For $\ell\in P$ write
$f(t)=\sum_{\pi\mid\ell}\log N\pi\cdot\min_j\sum_{\rho\in z_j}v_\pi(t-\rho)/2$ for $t\in\mathbb Q_\ell$. Each
$t\mapsto v_\pi(t-\rho)$ is non-decreasing as $t$ approaches any fixed $\sigma\in\mathbb Q_\ell$ (it equals
$\min(v_\pi(t-\sigma),v_\pi(\sigma-\rho))$ once $t$ is close to $\sigma$), so every value of $f$ is dominated by a
value near one of the finitely many local images $\sigma$ of the roots; near $\sigma$, $f$ is eventually constant,
because the minimum is attained by points omitting the corresponding root. So $f$ is bounded and its maximum is
attained on a ball. (3) The local maxima at the finitely many $\ell\in P$ are attained on balls; by the Chinese
remainder theorem there is $r\in\mathbb Z[1/P]$ in all these balls, and $t_0=r+mn$ ($m$ a large product of powers of
$\ell\in P$, $n\in\mathbb Z$) stays in them, is integral off $P$, and tends to infinity; apply (1). (4) For $t_0=p/q$,
$q^Dz_j(p/q)=F_j(p,q)$ with $F_j$ binary forms (after clearing a fixed denominator) with no common zero on
$\mathbb P^1(\overline{\mathbb Q})$ (primitivity and monicity), so the gcd of their values at coprime $(p,q)$ divides a
fixed nonzero ideal. Affine invariance: $C(t_0)$ is defined from the values only. ∎

**Which function-field clusters give integer clusters.** (a) $D$ even, $c=0$: bounded constant $C_\infty$ (Theorem
1.2), nonzero iff some pair has $\deg(z_i-z_j)=D/2$. (b) $D$ odd at the model's "$c=0$" ($K=(D-1)/2<D/2$), or $D$ even
with $K<D/2$: $C(t_0)\to0$, so such a cluster of size $M$ over $\mathbb Q(i)$ gives $M(C)\ge M$ for **every**
$C>0$. (c) $K>D/2$: $C(t_0)\to\infty$; useless. Conversely, every polynomial family $z_j\in\mathbb Z[i][u]$ of
integer clusters with bounded $C$ is a Q(i)(t)-cluster with $K\le D/2$ (`growth/families.md` Lemma P). This is the
precise form of (S6). It covers *only* specialisation along a rational branch ($t_0\to\infty$ in $\mathbb Q$); the
record's best families are of a different type (§3).

*Consequence for odd $D$.* The record has $M(3,0)=M(5,0)=3$, and $M(7,0)\ge4$ over $\mathbb Q(i)$ (exact). Any
Q(i)(t)-cluster of size $\ge6$ with $D$ odd and $K=(D-1)/2$ would give $M(C)\ge6$ for **every** $C>0$, beating the
record's 5 (cyclotomic family). Moon's Plotkin bound gives $M\le D+1$ here, and ff Prop 3.3 at $\theta=4/7$ gives
$\le5$ at $D=7$. None is known.

## 2. The cube family as templates: the smallest $C$, and explicit integer examples

**Method (exact up to the stated search ranges).** By Theorem 1.2,
$C_\infty(y)=\Gamma(y)/\sqrt{\prod_{\ell\in P}g^*_\ell(y)}$. Each local maximum $g^*_\ell$ is computed exactly
(`localgcd.py`):
* for odd $\ell$, by the nearest-candidate argument in the proof of Theorem 1.2(2): the maximum is attained at a
  local image of a root (split $\ell$) or at $\mathrm{Re}\,\rho$ (inert $\ell$). Only finitely many candidates;
* for $\ell=2$ (ramified), by enumerating $t\in2^{-j}\mathbb Z/2^{k}$ with $k\le14$, together with the candidates
  $\mathrm{Re}\rho\pm\mathrm{Im}\rho$. A truncated enumeration can only *underestimate* $g_2^*$, so every stated
  $C_\infty$ is a valid upper bound, realised along an explicit progression.

A brute-force scan over 12 denominators $q\mid72$ (`template.py`) reproduces the optimum at $y_{\rm J}$.

**Parameters.** The rational points of the cube surface $\sum y=0$, $\prod y=1$ of height $\le25$ were enumerated
(`ratpts.py`, 76 points; after $y\mapsto-y$ and ordering, 38 classes). The visible one-parameter family
$\{\pm a,\pm1/a\}$ is degenerate for every $\varepsilon$-assignment (repeated roots; proved by inspection:
$1/a+ia$ and $-1/a-ia$, or $1/a\mp ia$, coincide up to sign or conjugation). The three $\varepsilon$-assignments
of a point gave the same $C_\infty$ in every case tested.

**Results** (`run_family.py`; 28 of the 38 classes completed within the time limit, outputs `fam25a.txt`, `fam25b.txt`, merged in `famall.txt`; the other 10 classes have larger heights and were not computed):

| $y$ | height | $\Gamma$ | local factors $g^*_\ell\ne1$ | $C_\infty$ |
|---|---|---|---|---|
| $(-\frac94,-\frac23,\frac14,\frac83)$ | 9 | 118.31 | $2{:}\,1.73\cdot10^{-4}$, $3{:}\,1.52\cdot10^{-4}$, $5{:}\,25$ | **145877.66** |
| $(-\frac{25}6,\frac23,\frac{18}5,-\frac1{10})$ | 25 | 527.9 | $2{:}\,2$, $3{:}\,1.5\cdot10^{-4}$, $5{:}\,3\cdot10^{-6}$ | $1.89\cdot10^7$ |
| 26 further classes (of the remaining 36) | up to 441 | | | between $1.9\cdot10^7$ and $1.4\cdot10^{12}$ |

The denominators of the roots, $r_k=\varepsilon_kd_k/n_k+in_k/d_k$ for $y_k=n_k/d_k$, cost a factor
$\prod_k(n_kd_k)^2=(\prod n_k)^4$ in $|g|$. The gcd gains at the primes dividing root differences recover only a
bounded part of it. That is why $C_\infty$ grows quickly with height. **Smallest value found:
$C_{\rm cube}=145877.658383\ldots$**, a factor $2.2\cdot10^4$ below the unreduced $3.18\cdot10^9$ of ff Cor C′. The
gain comes entirely from the gcd reduction and from sampling $t_0\in\frac12\mathbb Z$.

**Sub-clusters** (same $y$, all subsets, `example.txt`). The best $k$-point sub-cluster constants are
$C_\infty^{(5)}=30.88$ (subset $\{0,A^c,B^c,P,P^c\}$), $C_\infty^{(6)}=$ 432.87 (subset $\{0,1,B,B^c,P,P^c\}$), $C_\infty^{(7)}=$ 7946.5 (drop $A^c$). Subsets of
size $\le4$ that share a root have $C\to0$: their reduced level is negative (4 points with $C\to0$, already known).

**Explicit integer example (exact, `example.py`).** Take $y=(-\frac94,-\frac23,\frac14,\frac83)$ and $t_0=13999/2$,
in the optimal class. Dividing the values by their fractional Gaussian gcd gives 8 distinct Gaussian integers with
gcd 1 on $|w|^2=R^2$, where
$R^2=76715325307615615132388329160620427624067222429804201398681082902106198205$ (74 digits, squarefree, 19 primes,
all $\equiv1\bmod4$). The constant is $C=145877.668212\ldots$. The points are
$(-6193356331893854325752232243998600469,\ 6193356331893549048567755974257051838)$ and its image under
$(x,y)\mapsto(-y,-x)$, and likewise for
$(\dots785479872085784293288446,\ \dots617894447902436590554517)$,
$(\dots639625293685896839599381,\ \dots763749026302324556008638)$ and
$(\dots701743192057945485163746,\ \dots701631127930276532348117)$; the full digits are in `example.txt`. The
symmetry is the conjugation symmetry of Lemma 3.2, after a unit rotation. All norms are checked equal, distinctness
and primitivity are checked, and $C$ is computed to 60 digits. At $t_0=(10^{12}+3999)/2$ and
$(10^{30}+423)/2$ (200 and 488 digits) the constant is $145877.658383\ldots$, equal to $C_\infty$ to 9 digits.

## 3. The record's 8- and 10-point clusters: specialisations or sporadic?

"Specialisation of a c=0 function-field cluster" means $W=\{z_j(t_0)/g(t_0)\}$ for a primitive Q(i)(t)-cluster with
$D\ge2$, $c=0$, and some rational $t_0$.

**Lemma 3.1 (proved).** If $W$ is such a specialisation, then $|W|\le M$, the size of the function-field cluster.
For every Gaussian prime $\pi$ outside the finite set of Theorem 1.2(2), the $\pi$-column of $W$'s exponent code
equals a column of the function-field code, up to complement.
*Proof.* A good $\pi$ divides $t_0-\rho$ for at most one root $\rho$, and then the points containing $\rho$ (resp.
$\bar\rho$) have $v_\pi=v_\pi(t_0-\rho)$ (resp. $0$). Distinct points of $W$ come from distinct $z_j$. ∎

**Lemma 3.2 (proved; conjugation symmetry forces the cube code).** Let $W$ be an 8-point cluster on a squarefree
normalised $n$, closed under $w\mapsto\bar w$ with no point on the axis. Then every layer cut separates each of the 4
conjugate pairs. So the cuts are among the $2^4/2=8$ "transversal" cuts of the 4 pairs, and the code has rank $\le4$.
The cube code (ff Thm C) is the code with all 8 transversal cuts. The cube cluster is itself closed under
$z\mapsto\bar z$ (pairs $0/1$, $A/A^c$, $B/B^c$, $P/P^c$), so all of its specialisations are axis-type clusters in the
sense of `families.md` §3.
*Proof.* $v_\pi(\bar w)=v_{\bar\pi}(w)=1-v_\pi(w)$ when $p=\pi\bar\pi\,\|\,n$. In the cube, $\bar z_F=z_{F^c}$ because
conjugating every chosen root gives the complementary transversal. ∎

So the cube templates are axis clusters: $C^2\to8d_{\max}+4\lim y_0^2/a$ (`families.md` §3.1), with integer depth set
$\{d_j\}$. Their huge constants mean depth sets of size about $C^2/8\approx2.7\cdot10^9$. Corollary S
($C\ge4\sqrt2-o(1)$ for every symmetric 8-point family) applies to them.

**The record clusters** (`record_codes.py`, exact factorisation in $\mathbb Z[i]$; output in `record_codes.txt`):

| cluster | $n$ | $C$ | rank | cuts | structure |
|---|---|---|---|---|---|
| 8 pts, $a=590749322$ | $5\cdot13\cdot17\cdot29\cdot53\cdot181\cdot457\cdot733\cdot3389$ | 5.7542 | 4 | 8 | **exactly the cube code**: the 8 transversal cuts of the conjugate pairs, unequal weights $\log65,\log17,\dots,\log3389$ |
| 8 pts, $a=5743586$ (record) | $5^3 13^2 17^2\cdot89\cdot109\cdot557$ | 5.6739 | 7 | 10 | non-squarefree; layers of $5,13,17$ break the pair separation |
| golden, $a=23186$ / $7465178$ | $5^2\cdot17\cdot61\cdot89\cdot233$ / $5^2\cdot13\cdots1597$ | 5.8978 / 5.8977 | 5 | 6 | 4 heavy cuts (the blocks $J_{n+s}$) are transversal cuts; the two layers of 5 ($F$ merged with a block) are not |
| 10 pts | $5^3 13^2\cdot17\cdot29\cdot113$ | 8.2393 | 7 | 8 | no perfect matching separated by all cuts |

**Answers.**
1. *10 points: not a $c=0$ specialisation of any known function-field cluster.* By Lemma 3.1 it would need a
   Q(i)(t)-cluster of size $\ge10$ at $c=0$. None is known: $\sup_D$ of known values is 8, and $M(D,0)\le6$ for
   $D\le6$ and $\le8$ for $D=7$ (exact or Plotkin). Under (R$_4$) the LP would give $\le8$ for large $D$. So the
   10-point cluster is either sporadic or a specialisation of a function-field cluster that would itself be a new
   record ($M(D,0)\ge10$). The referee's diagonal-symmetric search had already found it unique up to $n\approx5\cdot10^{19}$.
2. *The golden family is a specialisation, but of a different function-field object.* Put $X=\varphi^n$. Then
   $J_{n+s}=a\varphi^sX+\tau(a)\varphi'^s(\pm1/X)$, so $Z_w(X)=X^{-4}z_w(X)$ with $z_w\in\mathbb Q(\sqrt5,i)[X]$ of
   degree 8. The roots come in pairs $\pm\rho_s$, the same $t\mapsto-t$ symmetry as the cube. The family is evaluated
   at the units $X=\varphi^n$, not at integers. Because of the factor $X^{-4}$, $|Z|\asymp X^{4}=X^{D/2}$, so bounded
   $C$ needs $\deg(z_w-z_{w'})\le3D/4$: **level $c=D/4$** ($L=1$ at $D=8$), not $c=0$. (Proved by the computation in
   `families.md` §2: angular spread $\asymp\varphi^{-2n}=X^{-2}$, so $\deg(z_w-z_{w'})=6$.) Its asymptotic code is the
   cube code restricted to the 4 transversal cuts "one $e$-side". In general, Pell/unit templates, which specialise
   along a pair of conjugate real-quadratic places, are function-field clusters at level $K=3D/4$ plus an
   involution. They are cheaper than polynomial templates by a factor 2 in the contact. This is why they hold the
   records.
3. *The 5.754 cluster has exactly the cube code (exact), with unequal weights.* Combinatorially it is a cube, as is
   every squarefree symmetric 8-point cluster of rank 4 that uses all 8 transversal cuts (Lemma 3.2). It cannot lie
   in the asymptotic range of a cube template from the explored parameters, since those have $C_\infty\ge1.46\cdot10^5$
   (§2). Whether it equals a cube value $z_F(t_0)/g$ at a small $t_0$ for some rational $y$ of large height was
   **not decided**. A finite $t_0$ has no bearing on $M_\infty$ in any case.
4. *The 5.674 record* has non-squarefree $n$ and rank 7, so it is not even combinatorially a cube. Small primes with
   multiplicity carry its structure. The table indicates it is **sporadic** with respect to the function-field
   picture; this is evidence, not a theorem.

**Bottom line for (1).** The polynomial ($c=0$, rational-branch) templates are the *weakest* mechanism. The record's
8-point families arise from unit templates (level $D/4$) or from axis surfaces (del Pezzo $X_{0134}$). Their
combinatorics is the cube code (Lemma 3.2), which explains why the cube is the function-field optimum's shape, but
the arithmetic realisation is different.

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
    points of $W$, and $0.2$–$0.4$ on random points. It puts the cube curves off $W$ in every $S_5$-labelling, reduced to 1120 by the $S_3$-symmetry of $W$: the smallest residual is $7\cdot10^{-3}$, at the first test point where each labelling fails
    (`wtest.py`; **numerical evidence**: 6 random starts per point, so a missed solution is possible).
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
   $g\ge D\,\phi(M)$, $\phi(M)=\bigl(\binom{\lfloor M/2\rfloor}5+\binom{\lceil M/2\rceil}5\bigr)/\binom M5$ (each column constant on that fraction of 5-subsets at least, by convexity), so
   $\alpha\le4(1-\phi(M))\cdot D/(D-2c)$. Here $\phi(M)>0$ iff $M\ge9$, and $\phi\to1/16$, giving
   $\alpha\to15/4$. This is the function-field form of `gauss.md` Lemma 1.2. So $\sup_DM(D,c)\ge9$ produces
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

## 6. Lower bounds for $M_\infty(C)$, and what unbounded function-field clusters would imply

**Proposition 6.1 (proved; Theorem 1.2).**
(a) A primitive Q(i)(t)-cluster of size $M$ with $D$ even, $c=0$ and constant $C_\infty$ gives $M_\infty(C)\ge M$ for
every $C>C_\infty$. Sequences realising it: $t_0=r+mn$, $n\to\infty$.
(b) A Q(i)(t)-cluster of size $M$ with $K<D/2$ (for instance odd $D$ at the model's $c=0$) gives $M_\infty(C)\ge M$
for every $C>0$.
(c) Composition $z\mapsto z\circ h$ ($h\in\mathbb Q[t]$) never lowers $C_\infty$: its values are a subset of the
original's values. Rescalings $z\mapsto\lambda z$ and rational affine changes of $t$ leave $C_\infty$ unchanged.

**The cube family alone.** $M_\infty(C)\ge8$ for $C>C_{\rm cube}$, where the best value found is
$C_{\rm cube}\le145877.66$ (at $y=(-\frac94,-\frac23,\frac14,\frac83)$, $t_0\in\frac12\mathbb Z$ in the optimal
class). This is an explicit upper bound, proved by the example in §2. That it is the infimum over the family is
**evidence only**: it is the minimum over the 28 computed classes of rational $y$, and $C_\infty$ grows with height. The
record's $M_\infty(C)\ge8$ for $C>5.89771$ (golden family) is better by a factor of about $2.5\cdot10^4$. **No record
is improved.**

**Unbounded function-field clusters.** Suppose $\sup_DM(D,0)=\infty$.
1. *Over $\mathbb C(t)$ only, i.e. with no $\mathbb Q(i)$-rational members:* no integer consequence (Theorem 1.2
   needs $\mathbb Q(i)$-coefficients). Function-field consequence: clusters with $M\ge9$ and $D$ large have
   5-subsets with $\alpha\le4(1-\phi(M))+o(1)\to15/4$ (§5). If these do not lie in one fixed $Z\subset T_5$, then
   FF-$A_5$, and with it the truncated function-field Vojta conjecture for $(\mathrm{Bl}_eX,\pi^*\partial+E)$, fails.
2. *Over $\mathbb Q(i)(t)$, with constants $C_M$:* $M_\infty(C_M)\ge M$ for all $M$, so $\sup_CM_\infty(C)=\infty$. By
   the absolute form of the A-branch (round 23 J4, refereed), this **refutes $A_5(\delta)$ for every
   $\delta\in(0,1/30)$, and refutes $V_E(5)$**. It does **not** refute $(L_5^\circ)$ (J4's scope remark), and it does
   **not** refute Cilleruelo–Granville, which fixes $C$, unless $\sup_MC_M<\infty$. The function-field picture gives
   no control of $C_M$: $C_M=\Gamma_M/\sqrt{g^*_M}$ is governed by coefficient heights and denominators, which nothing
   in the model bounds. In the only family known, $C$ grows with the height of the parameters (§2).
3. *Bounded function-field clusters ((R$_4$) true)* imply nothing for integers. Integer clusters need not come from
   Q(i)(t)-clusters at all; the records come from unit templates and sporadic points (§3).

So the function-field line can produce at most refutations of the Vojta-type routes ($A_5$, $V_E(5)$), never a
proof, and never by itself a refutation of the conjecture at fixed $C$.

## 7. Status of the coordinator's observations used, and of computations

* **(S1)** Correct; used in §4.1. For $w=D/2$ pairs of the cube, $z_i-z_j=A\cdot2i\,\mathrm{Im}B$ with
  $\mathrm{Im}B$ constant.
* **(S2)** Not used.
* **(S3)** Not used.
* **(S4)** **Corrected** (§5). The right divisor is $\pi^*\partial X+E$, which gives $r\mu\le|S|-1+\varepsilon h$.
  At $c=0$ this allows $r=4$ and forbids $r\ge5$ off $Z_\varepsilon$ ($\varepsilon<1/4$). So (R$_4$) has the Vojta
  shape, and the cube (rank 4, $16=|S|$) violates the $\varepsilon=0$ form by exactly 1, i.e. it is
  Vojta-exceptional. The statement about the gcd theorems is confirmed in order of magnitude (Corvaja–Zannier shape
  $(\deg u\deg v\chi)^{1/3}\asymp D$; seen through search results only).
* **(S5)** Not used.
* **(S6)** **Proved in corrected and sharpened form** (Theorem 1.2):
  - $C_\infty=\Gamma/\sqrt{\prod_\ell g_\ell^*}$ is an invariant of the cluster;
  - only $c=0$ with $D$ even gives a positive constant;
  - odd $D$, or $K<D/2$, gives $C\to0$;
  - "polynomial templates" are only the rational-branch case, and the record's families are not of this type (§3).

**Computations** (all in `bridge_checks/`):

| script | content | status |
|---|---|---|
| `glib.py` | Gaussian-rational polynomials, $\mathbb Z[i]$ gcd | library |
| `cube_basic.py` | cube structure; pair degrees; 5-subsets: $g$, $\mu_0$, $\alpha$, rank, profiles; 15 cones | exact |
| `template.py` | brute-force $C(t_0)$ scan over 12 denominators | exact values, mpmath output |
| `localgcd.py` | exact local maxima $g_\ell^*$, $C_\infty$ | exact (2-adic part: lower bound for $g_2^*$) |
| `ratpts.py`, `ratpts25.txt` | rational points of the cube surface | exact |
| `run_family.py`, `fam25*.txt` | $C_\infty$ over the parameter classes | as `localgcd.py` |
| `example.py`, `example.txt` | explicit integer clusters, sub-cluster constants | exact checks |
| `record_codes.py`, `.txt` | exponent codes of the record clusters | exact factorisations |
| `sweep_rank.py` | swept dimension 3 of each 5-subset family | numerical (40 digits, SVD) |
| `wtest.py`, `wtest.txt` | membership in $W$ | numerical evidence |

Notes: `bridge_notes.md`.

