# Round 29, lane "clothes": the arc problem in other clothing

*Labels.* **proved** means a complete proof is given here; it is checked by the author only, so every new item is
**unrefereed**. **exact** means a finite exact computation (scripts in `clothes_checks/`). **evidence** means
data that proves nothing. **heuristic** means a count used by no proof. **open** means undecided.
**literature (snippet)** means a statement seen only in search-result snippets, not in the full text.
**memory** means quoted from memory and not checked. Record items are cited by their file in
`docs/claude_review/`.

## 0. Verdict

**The growth rate and the uniform bound are unchanged.** None of the six reformulations gives an
$o(\log R/\log\log R)$ statement, and none passes filter F1. What is new and correct:

1. **The residue-class twin (d) is a precise non-archimedean twin of the arc problem, and it is not
   easier** (§4).
   * *Clothing (proved, exact-checked).* For $q$ odd and coprime to $N$, the lattice points of $|z|^2=N$ in one
     class $z\equiv z_0\pmod q$ have three equivalent descriptions:
     * box points in a coset of the kernel of $a\mapsto\prod_p(\pi_p/\bar\pi_p)^{a_p}\in(\mathbb Z[i]/q)^*$;
     * the condition $q\mid a$ in the representation $2N=g(a^2+b^2)$ of form (a);
     * lattice points of a full circle of radius $\sqrt N/q$ about a centre of denominator $q$.
     Up to sign classes, the class count is the number of representations of $N$ by **one** positive binary
     quadratic form of discriminant $-4q^2$. At $q\asymp N^{1/4}$ it is the **Gaussian Lenstra problem at exactly
     the critical exponent $\alpha=1/4$**, where the Coppersmith–Howgrave-Graham–Nagaraj bound
     $O((\alpha-\frac14)^{-3/2})$ diverges.
   * *Same inequalities as the arc (proved).* The twin satisfies the arc's pair inequality, cut inequality and
     monomial-gap (Liouville) inequality, with $C\leftrightarrow2/\lambda$ for $q=\lambda N^{1/4}$. Theorem 1.1
     (`two_thirds/proof.md` Thm 7.5) transfers verbatim when $q$ has no prime factor $\le M$, so the twin obeys
     $M\le(\frac23+\varepsilon)\log R/\log\log R$.
   * *Smooth moduli (proved, a negative).* When $q$ is divisible by every prime $\le M$, which is possible
     once $M\lesssim\frac14\log N$, the pigeonhole engine of Theorem 1.1 is void. For such moduli the twin is
     *harder* to bound by known means.
   * *Ultrametric balls.* For $q=\ell^t$ the twin is the lattice points of the circle inside an $\ell$-adic ball.
     Ultrametricity makes the "arc" an exact subgroup, so $\ell$-adic analytic tools (Strassmann, Skolem) reduce
     to lattice counting in a coset and add nothing.
   * *Not easier (heuristic, labelled).* A Vojta-type prediction for the twin needs $S\ni q$ with $q$ varying, so
     no fixed-$S$ Vojta statement covers it.
   * *Boundedness: open.* Exact search gives small values: class size $\le5$ at $\lambda=1$, and $\le8$ at
     $\lambda=\frac12$, on 19 norms (**evidence**).
2. **Prime peeling (f) is decided negatively** (§5).
   * The coordinator's remark is correct (proved with an explicit constant, from Cilleruelo–Córdoba).
   * An elementary Plotkin form of peeling is proved: each peeled sub-arc is a code of relative distance
     $>\frac12$.
   * Together they give $M\le2W'/\log P_+(N')+O_C(1)$. That is the record's Theorem 3.10 (`round8/stability.md`)
     rediscovered, not a new bound.
   * Every Paley/Hadamard fake has $P_+\le M^{O(1)}$, and it satisfies the peeling inequality and the Plotkin
     inequality. The rough side of any smooth/rough dichotomy is therefore empty on fakes, and the smooth side is
     the original problem. **Peeling plus a dichotomy cannot beat the Hadamard profiles.**
3. **Forms (a), (b), (c), (e) are dictionary entries, not mechanisms** (§2). Each one is the pair condition
   plus the angle data.
   * The concentration theory for modular circles and hyperbolas reaches $H^{o(1)}$ only for boxes of side
     $H\le p^{1/4}$. That corresponds to arcs of length $R^{3/8}$, where Cilleruelo–Córdoba already gives $\le4$.
   * The divisor form is Erdős #887, which is open (record).
   * The strip form gives only the Jarník-type $O(C^{3/2}N^{1/8})$.

Filter outcomes for all six forms: **§6**.

## 1. Notation

$N=R^2$, $W=\log N$. After Lemma 1.1 of `two_thirds/proof.md` the cluster is normalised:
$N'=\prod_p p^{e_p}$ with all $p\equiv1\ (4)$, and
$z_x=\varepsilon_x\prod_p\pi_p^{a_{xp}}\bar\pi_p^{e_p-a_{xp}}$.
* $u_p:=\pi_p/\bar\pi_p$.
* For $v\in\mathbb Z^k$, $P_v:=\prod_{v_p>0}\pi_p^{v_p}\prod_{v_p<0}\bar\pi_p^{-v_p}$, so $|P_v|^2=e^{w(v)}$ with
  $w(v)=\sum|v_p|\log p$.
* The gcd of two points has norm $|\gcd(z_x,z_y)|^2=\prod_pp^{e_p-|a_{xp}-a_{yp}|}$, and $c_{xy}=(z_x-z_y)/\gcd$.
* The **pair condition** for an arc of length $\le C\sqrt R$ is
  $d_w(x,y):=\sum_p|a_{xp}-a_{yp}|\log p\ \ge\ W'/2-2\log C$, which follows from $1\le|c_{xy}|$ and
  $|z_x-z_y|\le CN'^{1/4}$.

## 2. The six forms, derived with all gcd and unit costs

**(a) $n(2N-n)=\square$** (**proved**; `dictionary.py`, 3792 pairs, 0 failures). Fix $z_1$. For $z\ne z_1$ on the
circle put $m=|z-z_1|^2/2$ and $z\bar z_1=u+iv$.
* $m\in\mathbb Z$, because $|z-z_1|^2=2(N-u)$.
* $u=N-m$ and $v^2=N^2-u^2=m(2N-m)$.
* With $g=\gcd(m,2N)$ one gets $m=ga^2$, $2N-m=gb^2$ and $2N=g(a^2+b^2)$.

Conversely, $(m,v)$ recovers $z=(N-m+iv)/\bar z_1$, which is a lattice point exactly when $\bar z_1\mid N-m+iv$. For
primitive $z_1$ this is the congruence $v\equiv\rho\,m\pmod N$, where $\rho=-x_1/y_1\bmod N$ is the root of
$\rho^2\equiv-1$ attached to $z_1$. So **cluster points correspond to $0<m\le C^2\sqrt N/2$ with $m(2N-m)$ a
square, plus one congruence mod $N$ that selects the sign class of $v$.**

The gcd cost is that $g$ carries the common Gaussian factor: $|z-z_1|^2=|\gcd|^2|c|^2$ with
$|\gcd|^2=N'\prod p^{-|\Delta a|}$. Smallness of $ga^2$ is therefore exactly the pair condition. The square
condition is pair data. All angle information sits in the congruence, that is, in the choice among the
$2^{\omega(N)}$ roots $\rho$. Form (a) is the pair inequality plus a root choice, and F1 applies, since fakes keep
the pair data and choose the angles.

**(b) Quadratic congruences / modular circles** (**proved** reduction; literature (snippet)).
* Put a box $B$ of side $H=CN^{1/4}$ around the arc. For $(x,y)\in B$, $|x^2+y^2-N|\le3CN^{3/4}$.
* So for any prime $p>6CN^{3/4}$, the cluster equals $\{(x,y)\in B: x^2+y^2\equiv N\ (p)\}$. Here
  $H\asymp p^{1/3}$, and no smaller $p$ helps, because spurious solutions only enlarge the count.
* Cilleruelo–Garaev (GAFA 2011, arXiv 1007.1526) give $M^{o(1)}$ for boxes of side $M<p^{1/4}$ on modular
  hyperbolas $xy\equiv\lambda$. Their general bound is $H^{4/3+o(1)}p^{-1/3}+H^{o(1)}$ (snippet). The circle
  analogue has not been checked (inference).
* At $H=p^{1/3}$ the general bound is $p^{1/9}$, which is useless. The $H^{o(1)}$ range $H\le p^{1/4}$ corresponds to
  arcs of length $\asymp R^{3/8}$, where Cilleruelo–Córdoba already gives $\le4$ points (the exponent
  $\frac12-\frac18=\frac38$ at $k=4$).

So nothing reaches the critical scale. The proofs are lattice-plus-divisor-bound arguments: pair data plus counting,
which fails F1.

**(c) Divisors in short intervals** (record). On the split form $xy=N$ this is Erdős–Rosenfeld (#887), which is
open in its absolute form. Its relation to Question B is in `round23/split.md` and `round4/hardness.md` Thm 4.
* The Gaussian version asks for Gaussian divisors $z\mid N$ with $|z|=\sqrt N$ in a sector of angle $CN^{-1/4}$.
  This is the problem itself.
* Cilleruelo–Jiménez-Urroz (JTNB 12 (2000), "The hyperbola $xy=N$") and Granville–Jiménez-Urroz give length
  bounds for intervals of $xy=N$ with $k$ lattice points (literature (snippet)). Nothing at the critical scale
  was found.

**(d) Rational centre, denominator $q\asymp\rho$.** See §4.

**(e) Norm form on a thin strip** (**proved**, elementary).
* The arc lies in a strip of length $L=CN^{1/4}$ and width $h\le L^2/(8R)+1\le C^2/8+1$.
* By Dirichlet, there is a primitive $v\in\mathbb Z^2$ with $|v|\le Q$ whose direction is within $1/(Q|v|)$ of the
  strip. The strip then meets at most $h|v|+L/Q+1$ lattice lines parallel to $v$.
* On each such line $t\mapsto|z_0+tv|^2$ is a quadratic polynomial, so the line holds at most 2 points.
* Taking $Q=\sqrt{L/h}$ gives $M\le2(2\sqrt{Lh}+1)=O(C^{3/2}N^{1/8})$.

The question "how often can one value be taken on the strip" is the problem itself. The Bombieri–Pila bound for
conics (memory) gives $L^{1/2+\varepsilon}$, which is no better. The only input is the 2-points-per-line count,
which fakes satisfy, so F1 fails.

**(f) Prime peeling.** See §5.

## 3. What the other communities know

| form | community result | reaches critical scale? | rate beyond $\log/\log\log$? | source status |
|---|---|---|---|---|
| (a), (c) | Erdős–Rosenfeld #887, Ruzsa #886, Cilleruelo–Jiménez-Urroz | no (#887 open) | no | record + snippet |
| (b) | Cilleruelo–Garaev, Bourgain–Garaev–Konyagin–Shparlinski, Chan–Shparlinski: concentration in boxes | only $H\le p^{1/4}$ (arcs $R^{3/8}$) | no | snippet |
| (d) | Lenstra 1984; Coppersmith–Howgrave-Graham–Nagaraj (Math. Comp. 77 (2008)): $O((\alpha-\frac14)^{-3/2})$ divisors in a residue class for $\alpha>\frac14$; Hales (Gaussian; record `round4/hardness.md`) | **diverges exactly at $\alpha=\frac14$** | no | snippet + record |
| (d), single forms | Blomer–Granville (Duke 135 (2006)): averages of $r_f(n)^\beta$ uniform in $D=o(x)$ | averaged only, so capped (`round10/analytic.md` Thm B) | no | snippet |
| (d), class groups | Ellenberg–Venkatesh: products of small split primes of norm $<\sqrt{\lvert D\rvert}/2$ are non-principal (the twin Liouville gap, §4.2) | gap only | no | memory + snippet |
| (e) | Jarník; Bombieri–Pila | $N^{1/8}$ | no | memory |
| arcs | Cilleruelo–Córdoba: an arc of length $\sqrt2R^{1/2-1/(4\lfloor k/2\rfloor+2)}$ has $\le k$ points; Cilleruelo–Granville, "Lattice points on circles, squares in APs and sumsets of squares" (CRM Proc. 43, 2007) | no | no | snippet (as quoted by Cilleruelo–Granville) |

No community has an $O(1)$ result at the critical scale or a rate beyond $\log/\log\log$. The divisor-in-residue-class
community states the criticality itself: its constant diverges at exactly the exponent the twin needs.

## 4. Deep dive 1: the residue-class twin (d)

Let $q$ be odd with $\gcd(q,N)=1$, and put $q=\lambda N^{1/4}$. Define
$T(N,q)=\max_{z_0}\#\{z:|z|^2=N,\ z\equiv z_0\ (q)\}$.

**Theorem 4.1 (dictionary; proved; `dictionary.py` 99 900 pairs and `ringclass.py` 25 cases, 0 failures).** For
$z\ne z_1$ on the circle, the following are equivalent:
* (i) $z\equiv z_1\pmod q$;
* (ii) $z\bar z_1\equiv N\pmod q$;
* (iii) $q\mid a$, where $m=|z-z_1|^2/2=ga^2$ as in (a);
* (iv) $\varepsilon\prod_pu_p^{a_p-a_{1p}}\equiv1$ in $(\mathbb Z[i]/q)^*$, where $\varepsilon=\varepsilon_z/\varepsilon_{z_1}$.

The class of $z_1$ is $\{z_1+q\delta\}$, with $\delta$ running over the lattice points of the circle
$|\delta+z_1/q|=\sqrt N/q$. The set $L=\mathbb Zz_1+q\mathbb Z[i]$ has covolume $q$, and
$\#\{z\in L:|z|^2=N\}$ is the number of representations of $N$ by one primitive positive form of discriminant
$-4q^2$. This count is the union of the classes $c z_1$ with $c^2\equiv1\ (q)$, so it lies between $T$ and
$2^{\omega(q)+1}T$.

*Proof.*
1. (i)⇔(ii): $z_1$ is invertible mod $q$, and $z_1\bar z_1=N$.
2. (ii)⇒(iii): Write $z\bar z_1=(N-m)+iv$. Then $q\mid m$ and $q\mid v$. Since $\gcd(q,2N-m)=\gcd(q,2N)=1$,
   $q^2\mid v^2=m(2N-m)$ gives $q^2\mid m=ga^2$. As $g\mid 2N$ is coprime to $q$, this gives $q\mid a$.
3. (iii)⇒(ii): $v=\pm gab$, so $q\mid v$, and $q^2\mid m$.
4. (i)⇔(iv): $z/z_1=\varepsilon P_{a-a_1}/\bar P_{a-a_1}=\varepsilon\prod u_p^{a_p-a_{1p}}$.
5. Lattice statement: $L$ is an $\mathcal O_q=\mathbb Z+q\mathbb Z[i]$-module. Its image in $\mathbb Z[i]/q$ has
   $q$ elements, so $L$ has index $q$, and the Gram form of a basis has discriminant $-4q^2$. A $z\in L$ has
   $z\equiv cz_1$, and $|z|^2=N\equiv c^2N$ forces $c^2\equiv1$.

$\square$

So the arc cluster, which is the box points in a **Bohr set** $\{\langle2a-e,\phi\rangle\in I\}$ of the circle group
with $|I|\asymp N^{-1/4}$, becomes the box points in a **coset of a subgroup** of index $\asymp q\asymp N^{1/4}$.
That subgroup is the kernel of the map to the ring class group $\mathrm{Pic}(\mathcal O_q)$. For $q=\ell^t$ the
pro-$\ell$ part is the condition $|\sum_p(a_p-a_{1p})\log_\ell u_p|_\ell\le\ell^{-t}$: an $\ell$-adic arc. It is an
exact subgroup because the metric is ultrametric. This is the task's "non-archimedean twin of the arc condition".

**Proposition 4.2 (twin Liouville gap; proved).** If $v\ne0$ (or $v=0$ with $\varepsilon\ne1$) and
$\varepsilon\prod u_p^{v_p}\equiv1\ (q)$, then $|P_v|\ge q/2$, that is, $w(v)\ge2\log(q/2)$.

*Proof.*
1. The congruence says $q\mid\varepsilon P_v-\bar P_v=:\beta$.
2. $\beta\ne0$. For $\varepsilon=1$, $P_v=\bar P_v$ forces $v=0$ by unique factorisation, since
   $\pi_p,\bar\pi_p$ are non-associate. For $\varepsilon\ne1$ and $v=0$ we have $\beta=\varepsilon-1\ne0$. For
   $\varepsilon\ne1$ and $v\ne0$, $P_v$ and $\bar P_v$ have different factorisations.
3. A rational $q$ dividing $\beta\ne0$ gives $|\beta|\ge q$, so $2|P_v|\ge q$.

$\square$

This is the arc's monomial gap with $e^{-w/2}$ replaced by $q^{-1}$. In class-group language it is the
Ellenberg–Venkatesh reduction lemma (memory): small-norm products of split primes are non-principal in
$\mathcal O_q$.

**Theorem 4.3 (transfer of the cut inequality and of Theorem 1.1; proved).** Let $z_1,\dots,z_M$ be distinct, with
$|z_j|^2=N$ and pairwise congruent mod $q$, where $q\ge\lambda N^{1/4}$, $q$ is odd and $\gcd(q,N)=1$. Normalise by
$g=\gcd(z_j)$, as in `two_thirds/proof.md` Lemma 1.1.
* (i) $\sum_{i<j}\log|c_{ij}/q|+\frac12Q\le\frac M4\log R'+\binom M2\log(2/\lambda)$. This is (2.5) of
  `two_thirds/proof.md` with $C=2/\lambda$.
* (ii) If moreover $q$ has no prime factor $\le M$, then $3M\log M\le W+(14+4\log^+(2/\lambda))M$. Hence the twin
  obeys $M\le(\frac23+\varepsilon)\log R/\log\log R$.

*Proof.*
1. $g$ is coprime to $q$, so the $z_j/g$ remain pairwise congruent mod $q$, and $q\ge\lambda N'^{1/4}$.
2. Identity (2.3) of the proof is purely algebraic: $\sum\log|z_i-z_j|=\frac{M(M-2)}4\log R'+\frac Q2+\sum\log|c_{ij}|$.
3. Twin bounds: $|z_i-z_j|\le2R'$, and $q\mid c_{ij}$, because $q\mid z_i-z_j$ and $\gcd(q,\gcd(z_i,z_j))=1$.
4. Hence $\sum\log|c_{ij}/q|+\frac Q2\le K\log2+\frac{M^2}4\log R'-K\log q$. Insert
   $\log q\ge\log\lambda+\frac12\log R'$ and $K\cdot\frac12=\frac{M(M-1)}4$, which leaves $\frac M4\log R'$. This
   proves (i).
5. For (ii): the lower bound (3.6) of the proof comes from collisions of the points modulo $\ell^a$ and $p^a$ for
   primes $\le M$, and from the levels. Since those primes are coprime to $q$, every forced factor divides
   $c_{ij}/q$. So (5.1) holds verbatim with $C=2/\lambda$, and Thm 7.5 and Cor 8.1 follow word for word. The
   chord bound was "the only geometric input" (proof.md §0), and step 4 replaces it.

$\square$

**Proposition 4.4 (smooth moduli defeat the engine; proved).** Suppose every prime $\ell\le M$ coprime to $N$
divides $q$. This is possible with $q\le\lambda N^{1/4}e^{(1+o(1))M}$, so with $q\asymp N^{1/4}$ whenever
$M\le(\frac14-o(1))\log N$.

Then every collision factor of (3.6) at those primes is already contained in $q$. The lower bound of Theorem 4.3(i)
keeps only the level terms $\frac12\sum_\tau(S_{p\tau}-M/2)^2$, which vanish on exactly balanced profiles. So the
method yields no bound at all on such clusters.

This is a statement about the method, not a construction. Whether smooth-modulus twin classes can be large is
**open**. Exact search over 19 norms found smooth-modulus classes of size $\le4$ (`twin_smooth.py`; **evidence**).

**F1 for the twin.** The twin's analogue of a fake keeps the primes and the profile and replaces the genuine
classes $[\mathfrak p_p]\in\mathrm{Pic}(\mathcal O_q)$ by arbitrary elements of an abelian group of order $\asymp q$.
* Every twin inequality used above holds for such *group fakes* as soon as the kernel contains the difference
  lattice $D$ of the profile and no nonzero kernel vector has $w<2\log(q/2)$. These are the pair, cut and gap
  inequalities, and the residue collisions at primes $\nmid q$.
* A random homomorphism of $\mathbb Z^k/D$ achieves the gap provided two things hold:
  * the saturation of $D$ has weighted minimum $\ge W/2-O(1)$;
  * the number $e^{O(k)}$ of short vectors is $<q$.
* **Heuristic, not checked:** for Paley profiles the record's real fakes satisfy the gap on $D$ only up to
  $2\log\lVert c\rVert_1$. So the transfer gives group fakes at $\lambda=M^{-O(1)}$, not at fixed $\lambda$.
* **Proved:** no input that sees only the box and the group structure can beat what group fakes satisfy. In
  particular, Strassmann- or Skolem-type $\ell$-adic analyticity reduces to lattice counting, since the
  $\ell$-adic arc is linear in $a$.

**F2 for the twin.** Lemma K1 transfers: dividing by the gcd scales the $q$-adic data as it scales the archimedean
data. The forced-divisibility (Coppersmith) margin at $\alpha=1/4$ is $0$ (`round4/hardness.md` Prop 6.1 and §5.3).
The divisor-in-residue-class literature's own constant $O((\alpha-\frac14)^{-3/2})$ diverges there.

**Is the twin easier? No** (heuristic). The arc problem is proximity to $e\in T_M(\mathbb Q)$ at the fixed place
$\infty$, which is what the record's $V_E(5)$ uses, with $S=\{\infty\}$. The twin is the same proximity at the
places above $q$, which vary with $N$. Every Vojta-type prediction therefore needs constants uniform in $S$, an
abc-strength uniformity. The arc problem needs none.

**Decision on (d).**
* Boundedness of $T(N,\lambda N^{1/4})$ is **open**. No proof, and no construction.
* It is **at least as hard** as the arc problem for every known method: the same inequalities hold, Theorem 1.1
  transfers, and on smooth moduli even Theorem 1.1 is lost.
* It is **equivalent** to uniform boundedness of $r_f(n)$ for single forms $f$ of discriminant $D=-4q^2$ with
  $n\le D^2/(16\lambda^4)$, up to the factor $2^{\omega(q)+1}$.
* Evidence (exact search, `twin.py`):

| $\lambda$ | max twin class, odd $q\ge\lambda N^{1/4}$ | max arc count at $C=2/\lambda$ |
|---|---|---|
| 1 | 5 | 3 |
| 1/2 | 8 | 5 |

  The norms are products of the first $2$–$9$ split primes, plus prime-power mixtures with $r(N)\le2048$. That the
  twin values exceed the arc values is expected: the twin maximises over about $N^{1/4}$ moduli as well as over
  base points.

## 5. Deep dive 2: peeling, Plotkin and the smooth/rough dichotomy (f)

**Proposition 5.1 (coordinator's remark, with constant; proved from Cilleruelo–Córdoba).** Let $p^e\Vert N'$ and
let $k$ be even with $p^e\ge(C/\sqrt2)^{2(k+1)/(k/2+1)}N'^{1/(k+2)}$. Then $M\le(e+1)k$.

*Proof.*
1. $z=\pi^a\bar\pi^{e-a}z'$ with $p\nmid z'$ and $|z'|^2=N_1=N'/p^e$. The map $z\mapsto z'$ is a similarity with
   ratio $p^{-e/2}$ and rotation $-(2a-e)\phi_p$.
2. So the points with a given $a$ lie on an arc of the same angular width $CR'^{-1/2}$ of the circle of radius
   $R_1=R'p^{-e/2}$.
3. With $k=2r$, Cilleruelo–Córdoba allows $\le k$ points if
   $R_1CR'^{-1/2}\le\sqrt2R_1^{r/(2r+1)}$, which is equivalent to $p^{e(r+1)}\ge(C/\sqrt2)^{2(2r+1)}R'$.

$\square$

**Proposition 5.2 (Plotkin peeling, elementary; proved).** Let $p^e\Vert N'$ with $e\log p>4\log C$. Then
$$M\le(e+1)\Big(1+\frac{W'-e\log p}{e\log p-4\log C}\Big).$$

*Proof.*
1. Points $x,y$ with the same $a_p$ have the pair condition $d_w(x,y)\ge W'/2-2\log C$, and the $p$-coordinate
   contributes $0$ to it.
2. Put $h_x(p',\tau)=\pm1$ according as $a_{xp'}\ge\tau$, over the levels of $N_1$, with weights $\log p'$. Then
   $\langle h_x,h_y\rangle_w=W_1-2d_w(x,y)\le-(e\log p-4\log C)=:-\delta<0$, where $W_1=W'-e\log p$.
3. $0\le|\sum_xh_x|^2\le mW_1-m(m-1)\delta$ gives $m\le1+W_1/\delta$ in each of the $e+1$ classes.

$\square$

**Corollary 5.3.** Take $p=P_+(N')$:
* $M\le 2+2W'/(\log P_+-4\log C)$ when $e=1$;
* $M\le1+e+2W'/(\log p-4\log C/e)$ in general.

This is the record's Theorem 3.10, $M<2W/(\log P_+-4\log^+C)$ (`round8/stability.md` Thm 3.1, attributed there to
a two-line consequence of Cilleruelo–Córdoba). **Nothing new about the rate.**

**Theorem 5.4 (the dichotomy cannot beat the Hadamard profiles; proved, given the record's fakes).**
1. Every recorded fake has $W\le(8+o(1))M\log M+8M\log(1/C)$ and $P_+\le M^{O(1)}$ (`round8/stability.md`:
   "every fake of size $O(M\log M)$ has $\log P_+=O(\log M)$").
2. The peeling inequalities 5.1–5.3 use only pair data on the sub-circles. Fakes keep the genuine primes and
   profile and satisfy every pair inequality, so they satisfy 5.1–5.3. For Paley fakes: $p\approx M^8$,
   $k\approx M$ primes, and 5.3 asks $M\le2+2M$.
3. Any dichotomy "large $P_+$ / smooth" therefore puts every fake in the smooth class.
4. On the smooth class the peeled sub-codes are exactly Plotkin-extremal (Hadamard rows at relative distance
   $\frac12$), and Plotkin's bound $M\le2\cdot\#\text{levels}$ is weaker than Theorem 1.1.
5. Recursive peeling over all primes of a smooth $N'$ only rewrites the cluster in exponent coordinates. The $e+1$
   sub-arcs in arithmetic progression of step $2\phi_p$ are a multi-arc angle condition, which fakes satisfy.

$\square$

The Plotkin code condition is the pair inequality. The profiles that are extremal for it are the ones the fakes
use. The remark "two cluster points differ on every prime set of product $>C^{-2}\sqrt N$" is F1-blind.

## 6. Filter table and conclusions

| form | F1 (fakes) | F2 (criticality) | outcome |
|---|---|---|---|
| (a) $n(2N-n)=\square$ | fails: pair data plus root choice | n/a (dictionary) | dictionary only |
| (b) modular circles | fails: lattice plus divisor bound | weaker than Cilleruelo–Córdoba | $H^{o(1)}$ only at arcs $R^{3/8}$ |
| (c) divisors | split form (#887, open); not a mechanism | record: Coppersmith margin $0$ at the vertex | open, record |
| (d) residue twin | group fakes satisfy every twin inequality (heuristic at fixed $\lambda$; proved for box-and-group-only inputs) | exactly critical: Lenstra/CHN constant diverges at $\alpha=1/4$; K1 transfers | **new clothing, not easier**; Theorem 1.1 transfers for $q$ coprime to small primes; void for smooth $q$; boundedness open |
| (e) strip | fails | n/a | $O(C^{3/2}N^{1/8})$ |
| (f) peeling | fails: pair-only; every fake is smooth | critical (CC/Liouville) | rediscovers Thm 3.10; **dichotomy decisively negative** |

**What to keep.** The single item worth an expert's time is Theorem 4.1 with Proposition 4.4: the twin is
representation numbers of single forms of discriminant $-4q^2$ at $n\asymp D^2$, and equally the Gaussian Lenstra
problem at exactly $\alpha=1/4$. Two questions follow, and both are open:
* Is $r_f(n)$ bounded for $n\le c\,D^2$, uniformly in $f$?
* Can smooth-conductor ring class groups carry large fibres, where the Theorem 1.1 engine is void?

A construction of large twin fibres would refute nothing about arcs, since no implication between the two is known.
A proof of twin boundedness would need an input outside every class mapped in the record, for the reasons in §4.

*Scripts* (`clothes_checks/`):
* `common.py`: point enumeration and arc maxima;
* `dictionary.py`: the identities of (a) and (d);
* `ringclass.py`: the form of discriminant $-4q^2$;
* `twin.py`, `twin_smooth.py`: twin counts.
