# Lattice points on short arcs: three Diophantine questions

*A short note for experts in Diophantine approximation, compiled on 2026-10-28 from the review record
`docs/claude_review/` of the repository `craigxchen/Jarnik` (branch `claude/upbeat-einstein-p0peju`). Each
statement below is the referee-corrected or killer-corrected version, and each carries a status label and a file
pointer. Paths are relative to `docs/claude_review/`.*

*Labels.* **Proved** means checked by an independent adversarial referee. **Proved (killer-checked)** applies to
round 20 only: the result was checked by an adversarial "killer" but not by a separate referee. **Proved by the
killer** marks the killer's own new results (round 20), which no second party has checked. **Conditional**
means proved from a named hypothesis. **Exact** means a finite exact computation. **Heuristic** means a count that
no proof uses. **Unrefereed** means checked by its author only. **Open** means neither proved nor refuted.

## 1. The problem and what is known

Let $N=R^2\in\mathbb N$ and $\mathcal E_R=\{z\in\mathbb Z[i]:|z|^2=N\}$. For $C>0$, let $M(C)$ be the supremum,
over all $R$ and over all arcs $I$ of $|z|=R$ of length $\le C\sqrt R$ **in any position**, of
$\#(I\cap\mathcal E_R)$.

> **Conjecture (Cilleruelo–Granville, critical case).** $M(C)<\infty$ for every $C>0$.

The scale $\sqrt R$ is critical because the sagitta $L^2/8R$ of an arc of length $L$ becomes of order $1$ there.
Arcs that touch an axis point are trivial, so the whole difficulty is uniformity in the position of the arc
(`overview/overview.md` §1.1).

**Known for every $R$.** The literature bound is $\ll_C\log R$ (Cilleruelo–Córdoba 1992; Bourgain–Rudnick).

> **Theorem 1.1** (**proved**, two referees; `two_thirds/proof.md` Thm 7.5, Cor 8.1). Let $M\ge2$ distinct lattice
> points lie on an arc of length $\le C\sqrt R$. Then $3M\log M\le\log R^2+(14+4\log^+C)\,M$. Hence for every
> $C,\varepsilon>0$ there is an explicit $R_0(C,\varepsilon)$ such that $M\le(\frac23+\varepsilon)\log R/\log\log R$
> for all $R\ge R_0$.

The proof compares an archimedean upper bound with an arithmetic lower bound for the product of the reduced chords.
After dividing by $\gcd(z_1,\dots,z_M)$ (the *normalised* cluster: $N$ odd, all $p\mid N$ split,
$z_j=\varepsilon_j\prod_p\pi_p^{a_{jp}}\bar\pi_p^{e_p-a_{jp}}$), put $c_{ij}=(z_i-z_j)/\gcd(z_i,z_j)$ (the reduced
chords) and $S_{p\tau}=\#\{i:a_{ip}<\tau\}$ for $1\le\tau\le e_p$. The link is an exact cut identity
(`overview/overview.md` §2.2, Step 2):
$\prod_{i<j}|z_i-z_j|^8=N^{M(M-2)}\prod_{p\mid N}p^{\sum_\tau(2S_{p\tau}-M)^2}\prod_{i<j}|c_{ij}|^8$.
The arithmetic side comes from pigeonhole collisions at primes $q\lesssim M$, together with a chain lemma on
exponent levels. The factor $\log M$ is the mass $\sum_{q\le M}\log q/q$. No bound $o(\log R/\log\log R)$ is
known for every $R$, either unconditionally or from any standard conjecture. The only standard-conjecture routes
give the uniform bound directly. They use Vojta's Main Conjecture for **one fixed** variety:
* $D=0$ on a resolution of $\mathrm{Bl}_YX_7$, where $X_7$ is the equal-norm variety and $Y$ its diagonal line
  (`growth/conditional.md` Thm 5.1);
* $D=0$ on $\overline M_{0,8}$ (§2 below);
* $D=E$, $S=\{\infty\}$ on a 4-fold (§3 below).

Explicit examples give $M(C)\ge10$ for $C\ge8.2394$ (one explicit 10-point cluster) and $M_\infty(C)\ge8$ for $C>5.89771$ (a
family) (`growth/families.md`, `growth/referee_families.md`).

**Why $\log R/\log\log R$ is a barrier for local inputs: fake circles** (`overview/overview.md` §§5.1, 5.4).
Divide a cluster by its Gaussian gcd. Then $N=\prod_jp_j^{e_j}$ with $p_j\equiv1\ (4)$, and
$z_x=\varepsilon_x\prod_j\pi_j^{a_{xj}}\bar\pi_j^{e_j-a_{xj}}$. The arc condition becomes a linear system in the
angles $\phi_j=\arg\pi_j$: for all $x$, $\langle2a_x-e,\phi\rangle\in\theta_0+[0,\Delta]+\frac\pi2\mathbb Z$, with
$\Delta\asymp R^{-1/2}$. A **fake circle** (angle model) keeps the genuine primes and the genuine exponent profile
but uses an arbitrary *real* $\phi$. That $\phi$ must solve the system and satisfy every monomial gap
$\mathrm{dist}(\langle v,\phi\rangle,\frac\pi2\mathbb Z)\ge\arcsin e^{-w(v)/2}$, where
$w(v)=\sum|v_j|\log p_j$; these gaps hold on actual circles. Since $W=\log N=cM\log M$ corresponds to
$M\approx(2/c)\log R/\log\log R$, an input class that admits fakes with $W=O(M\log M)$ cannot prove
$M=o(\log R/\log\log R)$. The record builds such fakes for several classes:
* Skew-Hadamard (Paley) profiles give fakes at $W\le(8+o(1))M\log M+8M\log(1/C)$. These satisfy every pair, cut,
  monomial-gap and rational-character inequality (`round8/lw.md` Thm 6.1; **proved**).
* Residue designs extend this to point-level residues at moduli $\le M^A$ with $A$ fixed in advance
  (`round5/construct.md`, `round5/lower.md`). Referee scope limits: all primes of $N$ must exceed the modulus bound
  $X$, and the twisted extension at $X=M^A$, $A>1$, is open (`round5/README.md`; `overview/overview.md` §5.2 item 2).
* Grafts cap robust averaged analytic inputs (`round10/analytic.md` Thm B; `round20/composite.md`).

Actual clusters satisfy $W\ge3M\log M-O(M)$. Not everything is capped. P3-all$^*$, the sharp-constant large
sieves and the low-height window $P_+(N')\le M^{2-\delta}$ remain uncapped; here $N'$ is the normalised norm and
$P_+(N')$ its largest prime factor (`round21/README.md`, status **open**),
and the residue caps are not uniform in $A$; see `overview/overview.md` §5.5. The record's synthesis is a judge's assessment, **unrefereed**
(`round17/README.md` §7.1; `overview/overview.md` §7.6). It says a new mechanism must be global, multi-place and
tied to the actual arithmetic of $N$. The only statements of that type in the record are Vojta-type, and their
sharpest arithmetic forms are Questions A and B below.

## 2. Question A: $(L_8)$, a statement about integer $2\times8$ matrices

**Data** (`growth/ptolemy.md` §7; `overview/overview.md` §7.1). The data are:
* a matrix $P=(P_0|\cdots|P_7)\in\mathbb Z^{2\times8}$ with primitive columns, $P_0=(1,0)^T$ and
  $P_j=(X_j,Y_j)^T$;
* positive integers $n_T$ for $\emptyset\ne T\subseteq[7]$, such that for all $0\le i<j\le7$
  $$\det(P_i,P_j)=t_{ij}\prod_{T\ni i,j}n_T,\qquad t_{ij}\in\mathbb Z\setminus\{0\}\quad(t_{0j}=Y_j);$$
* the **chain condition**: at each prime $p$, the sets $T$ with $p\mid n_T$ form two chains under inclusion, and
  the unions of the two chains are disjoint.

The *visible* blocks are those with $2\le|T|\le6$. There are $119$ of them, one for each boundary divisor
$D_T$ of $\overline M_{0,8}$. Singletons occur in no minor. The top block $[7]$ divides every minor with
$i,j\ge1$, but it can be inflated by a fresh prime ($z\mapsto nz$) without moving $\Phi$, so imposing balance
only on the 119 visible blocks is equivalent (`round18/README.md` §1, "Balance"). Column $0$ is the **anchor**: it
is normalised to $P_0=(1,0)^T$ and sent to $\infty$ in Kapranov's chart
$\mathbb P^5=\mathbb P(\mathbb Q^7/\mathbb Q\mathbf 1)$. The chain condition does not depend on the anchor, but $Z$
depends on the labelling, so a statement with an arbitrary anchor must symmetrise $Z$ (`round18/README.md` §1,
"Anchor"). Put
$H_\Sigma=\sum_j\log|Y_j|+\sum_{1\le i<j\le7}\log|t_{ij}|$. The datum is **$\eta$-balanced at scale $w$** if
$|\log n_T-w|\le\eta w$ for all 119 visible $T$. Finally, $\Phi\in M_{0,8}(\mathbb Q)$ is the point given by the
eight columns as points of $\mathbb P^1(\mathbb Q)$.

> **$(L_8)$** (**open**). There are $c>0$, $\eta_0>0$ and a proper Zariski-closed $Z\subset\overline M_{0,8}$ such
> that every such datum with $\Phi\notin Z$ and $\eta\le\eta_0$ satisfies $H_\Sigma\ge c\,w-O(1)$.

**Cleanest form** (`round18/theory.md` §1.3, Prop 2.1; `round18/census.md` Prop 1.3; **proved** reductions).
Every $\Phi\in M_{0,8}(\mathbb Q)$ has a primitive representative $x\in\mathbb Z^7$, taken modulo translation, with
point $0$ sent to $\infty$ and
$$|x_i-x_j|=\prod_{T\ \text{visible},\ T\ni i,j}d_T(\Phi).$$
Here $d_T=\prod_pp^{\ell_T(p)}$ is the **contact integer**, the finite part of the height $h_{D_T}(\Phi)$, and
$\ell_T(p)$ is the number of levels $t$ at which $T$ is a residue class of $x\bmod p^t$ (`round18/theory.md`
Lemma 2.2). Let $\mathrm{Bal}(\eta)$ be the set of $\Phi$ with $|\log d_T-w|\le\eta w$ for all 119 $T$. The
strongest natural form, $(L_8^c)^{\rm all}$, which implies $(L_8)$, says that **a rational point of $M_{0,8}$ off a
fixed $Z$ cannot have all 119 boundary contacts of the same size $e^{w}$, $w\to\infty$**; $(L_8)$ itself gives this
only at caterpillar points. Round 18 proved a sandwich, and no arrow is known to reverse:
$$(L_8^c)^{\rm all}\Rightarrow(L_8)\Rightarrow(L_8^f)\Rightarrow(L_8^c)^{\rm cat}\Rightarrow(L_8^{pc}).$$
* $(L_8^c)^{\rm all}(\eta_1)$: $\mathrm{Bal}(\eta_1)\cap\{w\ge w_1\}\subset Z$. It implies $(L_8)$ with
  $\eta_0=\eta_1/2$, $c=\eta_1/4$. Conversely $(L_8)$ implies $(L_8^c)(\eta_0,Z)$ only for $\Phi$ admissible from
  anchor 0 (any anchor, with symmetrised $Z$); this is not an equivalence (`round18/README.md` §2.1 item 1;
  `overview/overview.md` §7.1 item 2).
* "cat" restricts this to points that are **caterpillar** at every prime. A compatible split system is a
  *caterpillar* if one side of each split can be chosen so that the chosen sides form a chain; equivalently, all
  internal edges of its tree lie on one path. The chain condition is exactly the caterpillar condition on the block
  system at each prime (`round18/theory.md` Lemma 1.1; exact on all 39208 trees with 8 leaves). Caveat (referee
  item H): the block system equals the $p$-adic tree of $\Phi$ only at primes dividing no residue.
* "pc" restricts it to pairwise coprime $d_T$.
* $(L_8^f)$ replaces $cw$ by any $f(w)\to\infty$ with $f\le cw$. Every proof of $(L_8^f)$ must prove
  $(L_8^c)^{\rm cat}$, the $H_\Sigma=0$ content of $(L_8)$, so a slowly growing $f$ is no shortcut; but $(L_8^f)$ is
  not known to be equivalent to $(L_8)$ (`round18/referee_theory.md` §2.3, item A).

**Lattice form** (`round18/theory.md` Lemma 2.2, §5(2), as corrected by `round18/referee_theory.md` item G;
**open**; a consequence of $(L_8)$, not known to be equivalent to it).
* In Kapranov's model $\overline M_{0,8}\to\mathbb P^5$: $H\cdot\beta=31$, $E_T\cdot\beta=1$,
  $K=-6H+\sum_{3\le|T|\le6}(|T|-2)E_T$ and $K\cdot\beta=10$. Here $\beta$ is the balanced class, defined by
  $D_T\cdot\beta=1$ for all $T$.
* The congruence lattice $\Lambda(x)$ consists of the $y\in\mathbb Z^7$ with $y_i\equiv y_j\bmod p^t$ whenever
  $i,j$ lie in one residue class of $x\bmod p^t$ of size $\ge3$. It contains $x$ and has index $[\mathbb Z^7:\Lambda(x)]=\prod_{3\le|T|\le6}d_T^{|T|-1}$.
* A balanced point has $|x|\approx e^{31w}$. The Minkowski scale of its own lattice is $e^{49w}$, and the trivial
  (Liouville) bound is $e^{30w}$.

> **$(L_8^{pc})$, lattice form.** There are $\eta_1>0$, $w_1$ and a proper closed $Z$ with the following property.
> Let $e_T\in[e^{(1-\eta_1)w},e^{(1+\eta_1)w}]$, for $T\subseteq[7]$ with $3\le|T|\le6$, be pairwise coprime
> integers, and let $w\ge w_1$. Then the lattice
> $\Lambda_e=\{y\in\mathbb Z^7:e_T\mid y_i-y_j\ \text{for } i,j\in T\}$ has no primitive vector
> $y\notin\mathbb Z\mathbf 1$ such that all of the following hold:
> * $[y]\in M_{0,8}\setminus Z$;
> * $|y|\le e^{31w(1+\eta_1)}$;
> * the pair cofactors $(y_i-y_j)/\prod_{T\supseteq\{i,j\}}e_T$ are pairwise coprime of size $e^{w(1\pm\eta_1)}$.

In words, structured lattices have no primitive vector $e^{18w}$ below their Minkowski scale. Nothing between the
trivial bound and Minkowski is known for these lattices. The heuristic count gives $98-108=-10=-K\cdot\beta$ per
unit of $w$: $e^{98w}$ profiles against $e^{-108w}$ for shortness (**heuristic**). $(L_8^{pc})$ alone does **not**
give the uniform theorem; that needs the chain cases, that is, $(L_8)$ itself.

**Theorems around $(L_8)$.**
* **$(L_8)\Rightarrow$ uniform bound** (**proved**; `growth/ptolemy.md` §7 and Steps 0–5 of §6.2, refereed in
  `growth/referee_ptolemy.md`). The proof has four steps.
  1. Subdivide the arc, keep one Gaussian unit class, and divide by its gcd.
  2. A random balanced selection (item 532 of the notes on branch `agent/uniform-gcd-formalization`) shows that at least half of all ordered 8-tuples have
     balance $\eta\le2\sqrt2\cdot255/\sqrt{M-7}$ and $H_\Sigma\le3584\,w/(M-1)$, with $w=W/128$.
  3. At most $c(Z)M^7$ of the 8-tuples have $\Phi\in Z$ (Lemma 6.1).
  4. $(L_8)$ then bounds $M$ independently of $R$.
* **Vojta $\Rightarrow(L_8)$** (**conditional**; Thm D, `growth/ptolemy.md` §6). Assume Vojta's Main Conjecture with
  $D=0$ on $\overline M_{0,8}$ over $\mathbb Q$, with $B=\sum_TD_T$ the boundary and $A=K+B$, which is ample
  (Keel–McKernan). Then $(L_8)$ holds with $Z=Z_{1/43}$ and
  $\eta_0=1/381$. The key estimate is $h_K(\Phi)\ge10w-282\eta w-60H_\Sigma-O(1)$. The sign is
  $K\cdot\beta=2^{k-3}(k-8)+k+2$, which is $-7$ at $k=7$ and $+10$ at $k=8$. So $k=8$ is the first case in which the
  $D=0$ inequality says anything at balanced points.
* **Smaller $k$.** $(L_5)$ and $(L_6)$ are false (**proved**; `growth/referee_ptolemy.md` §3.3, `threshold/L8.md`
  Prop C and Thm B). $(L_7)$ is open and
  expected false. $(L_7)\Rightarrow(L_8)$ (`round18/census.md` Lemma 7.7).
* **Correction B** (`round18/referee_theory.md` §2.7). A refutation of $(L_8)$ refutes Vojta at balanced points, but
  not conversely: a refuter must produce caterpillar points.

**Round-18 facts** (`round18/README.md` §§0–3; all **proved** and refereed unless labelled).
1. **No feasible computation** (`round18/census.md` Thm E). Each prime divides at most $5$ visible blocks, so a balanced datum
   needs at least 24 primes. At $\eta=0.1$ the class is empty for $e^w<59.2$, and all 21 minors
   $\det(P_i,P_j)$ with $i,j\ge1$ are $\ge10^{49.4}$. Block-sampling searches return $H_\Sigma\asymp w$ even at
   $k=5,6$, where bounded-$H$ points are known to be abundant (**exact**). So a failed search says nothing.
2. **No family mechanism.**
   * Density, box and sieve arguments produce only $\kappa\le-\sum\log e_T$, where $\kappa$ is the finite part of
     $h_K$ (Prop 2.3).
   * Bounded-degree covering families of genus-$\le1$ curves whose members *each* carry infinitely many balanced
     small-$H$ points are excluded (`round18/census.md` Prop 7.5). Correction M1: a dense supply with finitely many points on
     each member, for example on the lines of Kapranov's $\mathbb P^5$, is **not** excluded.
   * Class-$\beta_7$ lifts of the $(L_6)$ families are excluded (`round18/lift.md` Thm 4.2). Classes $c\beta_7$
     with $c\ge2$, and point-level lifts, are open.
   * Restricted theorems are known only on special loci with locus-dependent constants, or for a fixed prime set
     (Ridout: $H_\Sigma\ge8(1-\eta)(1-\varepsilon)w-C(S_0,\eta,\varepsilon)$ when all blocks are $S_0$-units;
     theory Thm 2.5). The constant $C$ depends on $S_0$, so this is useless here.
3. **Chains carry no information on pairwise-coprime profiles** (theory Prop 2.1). On $\mathrm{Bal}^{pc}$ each prime
   divides one block, and the $p$-adic tree is a single edge. So the chain condition, the anchoring and the residues
   are vacuous, and **every proof of $(L_8)$ must prove $(L_8^{pc})$**. The gloss "nothing but
   $K\cdot\beta=10$ survives" is not a theorem (item G).

> **Question A.** Prove $(L_8)$, or at least its lattice form $(L_8^{pc})$, which every proof must cover. A proof of
> $(L_8)$ would be an unconditional case of $D=0$ Vojta on the rational 5-fold $\overline M_{0,8}$, for a structured
> but unbounded family of points. Alternatively, find, for every $\eta>0$ and every $w_1$, a Zariski-dense set of
> $\eta$-balanced points with $w\ge w_1$ that are caterpillar up to $o(w)$ of residue. That would refute both $(L_8)$
> and Vojta (Correction B; `round18/referee_theory.md` §2.7, `round18/theory.md` §3.1).

## 3. Question B: $(L_5^\circ)$, five lattice points on a circle

Every theorem in this section is **proved (killer-checked)** or **conditional**, except where labelled otherwise; the
statements $A_5$, $(L_5^\circ)$ and its arc form are **open**. Sources: `round20/gauss.md` §§1–4,
`round20/kill_gauss.md` §§1.2, 6–7, and `overview/overview.md` §§4.3, 7.1A.

**Setting.**
* Let $\zeta_1,\dots,\zeta_k$ be a **primitive** tuple of distinct lattice points of $|\zeta|=R''$, meaning
  $\gcd_{\mathbb Z[i]}=1$.
* Its similarity class $[\zeta]=(\zeta_j/\zeta_1)_j$ is a point of the norm-one torus
  $T_k=(R^1_{\mathbb Q(i)/\mathbb Q}\mathbb G_m)^k/\text{diag}$, which has dimension $k-1$. $T_k$ is the open torus of
  the twisted Losev–Manin space $\overline L{}^{\rm tw}_k$ (the permutohedral fan in $\mathbb Z^k/\mathbb Z\mathbf 1$,
  with Galois acting by $-1$ on cocharacters). $[\zeta]$ is invariant under $\zeta\mapsto\gamma\zeta$,
  $\gamma\in\mathbb Q(i)^*$, which is why the exceptional set $Z$ must live in $T_k$. For $A\subsetneq[k]$ nonempty,
  $E_A$ is the toric divisor of the ray $\mathbf 1_A$.
* **Unit class.** Fix the Gaussian primes $\pi_p$ once. Each point factors as
  $\zeta=\varepsilon\prod_p\pi_p^{a(p)}\bar\pi_p^{e_p-a(p)}$ with a unit $\varepsilon\in\{\pm1,\pm i\}$; points lie in
  one **unit class** when their units agree. Within one class, $\zeta_a/\zeta_b=P/\bar P$ with $P\in\mathbb Z[i]$
  conjugate-primitive (`growth/ptolemy.md` §1, following item 532 on branch `agent/uniform-gcd-formalization`;
  `docs/claude_review/` has no separate formal definition).
* Each threshold layer $(p,t)$ of the exponent profile gives a cut $\{A,A^c\}$ of $[k]$. Let $n_A$ be the product
  of the $p$ over the layers with that cut. Then $h_{E_A+E_{A^c}}([\zeta])=\log n_A$ exactly, and
  $\sum_A\log n_A=\log R''^2$ (`gauss.md` Lemma 1.1).
* A short arc means archimedean proximity of $[\zeta]$ to the identity $e$.
* Put $d_{ab}=\sum_{A\text{ separates }a,b}\log n_A$, define $t_{ab}$ by $|\zeta_a-\zeta_b|=2|t_{ab}|R''e^{-d_{ab}/2}$,
  and put $H_\Sigma=\sum_{a<b}\log|t_{ab}|$. Here $t_{ab}$ is a nonzero integer when $\zeta_a,\zeta_b$ lie in one Gaussian
  unit class. Across unit classes it need not be: $|t_{ab}|\in\sqrt{\text{odd}/2}$ occurs (killer correction, minor 1).
  The record allows either repair: restrict to one unit class (as below), or take $t_{ab}:=|c_{ab}|/2$ real, with
  $c_{ab}=(\zeta_a-\zeta_b)/\gcd(\zeta_a,\zeta_b)$, and read the conclusion as $H_\Sigma\ge cw-O(1)$
  (`round20/kill_gauss.md` §7 item 1; `round20/README.md` §1.1 item 1).

> **$A_k(\delta)$.** There are a proper Zariski-closed $Z\subset T_k$ and an $R_0$ such that every primitive
> $k$-tuple of distinct lattice points of a circle of radius $R''\ge R_0$, with all mutual distances
> $\le R''^{1/2-\delta}$, has $[\zeta]\in Z$.
>
> **$(L_k^\circ)$.** There are $c,\eta_0>0$, $w_0$ and a proper Zariski-closed $Z\subset T_k$ with the following
> property. Every primitive $k$-tuple **within one unit class** with $[\zeta]\notin Z$, $w\ge w_0$ and
> $|\log n_A-w|\le\eta_0w$ for all $2^{k-1}-1$ cuts has $H_\Sigma\ge cw$.

At $k=5$ there is also a sufficient **arc form** (**open**; it implies $(L_5^\circ)$ by the argument of `gauss.md`
Prop 3.1, and the converse is not proved; `overview/overview.md` §7.1A). *Take five primitive lattice points of a circle of radius $R$, in
one unit class, whose Gaussian profile is $\eta_0$-balanced over all 15 cuts and whose class lies outside $Z$. Then,
once $R$ is large, they do not all lie within $R^{7/15+c}$ of each other.* The exponent $7/15$ arises because each
pair is separated by 8 of the 15 cuts. So $|t_{ab}|\ge1$ (valid within one unit class) gives diameter $\ge2R''^{7/15}$ at exact balance, and
this is exactly the Liouville bound. For a balanced $k$-tuple taken from a large cluster on a critical $\sqrt R$-arc, dividing by the tuple gcd moves it
to the "subcritical" exponent $\frac12-\frac1{2^k-2}$ in the reduced radius (`gauss.md` Lemma 1.2).

**Theorems.**
* **Reduction** (`gauss.md` Thm 3.4′, Prop 3.1). Let $k\ge3$ and $0<\delta<1/(2^k-2)$. Then
  $A_k(\delta)\Rightarrow(L_k^\circ)\Rightarrow$ the uniform bound. For $\delta\le2^{-k-1}$, only a convexity
  averaging of the tuple gcd is needed. The converse $(L_5^\circ)\Rightarrow A_5$ is **not** proved.
* **The 4-point version is false** (`gauss.md` Thm 4.1, plus the killer's projection corollary; exact Jacobian rank
  reproduced at 6 rational points).
  * Take Gaussian linear blocks $H_T(s)=q_T(s-r_T)$, with roots $r_T$ on an explicit rational 2-parameter family.
    These give balanced 4-tuples with bounded residues at exponent $3/7=\frac12-\frac1{14}$, dense in $T_4$.
  * Hence $A_4(\delta)$ is false for every $\delta<1/14$, and $(L_4^\circ)$ is false for every $c,\eta_0,Z$.
  * Likewise $A_3(\delta<1/6)$ and $(L_3^\circ)$ are false (the killer's own projection corollary,
    `kill_gauss.md` §3, `k3_projection.py`; not independently checked).
  * So $k=5$ is the first case where the reduction can work.
* **Vojta predicts it** (**conditional**; `gauss.md` Thm 4.2, re-derived by the killer).
  * Let $Y_k=\mathrm{Bl}_e\overline L{}^{\rm tw}_k$, so $K_{Y_k}=-\pi^*\sum_\rho E_\rho+(k-2)E$.
  * $V_E(k,\varepsilon)$ is Vojta's Main Conjecture for $(Y_k,D=E,S=\{\infty\})$: off a proper closed
    $Z_\varepsilon$,
    $\lambda_{E,\infty}(y)+h_{K_{Y_k}}(y)\le\varepsilon h_A(y)+C_\varepsilon$.
  * Here $A=\pi^*L-\varepsilon'E$ is ample with $L=\sum_\rho a_\rho E_\rho$, $a_\rho>0$, and $c_A=\max_\rho a_\rho$, so
    that $h_A\le c_AW'+O(1)$ with $W'=\log R''^2$ (`gauss.md` Thm 4.2 proof).
  * If $k\ge5$, $\kappa=(k-1)(1+2\delta)/4-1>0$ and $\varepsilon<\kappa/c_A$, then $V_E(k,\varepsilon)\Rightarrow A_k(\delta)$.
  * At $k=5$, $\kappa=2\delta$, so $V_E(5)$ with one $\varepsilon<(1/32)/c_A$ gives the uniform bound (via
    $A_5(1/64)$), and $V_E(5)$ predicts $(L_5^\circ)$.
  * The instance is non-vacuous (`gauss.md` Thm 4.2, Remark (3)). The $D=0$ form would need $k\ge6$.
  * The literature's instance (Vojta on $X$ with a normal-crossings divisor of multiplicity 4 through $e$) gives the
    same $\alpha^{\rm ess}\ge4$ (`papers_check/README.md` §1.3).
* **Covering families are excluded** (`gauss.md` Thm 4.3(ii); $k\ge5$, $\delta>0$). Put $B=\pi^*\sum_\rho E_\rho$, so
  $h_B=W'$. For each $d$, the curves $C$ with $B\cdot C\le d$ that carry infinitely many "A-points" lie in a proper
  closed set. Here an A-point is the class of a primitive tuple with diameter $\le R''^{1/2-\delta}$. (Thm 4.3(i) was
  narrowed by the killer, `kill_gauss.md` §1.2 item 4.) The density mechanism is excluded at $k=5$ only
  **heuristically** (count $e^{-w}$ per scale).
* **Lemma K1: every recorded method is exactly critical** (`kill_gauss.md` §6.1; **proved by the killer**, not
  independently checked). Let
  $F$ be bihomogeneous of bidegree $(a,b)$.
  * Then $F(z)=g^a\bar g^bF(\zeta)$, and the forced divisor scales in the same way.
  * The angular width is unchanged, and so are the pair data and the residues.
  * Consequently, every forced-divisibility, Liouville, pair, cut, character or single-form inequality is the same
    before and after division by the gcd. The move from $\frac12$ to $\frac7{15}$ is a change of units.
  * Auxiliary forms are capped at ratio exactly $1$ (`threshold/construct.md` Cor C3, with Thm B), and single-form
    Lang–Waldschmidt is exactly tight.
  * Only Vojta's $D=E$ term and the heuristic count see the reduction.

**Approximation reading** (`kill_gauss.md` §6.2). Measure height by $h_B=\log R''^2$. An A-point satisfies
$\mathrm{dist}(y,e)\le H_B(y)^{-(1/4+\delta/2)}$, approximating the **rational** point $e$ of the 4-fold $T_5$.
* $\delta=0$ is counting-critical (**heuristic**).
* Liouville excludes balanced points for $\delta>1/30$, and arbitrary points for $\delta>1/10$ (Cilleruelo–Córdoba).
* What is needed is that, for some $\delta<1/30$, **the Liouville exponent at $e$ is not attained by a Zariski-dense
  set of balanced points**.
* At $k=4$ this fails (Thm 4.1). On $\mathbb P^1$, Liouville is attained at every rational point.
* Roth- and Schmidt-type theorems apply but are weaker. McKinnon–Roth (Invent. Math. 200 (2015)) Thm 5.1 gives,
  off a proper closed set, $\alpha^{\rm ess}_x\ge\beta_x(L)$, where
  $\beta_x(L)=\int_0^{\tau}\mathrm{Vol}(\pi^*L-\gamma E)/\mathrm{Vol}(L)\,d\gamma<\tau_x(L)$ (ibid. p. 25). Their Thm 6.3
  gives $\alpha_x\ge\varepsilon_x/2$, strict at rational $x$. At $e$ (on $X'$), $\beta_e(-K)\in[2.31,2.61]$ (inference),
  below Liouville's $10/3$. As McKinnon–Roth note (ibid. p. 24), at a $k$-rational point Liouville $\alpha_x\ge\varepsilon_x$
  is "a stronger inequality than $\alpha_x\ge\frac12\varepsilon_x$" (`papers_check/README.md` §2.3).
* So this is the McKinnon–Roth regime. According to the killer, the only predictions beyond Liouville/Seshadri are
  conjectural: Vojta, and McKinnon's conjecture that best approximations lie on rational curves. The killer also
  knew of no unconditional theorem that beats Liouville at a rational point of such a variety on a Zariski-dense
  set. Round 22 checked this at snippet level (two independent searches plus an offline arXiv mirror): no
  unconditional theorem with $\alpha^{\rm ess}>\tau_x$ at a rational point was found for any variety
  (`round22/README.md` §2, J3). Seven full texts were then read (McKinnon–Roth 2015, 2016; the McKinnon–Satriano
  erratum; Lehmann–McKinnon–Satriano 2024; Huang ANT 2021, Acta Arith. 2019, Bull. SMF 2020). Every unconditional lower
  bound on $\alpha^{\rm ess}$ in them is Liouville with a section, at most $\tau_x$. Every exact value of $\alpha^{\rm ess}$
  they compute equals $\tau$, with $\tau\ge\dim$ (`papers_check/README.md`).

**In the approximation-constant literature's own terms** (round 22; `round22/README.md` §§0, 4; `round22/toric.md`,
`round22/lit.md` with their killer reports). The 4-fold is the twisted Losev–Manin space $X=\overline{L}_5^{\rm tw}$. It
is smooth and **weak Fano**, not Fano: $-K$ is nef and big but not ample, with polytope the $A_4$ root polytope, and
$-K=\phi^*(-K_{X'})$ for the small crepant contraction $\phi$ onto the twisted Gorenstein toric Fano $X'$, with
$X(\mathbb Q)=X'(\mathbb Q)=T_5(\mathbb Q)$ and $(-K)^4=70$ (**proved**). With $H=H_{-K}=R''^2$, an A-point at
$\delta$ has McKinnon–Roth exponent $\gamma=4/(1+2\delta)$, and:
* $\varepsilon_e(-K)=\alpha_e(-K)=2$ (**proved**), attained only on diagonal conics, which carry repeated points;
* $\tau_e(-K)=10/3$ (**proved**, round 23 `width.md` Thm 2.5, refereed). For every $m$, the maximal order at $e$ of a
  section of $m(-K)$ is $\lfloor10m/3\rfloor$, for the split and the twisted form, and it is attained by Vandermonde
  products. The proof slices the lattice points of $a\Delta-b\Delta$ coordinatewise, which gives the same family again,
  and peels the slices. Consequently every finite family of auxiliary Laurent polynomials over $\mathbb Q(i)$, used
  through Liouville with toric forced divisibility, gives for arbitrary profiles exactly the Cilleruelo–Córdoba exponent
  $\delta=1/10$ and never more; for balanced profiles it gives $\delta=1/30$. McKinnon–Roth's base-locus Liouville at $e$
  is therefore exactly Cilleruelo–Córdoba. (McKinnon–Roth, Eur. J. Math. 2 (2016), Thm 3.3. It needs $L$ ample, so it is
  applied to $-K_{X'}$ on $X'$. The values are unchanged.) The open window $\alpha^{\rm ess}\in(15/4,4)$ is untouched.
* $A_5(\delta<1/30)$ is exactly $\alpha^{\rm ess}_e(-K)>15/4$ off a proper closed set, and $(L_5^\circ)$ is the same
  over balanced sequences; Vojta predicts $\ge4=\dim X$.

> **Question B, literature form.** Is the identity $e$ **essentially bounded** on the smooth weak Fano 4-fold
> $X=\overline L{}^{\rm tw}_5$ (so that $X(\mathbb Q)=T_5(\mathbb Q)$ for the anisotropic norm-one torus $T_5$), in the sense
> of Lehmann–McKinnon–Satriano (arXiv 2403.02480, Def 1.3: $\alpha_{\rm ess}(e,-K_X+\varepsilon A)\ge\dim X$ for every
> ample $A$ and rational $\varepsilon>0$)? Since $-K_X$ is big, this is equivalent to $\alpha^{\rm ess}_e(-K)\ge4$
> (inference). The same constant is obtained on $X'$, since $\phi$ is an isomorphism near $e$. Less ambitiously: is
> $\alpha^{\rm ess}_e(-K)>15/4$, at least along sequences whose 15 toric heights are asymptotically equal? A yes gives
> the uniform Cilleruelo–Granville bound (gauss Thm 3.4'). Vojta's Main Conjecture implies it. McKinnon–Satriano's
> erratum (correcting their TAMS 374 (2021) Prop 1.3) has Prop 0.1: "Let X be a smooth projective variety over a number
> field k, and let P ∈ X(k) … If −K_X is big, then Vojta's Main Conjecture implies that X is canonically bounded at
> P." Lehmann–McKinnon–Satriano Prop 1.4 gives the same for essential boundedness on any smooth projective $X$.
> Neither assumes splitness. Both use Vojta on $X$ with a normal-crossings divisor of multiplicity $\dim X$ through
> $P$, a different instance from $V_E(5)$ with the same conclusion. In Lehmann–McKinnon–Satriano Thm 1.5 essential
> boundedness is a hypothesis. No unconditional proof was found for any variety where $\tau_P(-K)<\dim X$. Full texts
> were checked for seven papers; McKinnon–Satriano (TAMS) and Monahan–Satriano remain unread (`papers_check/README.md`).

> **Split calibration** (round 23 `split.md`, refereed). The same window occurs on the *split* Losev–Manin 4-fold at
> $e$. There a primitive positive 5-tuple $u$ has $H_{-K}=(\max u/\min u)\,\mathrm{lcm}(u)$ and
> $\mathrm{dist}\asymp(\max u-\min u)/\min u$. Split Question B asks whether five distinct integers with gcd 1 in
> $[U,U(1+\mathrm{lcm}^{-(1+2\delta)/4})]$, $\delta<1/30$, have their ratio class in a fixed proper closed set once the
> lcm is large. Equivalently, it concerns five factorisations of $N$ within hyperbolic angle $N^{-(1+2\delta)/4}$ on
> $xy=N$. It is **not known**: it implies Erdős Problem #887 in its absolute form ($\exists K\,\forall C$: eventually at
> most $K$ divisors of $n$ in $[X,X(1+Cn^{-1/4})]$, for every $X$), and Ruzsa's #886 for
> $\varepsilon>(15/64)(1+2\delta)$. Both are open (formal-conjectures, October 2026). The proof uses a random 5-subset,
> averaging over threshold layers, and Schwartz–Zippel. The 2026 five-divisor construction for #887
> ($n=(x^2-1)(x^2-4)(x^2-16)(x^2-49)$, $z^2=56x^2+65$) is a rational curve hitting each of the 15 cuts once, with
> $\alpha=15/4$ exactly. It lies in a 2-parameter family sweeping a hypersurface $W$, which every admissible split
> exceptional set must contain. Such curves do not exist on the twisted form: for a curve over $\mathbb Q$ there,
> $E_A\cdot C=E_{A^c}\cdot C$, so a balanced twisted curve at $15/4$ has $B\cdot C\ge30$ and $\mathrm{mult}_e\ge8$. A
> proof of Question B valid for every $\mathbb Q$-form of the 4-fold would therefore solve #887; a proof for $T_5$ alone
> may need the anisotropy. Split Question B is undecided. It would fail if infinitely many such hypersurfaces existed,
> or if absolute #887 failed.

> **Anisotropy** (round 24 `theory.md`, `search.md`, refereed). *Pointwise, it does nothing.* Let $\mathcal L$ be the
> arguments "archimedean Taylor upper bound $\ge$ product-formula lower bound" whose lower bound may use integrality,
> toric forced divisibility at every place of $\mathbb Q(i)$, point parity $v_{\bar\pi}=-v_\pi$ and reality of values,
> with auxiliary objects that are isotypic sections for $\sigma_X=\iota\circ\overline{\phantom{x}}$, polynomials in
> group words of A-points (including $\iota y$), or polynomials in the Hilbert-90 vectors ($x_j=P_j/\bar P_j$). Every
> such inequality is implied by toric Liouville for one section (or its $n$-point form at the generating points), so
> on balanced simple profiles with disjoint primes nothing in $\mathcal L$ beats $15/4$ ($\delta=1/30$), and nothing
> beats $\delta=1/10$ for arbitrary profiles, for every family with $\log H(F)+\sum d_i\log d_i=o(\sum d_iw_i)$. The
> reason is that $\sigma_X$ fixes $e$ and orders at $e$, so Galois averaging is lossless; the isotypic parts of
> $H^0(m(-K))$ have maximal orders $\lfloor10m/3\rfloor$ and $\lfloor10m/3\rfloor-1$ in both forms. Twist-sensitive
> congruence divisibility at inert places contributes $O_F(1)$. *Globally, it is visible only through curves.* By
> the parity lemma a balanced twisted curve at $15/4$ has $B\cdot C\ge30$. Such curves exist over $\mathbb R$: an
> explicit one-parameter real family $M$ (two conjugate branches of contact 4, each ray hit once; Jacobian rank 28,
> so a component). Only finitely many curves of its complex component are defined over $\mathbb Q$ (geometric
> monodromy $S_4$ of the place configuration, Faltings; ineffective; none of height $\le300$). The same finiteness
> holds for the split form, so this is not an anisotropy effect. Whether the class has other components is open.
> Blind numerical solvers all missed $M$, so failed searches in this class are not evidence of emptiness. Consequently
> a proof of Question B for $T_5$ that is not also a proof for the split form (hence of Erdős #887) must use a global
> twist-sensitive input; the only one known is the parity lemma.
>
> The disjoint-primes hypothesis is essential (round 25 `intra.md`, `referee_intra.md`, **proved**). The
> $\binom M5$ gcd-reduced 5-subtuples $y_S$ of one cluster are the projections $\mathrm{pr}_S(Z)$ of a single point
> $Z\in T_M(\mathbb Q)$, $T_M=(R^1_{\mathbb Q(i)/\mathbb Q}\mathbb G_m)^M/\mathrm{diag}$. Every multi-point
> Liouville/Dyson comparison among them, with toric forced divisibility across shared primes, is one-point toric
> Liouville on $T_M$ at $Z$. Conversely, a star partition $T_M\cong T_5^{(M-1)/4}$ identifies route 5's
> forced-divisibility forms with intra-cluster forms. Two results follow.
> * If $Z$ is balanced on all $2^{M-1}-1$ cuts of $[M]$, nothing in this class with a Taylor-type archimedean bound
>   beats $15/4$, uniformly in the auxiliary form.
> * For 5-wise uniform cut designs (e.g. extended BCH, $M\ge8189$) every 5-subtuple is exactly balanced with
>   $\mu<4(1-\eta)$, yet route 5's LP supercriticality forces $F(y)=0$.
>
> So shared primes do beat the one-point cap, but only by producing exact algebraic relations. Using those relations
> requires route 5's open zero estimate at the actual point. This gives nothing towards Question B.

> **Absolute form of the A-branch** (round 23 judge J4; independent referee `round23/referee_J4.md`, **proved**). If
> $A_5(\delta)$ holds for one $\delta\in(0,1/30)$, there is $K=K(\delta,Z)$ such that for every $C>0$ and every
> $R\ge R_1(C)=\max(C^{2/\kappa},(CR_0)^2)$, every arc of length $\le C\sqrt R$ holds at most $K$ lattice points; so
> $\sup_C M_\infty(C)<\infty$ (here $\kappa=1-(15/16+\eta)(1+2\delta)>0$, $\eta=(1/(1+2\delta)-15/16)/2$). The same
> argument gives it from $A_k(\delta)$, $\delta<1/(2^k-2)$, for every $k\ge3$: convexity averaging of the tuple gcd alone
> reaches the sharp range, with no balanced selection and no unit-class restriction. Scope: this is proved from the
> $A_5$ branch only, not from $(L_5^\circ)$ (whose route subdivides the arc and gives $K$ linear in $C$). So a
> construction with $M_\infty(C)\to\infty$ would refute $A_5(\delta)$ for all $\delta<1/30$, and $V_E(5)$, but not
> $(L_5^\circ)$. The absolute form is only formally stronger than $M(C)<\infty$, and it also follows from
> Cilleruelo–Granville's own $R^{1-\varepsilon}$ conjecture ($\varepsilon=1/4$).

> **Question B.** Prove $(L_5^\circ)$, or $A_5(\delta)$ for some $\delta<1/30$. Or refute $(L_5^\circ)$ by a
> Zariski-dense supply of balanced bounded-residue 5-tuples at exponent $\ge7/15$ (for a rational family, one
> dominating the 4-dimensional $T_5$); that would refute $V_E(5)$. Unlike $k=4$, a 2-parameter family plus its
> curve parameter has image of dimension $\le3$ and refutes nothing. By the killer's dimension count, one-variable
> polynomial-block families cannot be dense at $k=5$ (**heuristic**), and bounded-degree covering families are
> excluded (Thm 4.3(ii)) (`kill_gauss.md` §§4, 6.3).
> The two-branch twisted Pell ansatz ($B\cdot C=30$, contact 4, round 24) yields real curves at exactly $15/4$ but only
> finitely many over $\mathbb Q$ on the known component; it can only enlarge $Z$.

**Round 27 additions to Question B** (`round27/README.md`; all **proved** or **exact** unless labelled).
* *Split exceptional set.* Every admissible split exceptional set contains, besides $W$ and its image under the
  inversion $u_a\mapsto\mathrm{lcm}(u)/u_a$, the 15 quadric cones $(u_2-u_1)(u_3-u_1)=(u_4-u_1)(u_5-u_1)$ (and their
  $S_5$-images), swept by the second 2026 five-divisor construction ($C=7/2$, $\alpha=7/2$), and the 15 cubic threefolds
  obtained from them by inversion. Only $Q\cap\sigma Q$ can carry rational points in the twisted form.
* *A 2-dimensional calibration.* At an explicit rational point $P$ of an explicit degree-2 del Pezzo surface with
  $\mathrm{Pic}=\mathbb ZK$: $\tau_P(-K)=3/2<2$, the $\tau$-curve has a real cubic branch field and $\alpha=2$, and no
  $\mathbb Q$-curve of degree $\le6$ through $P$ has $\alpha<2$. So "$\alpha^{\rm ess}>\tau$ with $\tau<\dim$" is open
  already there. The surface has no balance, so it does not model $(L_5^\circ)$; and the Losev–Manin tower has no
  lower-dimensional balanced analogue, since $\alpha_{\rm bal}(k)=(2^{k-1}-1)/2^{k-3}<k-1$ iff $k\ge5$ (judge check,
  **exact**, unrefereed).
* *Pointer* (literature seen through search results only). $\mathrm{Bl}_eLM_5$ is a contraction of $\overline M_{0,7}$
  (Castravet–Tevelev, Duke 2015), so $\tau_e(-K)$, Theorem B and the section caps are statements about
  $\mathrm{Eff}(\overline M_{0,7})$; Question A lives on $\overline M_{0,8}$.

**Compared with $(L_8)$** (`round20/README.md` §§5, 7). $(L_5^\circ)$ is smaller: a 4-fold and 5 points, against a
5-fold and 8 points. It needs the archimedean divisor $D=E$ rather than $D=0$, and its heuristic margin is much
thinner: $\kappa=2\delta\approx1/15$ and $e^{-w}$ per scale, against $K\cdot\beta=10$ and $e^{-10w}$. Neither
statement is known to imply the other. Fakes violate $(L_5^\circ)$, whereas $(L_8)$ is fake-blind.

## 4. Question C (smaller): positive-gain non-unit pencils on $\overline M_{0,8}$

**Gain** (`round18/theory.md` §2.4; `round19/theory.md` §1). Let $X=\overline M_{0,8}$, let $\partial=\sum_TD_T$
(119 components, snc), and let $\beta$ be the balanced class. For a pencil (a non-constant rational function) $f$,
write $F_p=f^*(p)\sim L$ and $R_p=(F_p)_{\rm red}\cdot\beta$. Define
$$\mathrm{gain}(f)=L\cdot\beta-(R_0+R_1+R_\infty)=\sum_{Y\subset F_0+F_1+F_\infty}\tfrac{m_Y-3}{3}\,Y\cdot\beta .$$
When the base locus lies in codimension-2 boundary strata, this equals $10-K^{\log}_{\tilde X/\mathbb P^1}\cdot\tilde\beta$.

**The abc criterion** (`round18/theory.md` Prop 2.6; **conditional** on abc; implication **proved**).
* A *unit pencil* is $f=\lambda A/C$ with $\lambda\in\mathbb Q^*$, where $A,C$ are coprime products of differences
  $x_i-x_j$ in Kapranov's chart. Put $B=C-\lambda A$.
* Claim: if $\mathrm{gain}(f)>0$, then $(L_8^{pc})$ holds with $Z=\{B=0\}$.
* Mechanism: on $\mathrm{Bal}^{pc}$ we have $\gcd(A(x),C(x))=\prod d_T^{\mu_T}$, so
  $h(f(\Phi))=(L\cdot\beta)w+O(\eta w)$. Divide $(\lambda A(x),B(x),C(x))$ by its gcd to get a coprime abc triple.
  The log of the radical of that triple is $\le(F_{\rm red}\cdot\beta)w+O(\eta w)$.
* Scope: the refereed abc barrier covers only relations among products of differences (`round8/weakest.md`
  Prop 4.2), so this was a genuine gap in the barrier map.

**Round 19 closed the unit part** (`round19/theory.md`, refereed in `round19/referee_theory.md`; summary in
`round19/README.md` §1).

> **Theorem 2.4** (**proved**; $\beta$-weighted log Riemann–Hurwitz). For every non-constant rational function $f$ on
> $\overline M_{0,8}$,
> $$\sum_{p\in\mathbb P^1}\big(L\cdot\beta-R_p^{\rm int}\big)\ \le\ 2L\cdot\beta+\Lambda_{k(f)}\cdot\beta\ \le\ 2L\cdot\beta+31.$$
> * $R_p^{\rm int}$ is the $\beta$-degree of the reduced non-boundary part of $F_p$.
> * $k(f)$ is the least $|S|$ such that $f$ factors through a forgetful map $\pi_S$.
> * $\Lambda_k=\pi_k^*\psi_k$, with $\Lambda_k\cdot\beta=16,24,28,30,31$ for $k=4,\dots,8$.
>
> Equality holds for the cross-ratio and for $4z(1-z)$.
>
> **Corollary 2.5.** $\mathrm{gain}(f)\le31-n'(f)$, where $n'(f)=\#\{T: f|_{D_T}\text{ constant}\}$.

*Proof idea.* The form $df$ spans a line subsheaf of $\Omega^1_X(\log\partial)$, and its divisor is read off the
fibres. The forgetful filtration
$F_k=\pi_k^*\Omega^1_{\overline M_{0,k}}(\log\partial_k)$ has line-bundle quotients $\Lambda_k$. Every line
subsheaf of $F_k$ that is not in $F_{k-1}$ has $\beta$-degree $\le\Lambda_k\cdot\beta$, because $\beta$ is positive on
every effective divisor (Keel–McKernan, via `threshold/L8.md` Thm A(a)). This is not a curve argument. Restricting to covering curves of
class $d\beta$ only gives gain $\le(2g-2)/d$, a bound that is $\ge10$ and so cannot exclude positive gain.

> **Theorem 3.3, Corollary 3.4** (**proved**). For every non-constant unit $u$ of $M_{0,8}$ and every
> $\lambda\in\mathbb C^*$,
> $$\mathrm{gain}(\lambda u)\le\Lambda_{k(u)}\cdot\beta-N(u)-n_1\le-2^{8-k(u)}\le-1\qquad(\le-16\ \text{if}\ k(u)=4),$$
> where $n_p=\#\{T:D_T\subset\operatorname{supp}F_p\}$ counts the boundary divisors in the fibre over $p$.
> The same bound $\mathrm{gain}\le-1$ holds for every composite $h\circ u$.
>
> The key input is $N(u)=\#\{T:\mathrm{ord}_{D_T}u\ne0\}\ge32$, from multilinear Schwartz–Zippel on $\{0,1\}^7$.

So Prop 2.6 never applies. Exact computation agrees: gain $\le-32$, attained, for all unit pencils of degree
$D\le6$ and all $\lambda\in\mathbb Q^*$ (`round19/search.md`, reproduced by the referee; `round18/theory.md` §2.4).
Morphisms have gain $\le-32$ (Bruno–Mella, **read from the abstract only**).

**What is open** (`round19/theory.md` §4; `round19/referee_theory.md` m2). The remaining case is non-unit pencils in
$$P^\circ=\{f\ \text{regular on}\ X\setminus\textstyle\bigcup_{T\ne T'}(D_T\cap D_T')\}.$$
Two things are needed, and both are open.
* **(i) A gcd extension of Prop 2.6** to $P^\circ$. This means bounding the excess of $\gcd(A(x),C(x))$ over
  $\prod d_T^{\mu_T}$ at balanced pairwise-coprime points. Round 18 asserted, without proof, that for base locus in
  codimension-2 boundary strata the excess gcd vanishes on $\mathrm{Bal}^{pc}$ (`round18/theory.md` §2.4, after
  Prop 2.6). Round 19 §4.1 records this as plausible but **unproved** and gives only an unrefereed sketch of the
  obstruction. (For base locus meeting the interior, outside $P^\circ$, round 18 calls the gcd "a Vojta-type gcd
  problem in its own right".)
* **(ii) A pencil with gain $\ge1$.** Any such pencil must satisfy these proved necessary conditions:
  * $n'(f)\le31-2^{8-k(f)}$ and $k(f)\ge5$;
  * $f$ is not of the form $h\circ u$ with $u$ a unit;
  * all three fibres are non-reduced, with $\beta$-average multiplicities $\mu_p=L\cdot\beta/R_p$ satisfying
    $\sum_{p\in\{0,1,\infty\}}1/\mu_p<1$;
  * $\sum_{p\notin\{0,1,\infty\}}(L\cdot\beta-R_p^{\rm int})\le30-n(f)$, where $n(f)=n_0+n_1+n_\infty$ counts the
    boundary divisors in the three special fibres ($n(f)\le n'(f)$).

  **No example is known**, and Theorem 2.4 alone allows gain $\le31$. Campana-type statements see only gcd
  multiplicities, so they do not decide the question.

Even with both (i) and (ii), the result would be only an abc-conditional $(L_8^{pc})$. Reaching $(L_8)$ also
needs a nested-chain extension. There, for the chain $\{1,2\}\subset\{1,2,3\}$, the gcd exceeds Prop 2.6's value
by a full block (`round18/README.md` §3.3, J2; **unrefereed**). The judge adds that the radical also counts $p$
once, which might compensate, and that this is not analysed; so J2 shows a gcd excess, not that the extension fails.

> **Question C.**
> * (a) Does some $f\in P^\circ$ have $K^{\log}_{\tilde X/\mathbb P^1}\cdot\tilde\beta<10$, that is,
>   $\mathrm{gain}(f)>0$?
> * (b) Does Prop 2.6 extend to $P^\circ$? This would be a gcd inequality for two forms at balanced points.
> * (c) Is Theorem 2.4 known, or useful elsewhere? The record notes that it holds verbatim for any movable class
>   $\gamma$, with $31$ replaced by $\max_k\Lambda_k\cdot\gamma$, the maximum taken over orderings (unused there;
>   `round19/theory.md` §7 item 4, referee m7). Is it known, or sharp, in that generality, and on $\overline M_{0,n}$?

## 5. What any proof must avoid

Each item below is a barrier result in the record (`overview/overview.md` §§5.4, 7.6). Statuses are as given there.
1. **Local inputs.** These are pairs, gcds, cuts, monomial (Liouville) gaps, rational characters, and point-level
   residues at moduli $\le M^A$ with $A$ fixed in advance. Fakes at $W=O(M\log M)$ satisfy all of them (**proved**,
   with the scope limits of `overview/overview.md` §5.2 item 2).
   Caveat: caps uniform in $A$ are not proved, and P3-all$^*$ is uncapped.
2. **Averaged or robust analytic inputs.** This covers large sieves with non-sharp constants, prime sums or box
   counts with power-saving error, and "almost all" statements. Grafted fakes satisfy them (`round10/analytic.md`
   Thm B, **proved**; composite-supported version `round20/composite.md`, killer-checked). Exceptions inside the
   low-height window are listed in `round21/README.md`.
3. **Auxiliary polynomials and single-place Diophantine approximation.**
   * On $X_M$: $\tau^*(M)=2-1/\lceil M/2\rceil<2$, so Liouville at $\infty$, and Ru–Vojta with $S=\{\infty\}$, never
     reach exponent $\frac12$ (`threshold/upper.md`, **proved**).
   * On $\overline M_{0,k}$ for any $k$: $\beta$ lies in the interior of the movable cone (`threshold/L8.md` Thm A,
     **proved** as a formal barrier).
   * On the Gaussian model: the forms are exactly critical (Lemma K1 with `threshold/construct.md` Cor C3 and Thm B).
   * Roth- and Schmidt-type bounds (McKinnon–Roth) give $\alpha^{\rm ess}\ge\beta<\tau$ at a rational target: weaker
     than Liouville.
   * Multi-point (Thue/Siegel/Dyson) arguments on the 5-subtuples of one cluster are one-point toric Liouville on
     $T_M$ at the cluster point, which is the determinant method with forced divisibility ("route 5"). In the growth
     regime they force exact algebraic relations and stop at a zero estimate at the actual point. With lacunary degrees
     they reduce to one-point caps (round 25 `intra.md`, **proved**; the comparable-degree dichotomy is **heuristic**).
   * *Zero estimates do not supply route 5's missing lemma* (round 26 `ze.md`, `referee_ze.md`; **proved** at the
     degree level, given the shape of Philippon's 1986 zero lemma, with margins $\ge e^{1.7\cdot10^4}$ that make its
     exact constants irrelevant). Route 5's Liouville step certifies, for one form $P_c$, vanishing at $Z^M$ only for
     integer endomorphisms $M$ with $\Lambda_M(c)<\mathrm{ord}(c)-O(E)$; for the powers $Z^j$ this means
     $|j|<\mathrm{ord}(c)/\Lambda(c)\le1.07\sqrt n/\lambda$, by the spectral cap
     $\mathrm{ord}(c)\le(\sqrt n/\lambda)\Lambda(\mathrm{supp}\,c)$, and only $|j|\le1$–$4$ on regime boxes. A
     Philippon estimate on $\mathbb G_m^d$ needs progressions of length about $d$, so its conclusion is met by the
     trivial subgroup. The dual framing (one exponential polynomial on $\mathbb G_a^d\times\mathbb G_m$ vanishing on
     the slab box) has obstruction subgroups whose counts do not depend on the point, is vacuous at the degree level,
     and is false in Hilbert-sharp form. Moreover it concerns $\mathrm{Bs}(V(A,m))$, while route 5 needs
     $Z\notin\mathrm{Bs}(V_H)$ for the height-cut forced space $V_H$. The global statement any group-theoretic
     argument would have to prove is therefore **ZE$^{\rm all}_H$: $\mathrm{Bs}(V_H)$ contains no nondegenerate point**,
     which is not group-structured. The Hilbert-level claim that no single-form estimate can help is proved only for
     powers-only data and well-conditioned profiles ($\lambda\ge\lambda^*(n,\mu)$, checked at four values of $n$).
     Profile automorphisms are certified too, and once they are counted the argument fails near $m_{\max}$ (the
     claim is then unproved, not refuted), so it should not be cited as general.
     The actual point lies on no proper algebraic subgroup, but that does not affect a count at the trivial subgroup.
   * *A linear-algebra sufficient condition for route 5's target* (round 27 `interp.md`, `kill_interp.md`; **proved**,
     new to the record; the hypothesis is **open**). Let $L_{\rm pair}=\{\pm e_x,\pm(e_x-e_y)\}$. If, for all $e<m$, the
     degree-$e$ least space of the slab box $A$ is the sum of those of its shift-intersections $A\cap(A-l)$,
     $l\in L_{\rm pair}$, then $\mathrm{Bs}\,V(A,m)$ lies in the braid-and-coordinate arrangement (round 11's ZE-braid);
     with a top-degree height condition this gives $M=O_C(\log\log R)$. The hypothesis is at least as hard as
     ZE$^{\rm all}_H$. Evidence on slab boxes is mixed: the criterion beats every product (zonotope) certificate on the
     box M8$'$, but stays below the zonotope order on 20 of 27 random boxes.
   * *No tight pencils* (round 27 `stab.md` Thm 3.3, Cor 3.4, `kill_stab.md`; **proved**, killer-checked). Mean-width
     extremal measures are unique up to scalar on their support ($\dim A_t(S)\le1$ for integer $t\ge MW_n(S)$, all $n$).
     Hence no pencil of tight sections exists on $\mathrm{Bl}_eX$ (both forms, every $k$), and "second-order integrality"
     at Liouville-tight points is void. This closes the two-section part of the equality-case loophole of Theorem B;
     the single-form part still rests on the structure half of `threshold/construct.md` Conjecture E.
   * *Integrality at the root arrangement adds nothing at $e$* (round 27 dp2 Thm C and transfer, split Lemma 3.5, the
     dp2 killer's line lemma; **proved**). Bounded pair cofactors make near-tight tuples integral with respect to the
     diagonals, but on $\mathrm{Bl}_e\overline L_k$ that pair has log-canonical degree $-1$ on the lines through $e$ for
     every $k\ge3$, so Lang–Vojta predicts nothing there, and poles along the diagonals amplify nothing beyond sections
     wherever the section cap holds. A genuine fixed-$S$ Subspace argument on $LM_5\setminus\bigcup\Delta_{ab}$ is not
     covered.
4. **Single linear forms in logarithms** (Baker–Matveev, single-form Lang–Waldschmidt). Fakes satisfy them
   (`round8/lw.md` Thms A, B; `round8/stability.md` Prop 6.1; **proved**).
5. **$S$-unit and Subspace arguments with thresholds depending on $S$ only through $|S|$.** These give at best
   $M\le\omega(N')$ (lit18 Props A, B; **unrefereed**). Fixed-$S$ theorems such as `round18/theory.md` Thm 2.5 have
   $S$-dependent constants, and $S$ is the unbounded set of primes of $N$.
6. **abc on single relations.** abc through relations among products of chords (pair, triangle, Ptolemy or
   $n$-term) is vacuous on squarefree profiles (`round8/weakest.md` Prop 4.2), which contain every recorded fake
   (Cor 4.3), so it cannot prove growth (**proved**). abc on unit pencils is closed as well (Theorem 3.3).
7. **gcd normalisation.** It is invisible to forced-divisibility, auxiliary-form, Liouville, pair, cut, character
   and single-form inequalities (Lemma K1, proved by the killer). Vojta's $D=E$ term and the heuristic count do
   see it.
8. **Genus-0 Arakelov theory and unlikely intersections in tori.** The first reduces to items 1 and 3 (overview
   §5.4 rows 2 and 8) and to the limitation "Lim" of known tools (Mordell–Weil packing gives $e^{O(\mathrm{rank})}$; bounded rank is equivalent
   to the uniform theorem). The second is vacuous, or capped at $M\le\omega(N')+O_C(1)$; the Rémond-cone case of
   that ceiling is only **heuristic** (`round20/arakelov.md`, `round20/wild.md`; killer-checked;
   `overview/overview.md` §5.4 rows 18–19).
9. **For refutations.**
   * No feasible computation exists (`round18/census.md` Thm E).
   * No density mechanism works: proved at $k=8$, heuristic at $k=5$ for circles.
   * Bounded-degree covering curve families whose members each carry infinitely many balanced small-$H$ points
     are excluded at $k=8$ (genus $\le1$ by `round18/census.md` Prop 7.5; genus $\ge2$ by Faltings), and, for
     circles, those whose members each carry infinitely many A-points are excluded at $k\ge5$, $\delta>0$
     (Thm 4.3(ii)). Families with infinitely many rational points per member are not excluded as such (the lines of
     Kapranov's $\mathbb P^5$ are one; correction M1).
   * What remains are dense supplies with finitely many points per curve, curves of unbounded degree, and sporadic
     points.

## 6. Pointers

Repository `craigxchen/Jarnik`, branch `claude/upbeat-einstein-p0peju`, directory `docs/claude_review/`. Each
write-up `X.md` has its referee report `referee_X.md` next to it, or, in round 20, a kill report `kill_X.md`. Its
exact checks are in `X_checks/`. The round READMEs give the referee-corrected statements.

| item | file |
|---|---|
| consolidated account (status labels, corrections) | `overview/overview.md` (§§4, 5.4–5.5, 7.1, 7.1A, 7.6, 8.2) |
| Theorem 1.1, the $(2/3+\varepsilon)$ bound | `two_thirds/proof.md` (Thm 7.5, Cor 8.1, §10) |
| Vojta on $\mathrm{Bl}_YX_7$ | `growth/conditional.md` (Thms 5.1–5.2), `growth/referee_conditional.md` |
| $(L_8)$; Vojta on $\overline M_{0,8}$ implies uniform | `growth/ptolemy.md` (§6 Thm D, §7), `growth/referee_ptolemy.md` |
| $(L_5)$, $(L_6)$ false; Liouville barrier on $\overline M_{0,k}$ | `threshold/L8.md` (Thms A, B, Prop C), `threshold/referee_L8.md`, `growth/referee_ptolemy.md` §3.3 |
| round 18: census, sandwich, lattice form, lifts, Prop 2.6 | `round18/README.md`, `round18/theory.md`, `round18/census.md`, `round18/lift.md` (with `referee_*.md`) |
| round 19: Thm 2.4, Thm 3.3, non-unit residue | `round19/README.md`, `round19/theory.md`, `round19/search.md` (with `referee_*.md`) |
| round 20: $(L_5^\circ)$, $V_E(5)$, Thm 4.1, Lemma K1 | `round20/README.md`, `round20/gauss.md`, `round20/kill_gauss.md` |
| round 21: low-height window (open) | `round21/README.md` |
| fake circles and capped classes | `round4/outside.md`, `round5/construct.md`, `round8/lw.md`, `round10/analytic.md` |
| mean-width Thm B, Cor C3 (auxiliary forms capped) | `threshold/construct.md` |
| round 26: zero estimates vs route 5's ZE; spectral LP cap; target ZE$^{\rm all}_H$ | `round26/README.md`, `round26/ze.md` (with `referee_ze.md`) |
| round 27: dP2 calibration, uniqueness of extremal measures, split cones, multi-shift criterion | `round27/README.md` (with `kill_*.md`) |
| papers check: seven full texts on approximation constants | `papers_check/README.md` |

## Appendix: notes on the check

All 23 items of the adversarial check were applied. Three were applied in a modified form after re-reading the
sources.
* *Item 6 (unit class).* The check's suggested "equivalently $\zeta_a/\zeta_b=P/\bar P$ with $P\in\mathbb Z[i]$" is
  not an equivalence as stated: $i=(1+i)/(1-i)$ and $-1=i/(-i)$ are themselves of the form $P/\bar P$. The note
  therefore defines the class by equal units for a fixed choice of the $\pi_p$ (as in `growth/ptolemy.md` §1) and
  states the $P/\bar P$ form, with $P$ conjugate-primitive, only as a consequence.
* *Item 8 (labels).* Both remedies were used: the legend gains the label "proved by the killer", and Lemma K1 and
  the $k=3$ corollary carry it.
* *Item 17 ($\hat\gamma$, $D_{PQ}$).* The sentence was dropped in favour of the pointer to `gauss.md` Thm 4.2,
  Remark (3), as the check allowed.
