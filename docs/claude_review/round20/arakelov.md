# Round 20 — Arakelov / arithmetic-intersection lens

Status: complete (author's write-up; unrefereed). Labels: proved / exact / evidence / conditional / interpretation.

## 0. Summary and verdict

**Lens:** arithmetic intersection theory. The cluster is read as $M$ integral sections of the regular model
$\mathcal X_N$ of the conic $uv=Nt^2$ over $\mathbb Z[i]$ (and, for multi-point questions, as the point
$\Phi\in\overline M_{0,M}$ or the tuple $P\in X_M$).

**Verdict: reduced to known caps; no growth improvement, no uniform bound.**

* **Dictionary (Theorem A, proved, checked exactly on 26 configurations).** The exponent profile records which
  components of $\mathcal X_N$ (chains $C_0,\dots,C_{e_p}$ resolving $A_{e_p-1}$ singularities) the sections cross; the finite
  intersection numbers of two sections are the reduced chords, $(P_i\cdot P_j)_{\rm fin}=\log\mathrm{Norm}\,c_{ij}$; the
  vertical (Néron–Hriljac) correction of $D_0=\sum_jP_j-\frac M2(O_u+O_v)$ is the cut slack $Q_p$. The cut identity of the
  $(2/3+\varepsilon)$ proof is the genus-0 Hodge-index equality for $D_0$, written place by place.
* **Adelic Fekete (Theorem B, proved).** The $(2/3+\varepsilon)$ finite inequality is the adelic $M$-point diameter
  inequality $\sum_v\epsilon_v\le0$ for $E=\prod_vE_v$ (arc at $\infty$, local conic points at finite places); the chain lemma
  is the $p$-adic diameter of the $A_{e-1}$ chain.
* **Pair-data reduction (Lemmas 2.1–2.3, Prop B.2, proved).** On a genus-0 arithmetic surface every global relation is a
  pair relation and every finite local configuration is ultrametric, so Hodge index, Faltings–Hriljac, adelic energy,
  capacity, equidistribution energies and slope inequalities on $\mathbb P^1$ all lie in an "adelic-capacity class" ADC, which
  every $A_X$ fake satisfies: **capped by checklist item 1** (witness: row-2 fakes, $W\le(27+o(1))M\log M$).
* **Zhang / Northcott / arithmetic bigness (Prop 3.3(b), proved).** On a fixed model, in the $(H,E)$-plane, reaching the arc
  would need $\tau^*(M)>2$; **capped by item 3**. Equidistribution and Bogomolov are inapplicable (bounded degree, large
  height). Néron–Tate packing on covers gives counts $e^{O(\rho)}$ (item-5 shape).
* **Multi-point term (§6).** (i) The finite multi-point contact vector $\Gamma^{>M}_T$ ($|T|\ge3$, primes $>M$ not dividing
  $N$) is the only intersection-theoretic multi-point quantity that fakes do not carry; it is **zero on every record
  cluster** (evidence; 1 triple in 15 windows, at a non-extremal one) and zero on fakes' pseudo-adelic completions, so it does
  not separate. (ii) Theorem C (proved): every height relation on $\overline M_{0,M}$ obtainable from local data and product
  formulas is a combination of pair product formulas. The uncapped content of Vojta/$(L_8)$ is **global coherence**
  (rationality of the adelic configuration), which intersection numbers do not see: the record's conditional lane.
* **Constant-level remark (numerics).** Using the archimedean capacity of the arc instead of the chord would lower $K_0$ by
  about $4\log4$. Constants do not matter here.

Files: write-up `round20/arakelov.md`; scripts and outputs in `round20/arakelov_checks/` (`dictionary.py`, `adelic.py`,
`multipoint.py`, `clusters.py`, `gauss.py`, outputs `*_out.txt`, derivation notes `notes.md`).

## 1. Setup: the cluster as sections of an arithmetic surface

**Conventions.** $K=\mathbb Q(i)$, $\mathcal O_K=\mathbb Z[i]$. Notation of `overview.md` §1.2: a *normalised*
cluster is $z_1,\dots,z_M\in\mathbb Z[i]$, $|z_j|^2=N$, $\gcd(z_1,\dots,z_M)=1$, so that $N=\prod_p p^{e_p}$ is odd with
only split primes, $z_j=\varepsilon_j\prod_p\pi_p^{a_{jp}}\bar\pi_p^{e_p-a_{jp}}$, $W=\log N$,
$g_{ij}=\gcd(z_i,z_j)$, $c_{ij}=(z_i-z_j)/g_{ij}$, $d_{ij}=\sum_p|a_{ip}-a_{jp}|\log p$, $\gamma_{ij}=\log\mathrm{Norm}\,c_{ij}\ge0$,
$S_{p\tau}=\#\{j:a_{jp}<\tau\}$, $Q_p=\sum_{\tau=1}^{e_p}(S_{p\tau}-M/2)^2$, $Q=\sum_pQ_p\log p$ (cut slack).
Nothing below uses the arc hypothesis except where it says so.

**The arithmetic surface.** With $u=x+iy$, $v=x-iy$ the circle is the conic
$$\mathcal C_N:\ uv=Nt^2\subset\mathbb P^2_{\mathbb Z[i]} .$$
Its generic fibre is $\mathbb P^1_K$ through the parameter $s=u/t$ (inverse $s\mapsto[s^2:N:s]$). The cluster gives
$\mathbb Z[i]$-sections $P_j=[z_j:\bar z_j:1]$ ($s(P_j)=z_j$), and there are two distinguished sections, the isotropic
points $O_u=[1:0:0]$ ($s=\infty$) and $O_v=[0:1:0]$ ($s=0$).
* At Gaussian primes $\varpi\nmid N$ (including $1+i$, because $N$ is odd) $\mathcal C_N$ is smooth over $\mathbb Z[i]_{(\varpi)}$
  and $s$ identifies it with $\mathbb P^1$.
* At $\pi\mid p\mid N$, $p^{e}\,\|\,N$, the chart $t=1$ is $uv=\pi^{e}w$ ($w$ a unit at $\pi$): an $A_{e-1}$ singularity
  for $e\ge2$. Its minimal resolution has the charts $V_k=\mathrm{Spec}\,\mathcal O[\alpha_k,\beta_k]/(\alpha_k\beta_k-\pi)$,
  $u=\pi^{k-1}\alpha_k$, $v=\pi^{e-k}\beta_k$ ($k=1,\dots,e$), each regular. The special fibre is a chain
  $C_0,\dots,C_e$ of $\mathbb P^1_{\mathbb F_p}$'s, all of multiplicity one, with $C_k\cdot C_{k+1}=1$,
  $C_k^2=-2$ ($0<k<e$), $C_0^2=C_e^2=-1$; $C_a$ is the closure of the torus $\{v_\pi(u)=a,\ v_\pi(v)=e-a\}$, with
  coordinate $u/\pi^a\bmod\pi$.

Call $\mathcal X_N$ the resulting regular model (resolve at every $\pi\mid N$ and every $\bar\pi\mid N$).

**Theorem A (dictionary; proved, and checked exactly).** For a normalised cluster:
1. *(components = exponent profile)* At $\pi=\pi_p$ the section $P_j$ meets the component $C_{a_{jp}}$, and at $\bar\pi_p$
   it meets $C_{e_p-a_{jp}}$. $O_u$ meets $C_0$ and $O_v$ meets $C_{e_p}$ at both places.
2. *(horizontal finite intersection = reduced chord)* The finite part of the arithmetic intersection of two cluster
   sections is
   $$(P_i\cdot P_j)_{\rm fin}:=\sum_{\varpi}i_\varpi(P_i,P_j)\log\mathrm{Norm}\,\varpi=\log\mathrm{Norm}\,c_{ij}=\gamma_{ij}.$$
   At $\pi\mid N$, $i_\pi(P_i,P_j)=0$ if $a_{ip}\ne a_{jp}$, and $=v_\pi((z_i-z_j)/\pi^{a})$ if $a_{ip}=a_{jp}=a$; this is the
   collision criterion of `two_thirds/proof.md` §3 ("different levels never collide"), now a statement about
   components.
3. *(vertical correction = cut slack)* For a divisor $D=\sum_jc_jP_j+c_uO_u+c_vO_v$ of degree $0$, let $\Phi_\pi$ be the
   vertical $\mathbb Q$-divisor supported on the chain at $\pi$ with $(D+\Phi_\pi)\cdot C_k=0$ for all $k$. Then
   $$-\Phi_\pi^2=\sum_{\tau=1}^{e}S_\tau(c)^2,\qquad S_\tau(c)=\sum_{a(\cdot)<\tau}c_\cdot ,$$
   (the charge to the left of the cut $\tau$). For $c=e_i-e_j$ this is $|a_{ip}-a_{jp}|$, so the vertical corrections
   of pairs are exactly the profile distances $d_{ij}$. For
   $$D_0:=\sum_jP_j-\tfrac M2(O_u+O_v)$$
   it is $Q_p$ at $\pi$ and at $\bar\pi$.
4. *(the identity)* Exactly,
   $$\prod_{i<j}\mathrm{Norm}(z_i-z_j)^4=N^{M(M-2)}\prod_pp^{4Q_p}\prod_{i<j}\mathrm{Norm}(c_{ij})^4,$$
   which is the cut identity of `two_thirds/proof.md` (2.4).

*Proof.* (1) $v_\pi(z_j)=a_{jp}$ and $v_\pi(\bar z_j)=e-a_{jp}$ by the normal form. In the chart $V_{a+1}$ ($a<e$) the
point has $\alpha=u/\pi^a$ a unit and $\beta=v/\pi^{e-a-1}\equiv0$, so it reduces into $\{\beta=0\}=C_a$ away from the
nodes; for $a=e$ use $V_e$, where $\alpha\equiv0$ and $\beta=v$ is a unit. At $\bar\pi$ exchange the roles of $u,v$.
$O_u$ has $u=1$, so it lies on the component where $u$ is a unit ($C_0$); symmetrically for $O_v$.
(2) Sections through different components meet the special fibre at different points of the open tori, so they do not
meet. If both are on $C_a$, the ideal of $P_j$ in $V_{a+1}$ is $(\alpha-\alpha_j,\beta-\pi/\alpha_j)$, and the length of
$\mathcal O/(\alpha_i-\alpha_j,\pi/\alpha_i-\pi/\alpha_j)$ is $v_\pi(\alpha_i-\alpha_j)=v_\pi((z_i-z_j)/\pi^a)$. Since
$v_\pi(g_{ij})=\min(a_{ip},a_{jp})$ and $v_\pi(z_i-z_j)=\min(a_{ip},a_{jp})$ when the levels differ, in all cases
$i_\pi=v_\pi(c_{ij})$. At $\varpi\nmid N$, $i_\varpi=v_\varpi(z_i-z_j)=v_\varpi(c_{ij})$ because $\varpi\nmid g_{ij}$. Summing,
$(P_i\cdot P_j)_{\rm fin}=\log\mathrm{Norm}\,c_{ij}$.
(3) Write $\Phi=\sum x_kC_k$ and $n_k$ for the charge on $C_k$. The conditions read
$x_{k-1}-2x_k+x_{k+1}=-n_k$ (with the end rows $x_1-x_0=-n_0$, $x_{e-1}-x_e=-n_e$), so $x_{\tau}-x_{\tau-1}=-S_\tau(c)$,
and $-\Phi^2=\sum_\tau(x_\tau-x_{\tau-1})^2=\sum_\tau S_\tau(c)^2$ (summation by parts on the path graph). For $D_0$ the
charge left of the cut $\tau$ is $S_{p\tau}-M/2$ ($O_u$ sits on $C_0$). At $\bar\pi$ the levels are reflected and
$O_u,O_v$ stay at the ends, which gives the same sum.
(4) Place by place: at $\pi\mid N$, $\sum_{i<j}v_\pi(z_i-z_j)=\sum_{i<j}[\min(a_i,a_j)+i_\pi(P_i,P_j)]$ and
$\sum_{i<j}[\min(a_i,a_j)+\min(e-a_i,e-a_j)]=\binom M2e-\sum_{i<j}|a_i-a_j|=e\big(\tfrac{M^2}4-\tfrac M2\big)+Q_p$
(the cut formula $\sum_{i<j}|a_i-a_j|=\sum_\tau S_\tau(M-S_\tau)$). Summing over all places with (2) gives
$\sum_{i<j}\log\mathrm{Norm}(z_i-z_j)=\frac{M(M-2)}4W+Q+\sum_{i<j}\gamma_{ij}$; multiply by 4 and exponentiate. $\square$

*Exact check* (`arakelov_checks/dictionary.py` → `dictionary_out.txt`; big integers and exact rationals): (1)–(4) on 26
normalised configurations: the referee's primitive 10-point cluster ($N=5^313^2\,17\cdot29\cdot113$, $C=8.2393$) and its two
9-point sub-clusters, the record 8-, 7-, 6-, 5-point clusters of `growth/families/check_families.py`, the best 5-, 6-,
7-windows on three highly composite norms, and nine random (non-arc) subsets. In (2) the intersection numbers are
computed from the model (component test, chart coordinate) and compared with $\mathrm{Norm}\,c_{ij}$ computed by the
Gaussian Euclidean algorithm; in (3) $\Phi$ is computed by solving the chain system, for 5 random zero-sum charges per
prime and for $2D_0$. **All pass.**

**Reading (interpretation of a proved identity).** $D_0$ has degree $0$ on $\mathbb P^1_K$, hence is principal
($D_0=\frac12\mathrm{div}\,f$, $f=\prod_j(s-z_j)^2/s^{M}$). The arithmetic Hodge index theorem (Faltings, Hriljac) for
$D_0+\Phi$ on $\mathcal X_N$ is an **equality** in genus $0$ (the Néron–Tate height on the trivial Jacobian is $0$).
Theorem A identifies its terms: horizontal finite part $\sum\gamma_{ij}$, vertical part $Q$ (the cut slack),
archimedean part $\sum\log|z_i-z_j|$ against the $W$-terms. So **the cut identity of the $(2/3+\varepsilon)$ proof
is the Hodge index identity on the regular model of the conic, and the exponent profile is the record of which
components of $\mathcal X_N$ the sections cross.** The global equality itself is just the product formula for the
discriminant; Theorem A is what makes each term an intersection number.

## 2. Pair-data reduction: what a genus-0 arithmetic surface can see

### 2.1 Three structural lemmas

**Lemma 2.1 (genus 0: every global relation is a pair relation; proved).** Let $f\in K(\mathbb P^1)^*$. Over
$\bar K$, $f=\kappa\prod_\alpha(s-\alpha)^{m_\alpha}$. Hence for distinct points $z_1,\dots,z_M\in K$ off the support of
$\mathrm{div}f$, the product formula for $f(z_j)$ is the sum, with multiplicities $m_\alpha$, of the product formulas
for $z_j-\alpha$ over $K(\alpha)$ (grouped by Galois orbits). So every multiplicative global relation among values of
rational functions at the cluster splits into **pair relations**: cluster–cluster ($\alpha=z_i$, the pair identity
$2\log|z_i-z_j|=W-d_{ij}+\gamma_{ij}$) or cluster–auxiliary ($\alpha\notin\{z_i\}$, a single point against a fixed
algebraic point). *Proof.* $\bar K(s)^*$ is generated by constants and linear polynomials. $\square$

**Lemma 2.2 (finite places are ultrametric; proved, checked).** At a finite place $\varpi$ the reduction tree of
$(\mathbb P^1;z_1,\dots,z_M)$ (the stable model, equivalently the convex hull of the points in the Berkovich line) is
determined by the pairwise valuations $v_\varpi(z_i-z_j)$, which satisfy the ultrametric inequality. In particular
every boundary contact $\mathrm{ct}_\varpi(T)$ of a subset $T$ (the length of the edge of the tree separating $T$ from
its complement, which is the local height of $\Phi=[z_1:\dots:z_M]\in\overline M_{0,M}$ along $D_T$) is a min/max
expression in pair valuations. *Proof.* Standard (Berkovich tree of a finite set). Exact check of ultrametricity at every
prime $q>M$, $q\nmid N$, for the 15 arc clusters: `multipoint_out.txt`, column `ultra`, all True. $\square$

**Lemma 2.3 (genus-0 Hodge index is an equality; proved).** On any regular arithmetic surface with generic fibre
$\mathbb P^1_K$, for every arithmetic divisor $D+\Phi$ of degree $0$ orthogonal to all vertical components, the
Hodge index theorem gives $-(D+\Phi)^2=2[K:\mathbb Q]\,\hat h_{NT}([D_K])=0$. By polarization, a quadratic form in the
charges $c$ on $\{c:\sum c=0\}$ is determined by its values at $c=e_i-e_j$, so this equality is equivalent to the
family of **pair** identities (one per pair of sections). *Proof.* The Jacobian of $\mathbb P^1$ is trivial;
$B(c,c)=-\frac12\sum_{i,j}c_ic_jB(e_i-e_j,e_i-e_j)$ on $\sum c=0$ for any symmetric $B$. $\square$

### 2.2 The adelic Fekete form of the $(2/3+\varepsilon)$ theorem

For a place $v$ of $K$ (normalised: $|x|_\infty=|x|^2$, $|x|_\varpi=\mathrm{Norm}\,\varpi^{-v_\varpi(x)}$) and a set
$E_v$ of local points, put
$$\epsilon_v(M;E_v)=\inf\Big\{-\sum_{i<j}\log|x_i-x_j|_v:\ x_1,\dots,x_M\in E_v\text{ distinct}\Big\}.$$
Take $E_\infty$ = the arc (chord $\le C\sqrt R$), $E_\varpi$ = the local integral points of $\mathcal C_N$ (all of them at
$\varpi\nmid N$; at $p\mid N$ the levels may be prescribed or free).

**Theorem B (proved; the computation is the $(2/3+\varepsilon)$ proof re-read).**
1. *(local minima are the pigeonhole and chain brackets)* For $q\nmid 2N$ split, resp. inert,
   $\epsilon_q=2\log q\sum_{a\ge1}E(M,(q-1)q^{a-1})$, resp. $2\log q\sum_aE(M,(q+1)q^{a-1})$ (sum over both places
   above a split $q$). For $p\mid N$ with free levels,
   $$\epsilon_p=\log p\Big[e_p\Big(\tfrac{M^2}4-\tfrac M2\Big)+\min_{(n_t)}\Big(Q_p(n)+2\sum_t\sum_aE(n_t,(p-1)p^{a-1})\Big)\Big]
   =\log p\Big[e_p\Big(\tfrac{M^2}4-\tfrac M2\Big)+2\Psi^{(e_p)}_p(M)\Big],$$
   where $\Psi^{(e)}_p(M)$ is the level-collision-plus-half-cut bracket of `two_thirds/proof.md` minimised over level
   vectors with $e$ levels ($\Psi_p=\min_e\Psi_p^{(e)}$ there).
   At $(1+i)$ the proof of the $2/3$ theorem uses only $\epsilon_{1+i}\ge\binom M2\log2$. At $\infty$,
   $\epsilon_\infty\ge-2\binom M2\log\delta_M(E_\infty)$ with $\delta_M$ the $M$-th diameter.
2. *(the inequality)* The product formula $\sum_v\sum_{i<j}\log|z_i-z_j|_v=0$ gives $\sum_v\epsilon_v\le0$. With the
   crude bound $\delta_M(E_\infty)\le C\sqrt R$ this is exactly the finite inequality
   $F^*_N(M)\le\frac M4\log R+\binom M2\log C$ of `two_thirds/proof.md` Thms 5.1–5.2 (for the given $N$; minimising over
   the patterns of $N$ gives $F^*(M)$). So **the $(2/3+\varepsilon)$ theorem is the adelic $M$-point (Fekete) diameter
   inequality for $E=\prod_vE_v$**, and the chain lemma computes the $p$-adic $M$-point diameter of the
   $A_{e-1}$ chain of Theorem A.

*Proof.* (1) At $q\nmid N$ the local conic is smooth, so by Hensel each class mod $q^a$ of local points has the same
number of lifts mod $q^{a+1}$, the class counts are $(q\mp1)q^{a-1}$ (`two_thirds/proof.md` §3), and
$-\log|z_i-z_j|_q=2\log q\cdot\#\{a:z_i\equiv z_j\bmod q^a\}$ (collisions at $\varpi$ and $\bar\varpi$ coincide for points of
equal norm prime to $q$). The minimum of $\sum_a\#\{\text{colliding pairs mod }q^a\}$ is attained by a nested
balanced allocation, which is balanced at every depth simultaneously, and equals $\sum_aE(M,Q_a)$. At $p\mid N$, by the
computation in Theorem A(4), the local energy is
$\log p[e(\frac{M^2}4-\frac M2)+Q_p+2\,\mathrm{coll}]$ where $\mathrm{coll}$ counts same-level collisions; minimising over
level vectors and, per level, over balanced allocations gives the bracket. (2) Insert (1) into $\sum_v\epsilon_v\le0$,
divide by 2 and compare with `two_thirds/proof.md` (5.1), using $\sum_pe_p\log p=W$. $\square$

*Exact check* (`arakelov_checks/adelic.py` → `adelic_out.txt`): on all 26 configurations the product formula
$\sum_v\mathrm{act}(v)=0$ holds to $10^{-13}$, and $\mathrm{act}(q)\ge\epsilon_q$ at every prime (assertion in the script).

**Remark B.1 (the archimedean capacity; a constant only).** Replacing the chord bound by the $M$-th diameter of the arc
(an arc of angular width $\theta$ has capacity $R\sin(\theta/4)\le L/4$) lowers the archimedean side by
$\binom M2\log(2/\delta_M([-1,1]))\to\binom M2\log4$. Numerically (Lobatto nodes, `adelic_out.txt`)
$\delta_M([-1,1])=0.7017,\,0.5495,\,0.5153$ at $M=10,50,200$ (limit $1/2$). In Theorem 2.1 of the overview this would lower
$K_0$ by about $4\log4\approx5.5$ plus an $O(\log M)$ correction. **Constant only; not pursued** (constants do not matter
here, and the $O(\log M)$ correction is stated from numerics, not proved).

### 2.3 The adelic-capacity class is capped

**Definition.** The *adelic-capacity class* ADC consists of all inequalities obtained from
(i) the per-pair product formula, used in the form "finite contact of a pair $\le$ its archimedean budget",
(ii) at each place separately, any lower bound on a nonnegative combination $\sum_{i<j}\lambda_{ij}(\text{local contact})$
valid for all $M$-tuples of local points with the cluster's local data (levels at $p\mid N$; arc at $\infty$), and
(iii) the exponent profile. Every Arakelov-intersection inequality on a genus-0 arithmetic surface model of $\mathbb P^1_K$
(Hodge index, Faltings–Hriljac, adelic energy, Fekete–Szegő–Cantor–Rumely capacity, energy forms of quantitative
equidistribution, Bost's slope inequality for evaluation at the cluster points on $\mathbb P^1$) lies in ADC, by
Lemmas 2.1–2.3: its global content is per-pair product formulas (Lemma 2.1, 2.3), its local content at finite places is
a function of pair valuations (Lemma 2.2), and its local content at $\infty$ is a function of the real positions.
(Local intersection numbers at genuine places are integers, and ADC may use that. What ADC does not include is the
exact equality "total pair number $=\log$ of an integer" for all pairs at once; see Theorem C(3), *Scope*.)

**Proposition B.2 (ADC is capped; proved, given the refereed row 2 of the overview's §5.4).** Let $X\ge M$. Every member
of the class $A_X$ (angle model with point-level residues at moduli $\le X$, each character with its actual collision
modulus) satisfies every ADC inequality whose local inputs are taken at places $q\le X$, $q\mid N$ or $q=\infty$.
Hence ADC cannot prove $M=o(\log R/\log\log R)$; the fully flipped Paley fakes of `round5/construct.md` give the witness
$W\le(27+o(1))M\log M$.

*Proof.* (a) At $\infty$ an angle model has genuine real points on an arc, so it satisfies every archimedean input.
(b) At $q\le X$ its residue design is a local configuration with the actual class counts, so every local lower bound
of type (ii) holds for the designed contacts: such a bound depends only on the tree of residue classes of local points,
which by Hensel is the regular tree with the actual class counts (proof of Theorem B(1)). (c) At $q>X\ge M$,
$q\nmid N$: the local minimum of every nonnegative combination is $0$, since there are at least $q-1\ge M$ classes
(split) or $q+1>M$ (inert), so ADC gains nothing there. (d) At $p\mid N$ the profile terms are genuine; the record's
witnesses have all $p_j>X\ge M$, so level collisions have local minimum $0$. (e) The pair budget: the fake's
"reduced chord" $\gamma^{\rm fake}_{ij}:=2\log|z_i-z_j|-(W-d_{ij})$ is $\ge$ its designed small-prime contacts; this is
the pair case $c=e_x-e_y$, $u=1$ of the residue axiom (R) of `round5/construct.md` §1 ("the collision-strengthened pair
gap of `two_thirds`, stated for angles"), which the witnesses satisfy. Their residues $\rho_x=s\,\omega_x/\bar\omega_x$
lie in one norm class mod $q^a$, i.e. they are local points of the conic, with the actual class counts. The remainder plays the role of
the contacts at primes $>X$, where (c) shows ADC asks for nothing. $\square$

So every genus-0 Arakelov argument is capped by item (1) of the checklist, and its best possible output is the
$(2/3+\varepsilon)$ inequality up to the archimedean constant of Remark B.1. (Whether pseudo-configurations attain all
local minima simultaneously, which would make $2/3$ the exact optimum of ADC, is the question the overview's §2.1 records
as "not examined"; the cap does not depend on it.)

## 3. Candidate inequalities, one by one

Each candidate is classified as **pair** (reduces to pair data and is in ADC, hence capped by Prop B.2),
**$\tau^*$** (reduces to geometric bigness on a fixed model, capped by checklist item 3), **Lim** (output exponential
in a rank or in $|S|$, checklist item 5 type), **inapplicable**, or **multi-point** (the term is stated precisely and
tested in §6).

### 3.1 Arithmetic Hodge index / Faltings–Hriljac

* *On $\mathcal X_N$ (genus 0).* An equality (Lemma 2.3), identified with the cut identity by Theorem A. **Pair.**
* *On covers.* Let $\phi:Y\to\mathbb P^1_K$ be a finite cover of genus $g\ge1$ and $Q_j\in Y(L)$ lifts of the $P_j$
  (over some finite $L/K$; in general the field of definition grows with $j$, e.g. $L=K(\sqrt{f(z_j)})$ for a double
  cover).

  **Proposition 3.1 (proved, modulo standard Faltings–Hriljac).** (a) $\hat h_{NT}(Q_i-Q_j)=-(Q_i-Q_j+\Phi_{ij})^2/[L:\mathbb Q]$
  on a regular model of $Y_L$. Every term is a function of the pair: the archimedean Green function of $Y$ at
  $(Q_i,Q_j)$ is a function of the archimedean positions of $z_i,z_j$ and of the sheet choices (fake-available); at a
  finite place the intersection and the vertical correction depend only on the reductions of $Q_i,Q_j$, i.e. on
  $v(z_i-z_j)$, on $v(z_\cdot-b)$ for the branch points $b$ (single point against a fixed point) and on the sheet of
  $Q_\cdot$ modulo the place (a residue/character datum of $f(z_\cdot)$). So the pairwise Néron–Tate data are pair plus
  residue/character data: **checklist item 1**.
  (b) The one genuinely multi-point content is positivity: the matrix $(\hat h(Q_i-Q_j))$ is of negative type
  (Schoenberg), with Euclidean rank $\rho\le\mathrm{rk}\,\langle Q_j-Q_1\rangle$. Its only known use for counting is the
  Mumford–Bombieri–Vojta packing argument, whose output is $M\le e^{O(\rho)}$ for points of nearly equal height that
  are pairwise close at one place. Here $\rho$ is not bounded by anything smaller than $M-1$; in the torus analogue
  ($u_j=z_j/z_1$ in the norm-one torus, height $h(u_j)=d_{1j}/2$, an $\ell^1$ norm on exponent vectors) the rank is
  the rank of the exponent lattice, and bounding it is **equivalent** to the uniform theorem (overview Thm 3.12). **Lim**
  (same shape as checklist item 5: a count exponential in a rank that the fakes, and possibly actual clusters, make
  $\asymp M$).

### 3.2 Arakelov adelic energy (capacity, transfinite diameter, Fekete–Szegő–Rumely)

Theorem B: the adelic $M$-point diameter inequality **is** the $(2/3+\varepsilon)$ finite inequality; the chain lemma is
the $p$-adic diameter of the $A_{e-1}$ chain. The energy is a quadratic form in point masses, so by polarization it is
pair data; its local minima depend only on class counts. **Pair; capped (Prop B.2)**, and its best output is the
$2/3$ theorem up to the archimedean constant of Remark B.1. The adelic Fekete–Szegő *existence* direction (capacity
$>1$ gives infinitely many points) concerns Galois orbits of growing degree and says nothing about a single $N$.

### 3.3 Zhang successive minima, arithmetic Bogomolov, equidistribution of small points

**Proposition 3.3 (proved).**
(a) *Equidistribution is inapplicable or pair.* The theorems of Bilu, Szpiro–Ullmo–Zhang, Chambert-Loir, Thuillier,
Favre–Rivera-Letelier, Baker–Rumely and Yuan concern Galois-stable $0$-cycles of degree $\to\infty$ whose normalised
height tends to the essential minimum. The cluster points are $K$-rational (Galois orbits of size $\le2$) and their
heights are $\ge W/4-\log C$ (pair inequality), not small. The version for arbitrary finite sets (Fekete-type
equidistribution, quantitative discrepancy bounds) is controlled by the pair energy $\frac1{M^2}\sum_{i\ne j}g_\mu(x_i,x_j)$,
which is in ADC. **Inapplicable / pair.**
(b) *Zhang/Northcott/arithmetic bigness on a fixed model cannot reach the arc, in the $(H,E)$-plane.* Let $X'$ be a fixed
smooth projective model of $\mathrm{Bl}_YX_M$ and $\bar L$ an adelic $\mathbb Q$-line bundle with class $L=aH-bE$,
$a>0$, $b\ge0$, and fixed model metrics at almost all places. For short-arc tuples,
$h_{\bar L}(P^\#)\le(a-\tfrac b2)\log R'+b\log C+O(1)$ (because $h_E\ge\lambda_{E,\infty}-O(1)\ge\frac12\log R'-\log C-O(1)$,
`growth/conditional.md` Lemma 4.1). A non-density conclusion from smallness of $h_{\bar L}$ (Northcott for big $L$,
Zhang's lower bound on the essential minimum, or arithmetic bigness $\widehat{\mathrm{vol}}(\bar L)>0$, which implies
that $L$ is big) works as follows: for some $\varepsilon>0$ the class $L-\varepsilon H=(a-\varepsilon)H-bE$ is big, so off a fixed proper closed set
$h_{\bar L}\ge\varepsilon h_H-O(1)=\varepsilon\log R'-O(1)$, and this must contradict the upper bound for large $R'$. That needs
$a-\frac b2<\varepsilon$, while bigness of $(a-\varepsilon)H-bE$ needs $b/(a-\varepsilon)<\tau^*(M)$, where
$\tau^*(M)=2-1/\lceil M/2\rceil$ is the pseudo-effective threshold of $H-tE$ (`threshold/upper.md`). Together:
$b/\tau^*<a-\varepsilon<b/2$, i.e. $\tau^*>2$ (if $b>0$; if $b=0$ the first condition fails). Impossible. **$\tau^*$.**
(Zhang's successive-minima bound and arithmetic bigness, $\widehat{\mathrm{vol}}(\bar L)>0\Rightarrow L$ big, enter
only through this geometric bigness, the fixed metrics contributing $O(1)$.)
Classes that also involve exceptional divisors over the diagonals $\Delta_{ij}$ or the loci $\{u_j=0\}$ pick up finite
contacts (collisions, levels), i.e. auxiliary forms with pigeonhole-forced finite vanishing: this is the determinant
method of checklist item 6 / overview row 16, and is not re-capped here. Finite metrics that depend on $N$ or on the
cluster would make the essential minimum, and hence the threshold, depend on the template (**excluded by the project
constraint**).
(c) *Bogomolov on tori is inapplicable.* For $u_j=z_j/z_1$ in the norm-one torus, $h(u_j)=d_{1j}/2\ge W/4-\log C$, while
Bogomolov/Zhang concern points of height $\to0$. The relevant quantity is the ratio of archimedean proximity to
height: per pair it is Liouville's bound on the torus, attained infinitely often; jointly it is the exact Dirichlet
exponent $1/2$ for the diagonal line of $X_M$ (`round4/hardness.md`; §5 below), where Roth/Schmidt-type and
successive-minima methods give nothing.

### 3.4 Arithmetic Riemann–Roch and Bost's slope inequalities

**Lemma 3.4 (proved, elementary).** On $\mathbb P^1_{\mathbb Z[i]}$ (or $\mathcal X_N$ away from $N$) the sub-lattice of
$H^0(\mathcal O(D))$ of sections vanishing to order $\ge m_j$ at $z_j$ is $V\cdot H^0(\mathcal O(D-\sum m_j))$ with
$V=\prod_j(s-z_j)^{m_j}$, which is monic, hence primitive. So the slope filtration by vanishing at the cluster is
explicit, and every slope inequality for it reduces to the size of $V$ at $\infty$ (a product of the
$|s-z_j|$) and the product formula: **pair**. With derivatives (confluent evaluation) the evaluation determinant is a
confluent Vandermonde times a symmetric Gaussian integer (`growth/auxiliary.md` Prop A, unrefereed; overview row 16).

*Proof.* A polynomial vanishing at $z_j$ to order $m_j$ is divisible by $(s-z_j)^{m_j}$ in $K[s]$, and by Gauss's
lemma in $\mathbb Z[i][s]$ since $V$ is monic. $\square$

On $X_M$ (the tuple as one point) the slope method is Liouville at one place (**$\tau^*$**, checklist item 3); with
finite-place vanishing filtrations it is the determinant method (**checklist item 6**).

### 3.5 The blow-ups of `growth/conditional.md`

In arithmetic-intersection language the conditional theorem is the statement $h_{K_{X^\#}}(P^\#)\le\varepsilon h_H(P^\#)+O(1)$
off a fixed $Z$, i.e. Vojta's $D=0$ inequality on one fixed variety; Prop 3.3(b) shows that no positivity statement
in the $(H,E)$-plane proves it, and arithmetic intersection theory supplies no proved $D=0$ inequality on
$\mathrm{Bl}_YX_7$ or $\overline M_{0,8}$ (the arithmetic Hodge index theorem on higher-dimensional arithmetic varieties,
Moriwaki/Yuan–Zhang, concerns divisors numerically trivial on the generic fibre and gives no height inequality for
points). **Conditional** (unchanged).

## 4. Explicit computations (extremal clusters and fake circles)

All scripts are in `arakelov_checks/`; integer arithmetic is exact, logarithms are floats (labelled).

**4.1 Dictionary (exact; `dictionary.py` → `dictionary_out.txt`).** Theorem A (1)–(4) on 26 configurations; all pass
(§1).

**4.2 Place-by-place slack of the adelic inequality (`adelic.py` → `adelic_out.txt`; floats; evidence).** For each
cluster, $\mathrm{slack}=-\sum_v\epsilon_v$ (local minima, archimedean term = Fekete diameter of the cluster's own chord) is
split into the excess of the actual local energy over the local minimum at $\infty$, at small primes ($q\le M$ or
$q=2$, $q\nmid N$), at $p\mid N$, and at large primes ($q>M$, $q\nmid N$). The product formula holds to $10^{-13}$ on every
row. Arc clusters from the record:

| cluster | $M$ | $W$ | $C$ | slack | $\infty$ | $q\le M$ | $p\mid N$ | $q>M$ | Fekete gain over chord |
|---|---|---|---|---|---|---|---|---|---|
| ten (1176852625) | 10 | 20.89 | 8.239 | 87.97 | 7.70 | 27.73 | 14.79 | 37.76 | 94.26 |
| eight axis 5743586 | 8 | 31.13 | 5.674 | 76.30 | 1.38 | 29.43 | 14.02 | 31.47 | 54.97 |
| eight axis 590749322 | 8 | 40.39 | 5.754 | 95.62 | 0.25 | 13.86 | 0.00 | 81.51 | 54.97 |
| seven sporadic | 7 | 19.28 | 5.270 | 38.26 | 7.72 | 16.75 | 3.22 | 10.57 | 39.31 |
| six axis | 6 | 26.88 | 4.000 | 40.83 | 6.45 | 13.98 | 0.00 | 20.40 | 26.28 |
| five sporadic | 5 | 27.25 | 1.996 | 16.08 | 4.75 | 11.33 | 0.00 | 0.00 | 15.87 |

*Readings (evidence).* (i) The record clusters are close to Fekete-distributed on their own chord: the archimedean
excess is $0.25$–$7.7$ against a Fekete gain of $16$–$94$. (ii) The slack sits at the finite places, and the share at
large primes $q>M$ ranges from $0$ to $85\%$ ($81.5/95.6$ for the second non-golden 8-point cluster). Large-prime
contacts are exactly the data on which ADC (and so every genus-0 Arakelov inequality) asks nothing (Prop B.2(c)).
(iii) The non-arc random subsets in the same file have negative "archimedean excess" because their chord is not an arc
of length $O(\sqrt R)$; they only exercise the finite-place code.

**4.3 Multi-point contacts at finite places (`multipoint.py` → `multipoint_out.txt`; exact counts, float expectations;
evidence).** For each arc cluster and each prime $q>M$, $q\nmid 2N$: ultrametricity of $v_q(z_i-z_j)$ (all True), the
number of pairs and of triples lying in one class mod $q$, and the expectation in the model "independent uniform classes
among $m_q=q\mp1$, restricted to $q\le|c_{ij}|$" (a prime $q\nmid N$ can divide $z_i-z_j$ only if it divides the reduced
chord, so only if $|c_{ij}|\ge q$). Totals over the 15 arc clusters: **pairs 28 observed vs 25.5 expected; triples 1
observed vs 0.02 expected** (the one triple is in a non-extremal 6-window, `w6_2576450045`). On the 10-point cluster:
4 vs 2.8 pairs, 0 vs 0.000 triples.

*Readings (evidence).* Given the sizes of the reduced chords, large-prime pair contacts are as random, and
large-prime multi-point contacts ($|T|\ge3$) are absent on every record cluster. The reason is size: by (UB),
$\sum_{i<j}\log|c_{ij}|\le\frac M8W+\binom M2\log C$, so at $W=cM\log M$ the average reduced chord is
$|c_{ij}|\lesssim C\,M^{c/4}$, and a prime $q$ can be a common contact of a triple only if three reduced chords exceed
$q$. **So on actual extremal clusters the finite multi-point term outside $N$ is (empirically) zero, and an adelic
pseudo-configuration built from a fake by giving each pair a private large place has it zero as well. It does not
separate actual clusters from fakes.**

**4.4 What was not computed, and why.** (a) No fake circle was constructed here: Prop B.2 is a proof that every
$A_X$ member satisfies every ADC inequality, so a numerical fake test would test nothing new. (b) No Vojta/$(L_8)$
evaluation at $k=8$: census Thm E (round 18) proves balanced data need minors $\ge10^{49}$; the record clusters have
$w=W/2^7<1$, where every height is $O(1)$-dominated.

## 5. Barrier checklist, item by item

The lens produced no proposal that escapes the checklist; this section records, for each item, where the lens's
candidates land.

1. **Inputs fake circles satisfy (pairs, gcds, cut, residues, characters, GAP-Fourier, energy, $L^p$).** Every
   genus-0 Arakelov inequality is in ADC (§2.3) and ADC is capped by the row-2 fakes (Prop B.2). The Hodge-index identity
   on $\mathcal X_N$ *is* the cut identity, and the adelic energy *is* the pigeonhole/chain bracket (Theorems A, B).
   Néron–Tate data on covers are pair plus residue/character data (Prop 3.1(a)). **Capped by item 1.**
2. **Robust/averaged analytic statements.** Quantitative equidistribution is an energy (pair) statement and
   Galois-orbit equidistribution is inapplicable (Prop 3.3(a)). Nothing here is an ensemble statement. **Not used.**
3. **Single-place auxiliary polynomials on $X_M$; auxiliary polynomials on $\overline M_{0,k}$.** Zhang/Northcott/
   arithmetic bigness on a fixed model in the $(H,E)$-plane needs $\tau^*(M)>2$ (Prop 3.3(b)); slope inequalities on $X_M$
   are Liouville. **Capped by item 3.**
4. **Single linear form in logarithms; SLW.** Not used. (The torus heights of §3.1(b) are $\ell^1$ profile norms, not
   linear forms.)
5. **$S$-integer arguments with thresholds depending only on $|S|$.** Mordell–Weil packing gives counts $e^{O(\rho)}$ in
   a rank $\rho$ that is not bounded below $M-1$ (Prop 3.1(b)); the Vojta/Subspace form of global coherence for fixed $S$
   has $S$-dependent exceptional sets (§6.2). **Same shape as item 5; not used.**
6. **Route 5's box statement, CIRCLE$(q)$.** Classes in Prop 3.3(b) that involve exceptional divisors over the
   diagonals, and slope filtrations with finite-place vanishing, are the determinant method. **Not reopened.**
7. **$(L_5)/(L_6)$ family mechanisms at $k=8$; computational tests of $(L_8)$.** Not attempted; §4.4(b) explains why no
   $k=8$ height was evaluated.
8. **The abc unit-pencil lane.** Not used. (Lemma 2.1 shows that on $\mathbb P^1$ every multiplicative relation is a
   product of pair relations; abc on the three-term relation $(z_1-z_2)+(z_2-z_3)=(z_1-z_3)$ is the vacuous squarefree case
   of `round8/weakest.md` Props 4.1–4.2.)

**Dirichlet exponent (`round4/hardness.md`).** In the lens, per pair: for $u_{ij}=z_i/z_j$ in the norm-one torus,
$\lambda_\infty(u_{ij},1)=\log R-\log|z_i-z_j|=\frac12(d_{ij}-\gamma_{ij})$ and $h(u_{ij})=d_{ij}/2$, so
$\lambda_\infty\le h$ is Liouville's inequality on the torus, attained up to $O(1)$ by infinitely many $u$ (those
$P/\bar P$ with $|\mathrm{Im}\,P|=1$); a pair statement has no Diophantine content beyond it, and the arc condition
$\lambda_\infty\ge W/4-\log C$ is just the pair inequality. Jointly, on $X_M$ the arc condition is proximity of exponent
exactly $1/2$ to the diagonal line (overview §4.1, non-density equivalence), the Dirichlet exponent, where Roth/Schmidt
methods (which need $1/2+\varepsilon$) and successive-minima methods give nothing. Nothing in the lens moves this.

**Project constraints.** No exceptional set of a fixed-$S$ theorem is used; Gaussian gcd and denominator costs are kept
exactly (Theorem A(2): $\gamma_{ij}=\log\mathrm{Norm}\,c_{ij}$ with the gcd removed, the $(1+i)$ contact and the vertical
corrections all kept); no phase is set to zero (the archimedean data are the real positions); no template-dependent
threshold is treated as uniform (Prop 3.3(b) explicitly excludes $N$-dependent metrics).

## 6. What survives: the precise multi-point term

The task asks: if a candidate does not reduce to pair data, state the multi-point term and test it. Two candidates
for a multi-point term arise in this lens. The first is precise and tested; the second is the record's known
conditional lane.

### 6.1 The finite multi-point contact vector (precise; tested; does not separate)

For $T\subseteq[M]$, $2\le|T|\le M-2$, and a finite place $\varpi$, let $\mathrm{ct}_\varpi(T)\ge0$ be the contact of $T$ in the
$\varpi$-adic stable tree of $(\mathbb P^1;z_1,\dots,z_M)$. With the smooth integral model $\overline M_{0,M}/\mathbb Z$
(Keel), whose boundary is a normal-crossings divisor, the finite local height of the point
$\Phi=[z_1:\dots:z_M]$ along $D_T$ is $\mathrm{ct}_\varpi(T)\log\mathrm{Norm}\,\varpi$, so
$$h_{D_T}(\Phi)=\sum_\varpi\mathrm{ct}_\varpi(T)\log\mathrm{Norm}\,\varpi+\lambda_{D_T,\infty}(\Phi)+O(1).$$
The **multi-point term** is the vector $\Gamma_T=\sum_{\varpi}\mathrm{ct}_\varpi(T)\log\mathrm{Norm}\,\varpi$ for $|T|\ge3$.

**Proposition 6.1 (proved).** (a) $\Gamma$ is pair-determined place by place (Lemma 2.2), but it is **not** a linear
function of the global pair numbers $\gamma_{ij}$: at each place $v_\varpi(z_i-z_j)$ is the depth of the meet of $i,j$ in the
affine (rooted) tree, i.e. a sum of edge lengths over the clusters containing $\{i,j\}$, and the map from pair
valuations back to contacts is piecewise linear (min/max), place by place. ADC sees only the sums over places of the
pair valuations (Prop B.2), i.e. a linear projection of the tree data that forgets the place-wise min/max structure.
(b) Split $\Gamma=\Gamma^{N}+\Gamma^{\le M}+\Gamma^{>M}$ by places dividing $N$, places over $q\le M$, and places over
$q>M$, $q\nmid N$. Then $\Gamma^N$ is profile data (levels; the Ptolemy blocks $n_T$ of `growth/ptolemy.md` are built
from it), and $\Gamma^{\le M}$ is point-residue data at moduli $\le M$; both are carried by $A_X$ fakes ($X\ge M$).
Only $\Gamma^{>M}$ is not fake data.

**Test (§4.3; evidence).** $\Gamma^{>M}_T=0$ for every $|T|\ge3$ on all 6 record arc clusters (10, 8, 8, 7, 6, 5 points)
and on 8 of the 9 best windows (one triple, at the non-extremal window `w6_2576450045`), in line with the random model given the reduced-chord sizes (expected
total $0.02$). For fakes, the canonical pseudo-adelic completion (each pair's unused budget placed at a private large
place) also has $\Gamma^{>M}_T=0$ for $|T|\ge3$. **So the only multi-point term that the intersection numbers carry and
fakes do not, is zero on both sides in every available example. It cannot be the separating input.** (A proof that
it vanishes is not claimed; only the evidence and the size argument of §4.3.)

### 6.2 Global coherence (the Vojta / $(L_8)$ content)

**Theorem C (proved, from standard facts on $\overline M_{0,M}$).** Let $M\ge4$.
1. $\mathrm{Pic}(\overline M_{0,M})$ is generated by the boundary divisors $D_T$ (Keel). Hence every Weil height on
   $\overline M_{0,M}$ is, up to $O(1)$, a rational combination of the $h_{D_T}(\Phi)$, i.e. of local contacts
   (finite places: §6.1; $\infty$: a function of the real positions).
2. Every linear equivalence among boundary divisors is the divisor of a unit of $M_{0,M}$, and the units are, up to
   constants, products of cross-ratios, i.e. of differences $z_i-z_j$ (in Kapranov's chart $z_M=\infty$, $M_{0,M}$ is the
   complement of the braid arrangement). So the product formula for any such unit at $\Phi$ is a $\mathbb Z$-combination of
   the pair product formulas $\sum_v\log|z_i-z_j|_v=0$.
3. Consequently every inequality among heights of $\Phi$ that follows from local inequalities at each place and the
   product formulas of rational functions on $\overline M_{0,M}$ holds for every **pseudo-adelic configuration**:
   independent local configurations $E_v^M$ at every place of $K$, plus, for each pair, at most one *formal place* of
   positive real weight at which only that pair collides, tied only by the per-pair product formulas. Every $A_X$ fake
   ($X\ge M$) has such a completion (Prop B.2(e): no collisions at the genuine places $q>X$, which is possible since
   there are $\ge M$ classes, and each pair's unused budget, which is $\ge0$ by axiom (R), at its formal place).
   *Scope.* Formal places carry no integrality: arguments that use, with equality, that the global pair number
   $\mathrm{Norm}(z_i-z_j)$ is an **integer** are not covered. That integrality is not a pair *inequality*; for all pairs at
   once it is exact global (additive) data of the kind the overview lists as uncapped in §5.5 ("additive data [...] not
   modelled by angle models. Used in full they are the problem itself"). If one
   insists on genuine places only, replace each formal place by a private large prime of weight $\le$ the budget, and
   the per-pair product formula becomes the inequality "finite $\le$ archimedean", which is the direction every
   capacity or height-lower-bound argument uses.

*Proof.* (1) Keel, *Intersection theory of moduli space of stable n-pointed curves of genus zero* (Trans. AMS 1992):
$\mathrm{Pic}$ is free on the boundary modulo the Keel relations; heights are additive up to $O(1)$. (2) The Keel relations
are pull-backs of $D_0+\dots=D_1+\dots=D_\infty+\dots$ on $\overline M_{0,4}=\mathbb P^1$ under forgetful maps, i.e. the
divisors of cross-ratios; a unit on the arrangement complement is a monomial in the defining linear forms. (3) Sum (2)
over the relations used. $\square$

**What this means.** Arakelov intersection theory on $\overline M_{0,M}$ (or on any blow-up of $X_M$ along diagonals,
by the same argument with the linear forms $u_i-u_j$, $v_i-v_j$) produces no inequality that pseudo-adelic
configurations, and hence fakes, violate. The content of Vojta's $D=0$ inequality on $\overline M_{0,8}$, and of $(L_8)$,
beyond this is **global coherence**: the statement concerns *rational* points off a fixed proper closed $Z$, i.e. adelic
configurations whose local pieces are reductions of one global configuration (in particular with integral global pair
numbers, Theorem C(3), *Scope*). This is the multi-point term in its
precise form, and it is exactly the record's conditional lane: in its full strength (lattice gaps at every prime) it is
the problem itself (`round10/analytic.md` Thm D); in a fixed-$S$ form it is Subspace with $S$-dependent exceptional
sets (checklist item 5, project constraint); in the Vojta form it is open, untestable at $k=8$ (census Thm E) and not
provable by auxiliary polynomials on $\overline M_{0,k}$ (checklist item 3). The additive Ptolemy/Plücker relations,
which $(L_8)$ uses through its exact residues, are additive relations among values and are invisible to intersection
numbers altogether.

**Test on clusters and fakes.** For actual clusters global coherence holds by definition; for fakes it fails by
definition (angle models have no global points). A numerical test is therefore empty; the only testable proxy is the
finite multi-point vector of §6.1, which does not separate.

## 7. Conclusion

**No growth improvement and no uniform bound.** The lens is **reduced to known caps**:

| candidate | reduces to | status |
|---|---|---|
| arithmetic Hodge index / Faltings–Hriljac on $\mathcal X_N$ | the cut identity (Theorem A) | pair; capped (Prop B.2) |
| Faltings–Hriljac on covers of genus $\ge1$ | pair + residue/character data; packing $e^{O(\rho)}$ | item 1; Lim (Prop 3.1) |
| adelic energy, capacity, Fekete–Szegő–Rumely | the $(2/3+\varepsilon)$ finite inequality (Theorem B) | pair; capped (Prop B.2); best output $2/3$ up to a constant |
| equidistribution of small points | inapplicable (bounded degree, large height) or energy | Prop 3.3(a) |
| Zhang successive minima, Northcott, arithmetic bigness | geometric bigness of $aH-bE$ | $\tau^*<2$; item 3 (Prop 3.3(b)) |
| arithmetic Bogomolov on tori | inapplicable; Dirichlet exponent | Prop 3.3(c) |
| arithmetic Riemann–Roch / Bost slopes on $\mathbb P^1$ | Vandermonde, product formula | pair (Lemma 3.4) |
| slopes / bigness with finite-place vanishing on $X_M$ | the determinant method | item 6 |
| heights on $\overline M_{0,M}$ from local data + product formulas | pair product formulas (Theorem C) | satisfied by pseudo-adelic fakes |
| finite multi-point contacts $\Gamma^{>M}_T$, $\lvert T\rvert\ge3$ | precise multi-point term | zero on all record clusters and on fakes' completions (evidence); does not separate |
| global coherence (Vojta $D=0$ on $\overline M_{0,8}$ or $\mathrm{Bl}_YX_7$, $(L_8)$) | the record's conditional lane | conditional; open |

**What is new (all proved unless marked).**
1. *Theorem A* (exactly checked on 26 configurations): the exponent profile is the set of components of the regular model
   $\mathcal X_N$ of the conic met by the cluster sections; the reduced chords are their finite intersection numbers; the
   cut slack is the Néron–Hriljac vertical correction of $D_0=\sum P_j-\frac M2(O_u+O_v)$; the cut identity is the genus-0
   Hodge-index equality written place by place.
2. *Theorem B*: the $(2/3+\varepsilon)$ finite inequality is the adelic $M$-point (Fekete) diameter inequality; the chain
   lemma computes the $p$-adic $M$-point diameter of an $A_{e-1}$ chain. *Remark B.1* (numerics, constant only): the
   archimedean capacity would lower $K_0$ by about $4\log4$.
3. *Prop B.2*: the adelic-capacity class, which contains every genus-0 Arakelov inequality, is capped by the refereed
   row-2 fakes.
4. *Prop 3.3(b)*: no Zhang/Northcott/arithmetic-bigness argument on a fixed model in the $(H,E)$-plane reaches the arc,
   because it would need $\tau^*(M)>2$.
5. *Theorem C*: every height relation on $\overline M_{0,M}$ obtainable from local data and product formulas is a
   combination of pair product formulas, so it holds for pseudo-adelic completions of fakes; the uncapped content of
   Vojta/$(L_8)$ is global coherence of the adelic configuration (rationality), which intersection numbers do not see.
6. *Evidence*: on the record clusters the adelic slack sits at finite places (large primes carry $0$–$85\%$), the
   clusters are nearly Fekete on their own chord, and multi-point contacts at primes $>M$ are absent (1 triple in 15
   windows, at a non-extremal window; random-model expectation $0.02$).

**Statement needed** (not supplied by this lens): an inequality that uses the global coherence of the local
configurations across places in a form weaker than "lattice gaps at every prime" (which is the problem itself) and not
via a fixed-$S$ exceptional set. The only standard statement of that type in the record remains Vojta's $D=0$
inequality on $\mathrm{Bl}_YX_7$ or $\overline M_{0,8}$, i.e. $(L_8)$.
