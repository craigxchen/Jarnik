# Lattice points on short arcs of circles: what the record proves, what it does not, and where it stops

> **Abstract.** Cilleruelo–Granville conjecture that an arc of length $C\sqrt R$ on $x^2+y^2=R^2$ holds at most $M(C)<\infty$ lattice points, uniformly in $R$.
> *Growth rate:* the best unconditional bound is $M\le(2/3+\varepsilon)\log R/\log\log R$ (proved, two referees); the rate $\log R/\log\log R$ is unchanged.
> No improvement short of the uniform bound is known, unconditionally or from any standard conjecture.
> *Uniform bound:* open; it follows **conditionally** from one instance of Vojta's conjecture ($D=0$) on one fixed variety ($\mathrm{Bl}_YX_7$, or $\overline M_{0,8}$ via a lemma $(L_8)$).
> Restricted unconditional results: pure Hadamard-core clusters are uniformly bounded; one split prime of bounded cotangent valuation gives $M=O(\sqrt{\log R})$.
> Lower bounds: explicit clusters give $M(8.2394)\ge10$ and families give $M_\infty(C)\ge8$ for $C>5.89771$; computations to $R^2\le10^{12}$ are evidence only.
> Barriers: most known inputs are capped at $\log R/\log\log R$, mainly by *fake circles* (genuine primes and profiles, real angles) at $W=O(M\log M)$.
> One cap (the Subspace route QS) rests on unrefereed propositions; P3-all$^*$, H$^{\rm gen}$ and several analytic inputs remain uncapped.
> Route 5 (determinant method) reduces growth to a box inequality whose decisive fake test is undecided; other routes stop at precisely stated open lemmas.
> §7 states the smallest open problems; every statement carries a status label and a file pointer, and §8.2 lists the referee corrections.

*A consolidated account of a research program on the Cilleruelo–Granville conjecture, written for
experts in Diophantine approximation and analytic number theory. Compiled 2026-10-28 from the
review record in `docs/claude_review/` of the repository `craigxchen/Jarnik` and the research notes
on the branch `agent/uniform-gcd-formalization`.*

**How to read this document.** Every statement below is the version a referee accepted. Where a
referee corrected a statement, the corrected one is given and the correction is listed in §8.2.
Each statement carries a label:

| label | meaning |
|---|---|
| **proved** | proved with all quantifiers, and checked by at least one independent adversarial referee |
| **conditional** | proved assuming a named hypothesis, which is stated with it |
| **model** | proved, but in a probability or relaxation model, not for actual circles |
| **exact** | a finite exact computation, reproduced independently |
| **evidence** | a finite computation that no proof uses |
| **heuristic** | no proof |
| **unrefereed** | checked by no one but its author |
| **open** | neither proved nor refuted |

File pointers are relative to `docs/claude_review/` unless they say otherwise. "Item `N`" and bare
note names (such as `inert_prime_cofactor_bound.md`) refer to `docs/` on the branch
`origin/agent/uniform-gcd-formalization`.

**Summary of the state of the problem.**

1. *Growth rate.* The best unconditional bound is $M\le(2/3+\varepsilon)\log R/\log\log R$ (§2).
   **No growth improvement short of the uniform bound is known, either unconditionally or from a
   standard conjecture**: no intermediate rate $o(\log R/\log\log R)$ is derived from any standard
   conjecture. The only standard-conjecture route (item 2) gives the uniform bound directly
   (`README.md` growth-rate status; `round7/README.md` §0; `round8/README.md` §2.3(a)). Non-standard
   hypotheses that would give intermediate rates are listed in §6.8.
2. *Uniform bound.* It is open unconditionally. The only route from a standard conjecture is Vojta's
   Main Conjecture with $D=0$ for **one fixed variety**: the blow-up of the equal-norm variety $X_7$
   along its diagonal line, or $\overline{M}_{0,8}$ through an arithmetic lemma $(L_8)$ (§4;
   **conditional**).
3. *Restricted results.* Pure Hadamard-core clusters are uniformly bounded. One split prime with
   bounded cotangent valuation forces $M=O(\sqrt{\log R})$. Any five lattice points have diameter
   $\ge 1.9433\,R^{2/5}$ (§3).
4. *Barriers.* Most known inputs are capped at $\log R/\log\log R$. The usual witnesses are *fake
   circles*: angle models with genuine primes, a genuine exponent profile and real angles, which
   satisfy every inequality of the input class at $\log R=O(M\log M)$; a few caps are instead
   theorems on a variety, and one (the Subspace route QS) rests on unrefereed propositions (§5.4).
5. *Routes.* Several routes were stopped at a precise missing statement (§6). The smallest
   self-contained open problems are collected in §7.

---

## 1. Problem, history, notation, known bounds

### 1.1 The problem

Let $R>0$ with $N:=R^2\in\mathbb N$. Let $\mathcal E_R=\{(x,y)\in\mathbb Z^2: x^2+y^2=N\}$, and
identify $(x,y)$ with $z=x+iy\in\mathbb Z[i]$. For $C>0$ put

$$M(C):=\sup_{R}\ \max_{I}\ \#(I\cap\mathcal E_R),\qquad M_\infty(C):=\limsup_{R\to\infty}\ \max_{I}\ \#(I\cap\mathcal E_R),$$

where $I$ ranges over the arcs of $x^2+y^2=R^2$ of length at most $C\sqrt R$, in **any position**.
Since an arc of length $L$ holds at most $L+1$ lattice points, $M(C)<\infty$ iff $M_\infty(C)<\infty$.

> **Conjecture (Cilleruelo–Granville).** For every $C>0$, $M(C)<\infty$.

Cilleruelo and Granville state this for arcs of length $R^{1-\varepsilon}$ [CG07, Conj. 14]. In
coordinates they state it as boundedly many representations $N=a^2+b^2$ with $|b|$ in a window of
length $N^{\alpha}$, $\alpha<1/2$ [CG07, Conj. 13]. They argue that these forms and the diagonal
form (Conj. 15) are equivalent [CG07, p. 11; card in `lit18/gh/cg_index.md`]. This document is
about the critical case: arcs of length $C\sqrt R$.

**Two questions.** The document is organised around two questions:

* **the growth rate:** the best function $f(R)$ with $\max_I\#(I\cap\mathcal E_R)\ll_C f(R)$;
* **the uniform bound:** whether $f$ can be taken bounded.

The leading constant in front of $\log R/\log\log R$ is secondary. It is recorded, but no route is
judged by it.

**Why $\sqrt R$ is critical.** An arc of length $L$ on a circle of radius $R$ has sagitta about
$L^2/(8R)$. At $L\asymp\sqrt R$ the sagitta is of order $1$, so this is the first scale at which
the curvature of the circle is visible on the lattice. Below it, Liouville-type arguments work.
Arcs of length $C\sqrt R$ that touch one of the four axis points are trivial: if $|b|\le\sqrt R$
then $a=\sqrt{N-b^2}$ lies in an interval of length at most $1$ [CG07]. The difficulty is
uniformity over the **position** of the arc.

### 1.2 Notation used throughout

* $N=R^2$, $W:=\log N=2\log R$. An arc of length $\le C\sqrt R$ has angular width
  $\Delta=Ce^{-W/4}$.
* A *cluster* or *configuration* is a set of $M$ distinct $z_j\in\mathbb Z[i]$ with $|z_j|^2=N$ on
  one such arc. The chord bound $|z_i-z_j|\le C\sqrt R$ is the only geometric input of the main
  theorem.
* $p$ denotes a prime $\equiv1\pmod 4$ ("split"), $\ell$ a prime $\equiv3\pmod4$ ("inert"), and $q$ an
  arbitrary prime. For each split $p$ fix $\pi_p$ with $\pi_p\bar\pi_p=p$.
* *Normalised cluster*: one with $\gcd(z_1,\dots,z_M)=1$ in $\mathbb Z[i]$. Then $N$ is odd and
  has no inert factor, $N=\prod_p p^{e_p}$, and
  $z_j=\varepsilon_j\prod_p\pi_p^{a_{jp}}\bar\pi_p^{\,e_p-a_{jp}}$ with units $\varepsilon_j$
  (`two_thirds/proof.md` Lemma 1.1). The vectors $a_j=(a_{jp})_p$ form the **exponent profile**.
* $w(v)=\sum_j|v_j|\log p_j$ for $v\in\mathbb Z^r$; $d_{xy}=w(a_x-a_y)$; the **pair margin** is
  $\mu_{xy}=d_{xy}-W/2$.
* *Growth dictionary* (`round7/README.md`, Dictionary). For near-extremal clusters
  $\log M=(1+o(1))\log\log R$, so $W=cM\log M$ corresponds to
  $M=(2/c+o(1))\log R/\log\log R$. **A growth improvement means $W/(M\log M)\to\infty$ along
  near-extremal actual clusters.** An input class is **capped** when a barrier result shows that no
  argument using only those inputs proves $M=o(\log R/\log\log R)$. Most caps are *fake caps*: the class
  contains fakes (§5.1) with $W\le cM\log M$ for some fixed $c$. The other kinds of cap, and the status
  of each, are listed in §5.4.

### 1.3 History and known bounds

The literature citations in this subsection are taken from the repository paper
`paper/jarnik_sublog.tex` (§related, main branch) and from the 2015–2026 literature sweep `lit18/`.
They were not re-derived in this program.

1. **Jarník (1926).** A strictly convex arc of length $L$ passes through $\ll L^{2/3}$ lattice
   points. For $L=C\sqrt R$ this gives $\ll_C R^{1/3}$. It uses the curvature, not the arithmetic.
   The record notes that the minimal constants $C_3$ found for three points decay like $n^{-1/12}$,
   which is Jarník's exponent (`growth/families.md` line 368, evidence).
2. **Divisor bound.** $\#\mathcal E_R=r_2(N)=R^{o(1)}$ for the whole circle. For arcs of length
   $R^\alpha$ with $\alpha>1/2$ this is still the best known uniform bound.
3. **Cilleruelo (1991).** An arc of length $(16R)^{1/3}$ contains at most two lattice points. This
   is sharp: the triples $(4n^3-1,2n^2+2n)$, $(4n^3,2n^2+1)$, $(4n^3+1,2n^2-2n)$ realise it.
4. **Cilleruelo–Córdoba (1992).** For every $k\ge1$, an arc of length
   $\sqrt2\,R^{1/2-1/(4\lfloor k/2\rfloor+2)}$ contains at most $k$ lattice points [CC92; CG07
   Thm 13]. The exponent is sharp for $k\le3$ and has not been improved for any $k\ge4$ since 1992.
   Taking $k\approx\frac12\log R$ and cutting an arc of length $C\sqrt R$ into $O(C)$ pieces gives
   $\ll_C\log R$ at the critical scale.
5. **Cilleruelo–Granville (2007, 2009).** The conjectures above, the link to squares in arithmetic
   progressions (arcs $\Rightarrow$ Rudin-type statements, not the reverse), the sharp constant for
   three points at scale $R^{1/3}$, and infinite families of four points at that scale (one built
   from Fibonacci numbers) [CG07, CG09]. According to the project notes
   (`endpoint_continuation_status.md`, branch), their open question at the next step concerns **five**
   points on arcs of length $O(R^{2/5})$ (see Corollary 3.8). This attribution is **unverified against
   the CG papers**: the referee could not fetch them, and reports only that the reading is consistent
   with the four-point families at scale $R^{1/3}$ (`growth/referee_cotangent.md` §4). (The repository paper's line saying that four
   points on arcs $o(R^{2/5})$ are open is contradicted by the Pell family at scale $R^{1/3}$; §8.2.)
6. **Bourgain–Rudnick (2011–2013).** $O(\log m)$ points on arcs of length $\le m^{1/4}$ of the circle of
   radius $\sqrt m$ (so arcs of length $\le\sqrt R$ in our notation, with $m=N$), in the form "known
   for almost all $\lambda$, individually only $B_\lambda\ll\log\lambda$". Outside a
   sparse exceptional set of radii, distinct lattice points are $R^{1-\varepsilon}$ apart. **The best
   bound for every $R$ in the 2015–2026 literature is $O(\log R)$** (`lit18/README.md` §0, §1.1,
   row 1; status W: confirmed by an independent checker's snippet).
7. **Almost-all results.** Kurlberg–Lester (Invent. Math. 2025: Poisson spacing along a density-one
   subsequence), Granville–Wigman, Sartori, Klurman–Sartori. They are silent on the density-zero set
   of $N$ that carries large clusters (`lit18/README.md` §1.1).
8. **Withdrawn claim.** Oganesyan, arXiv:2107.09991, claimed unboundedness on arcs $R^\alpha$,
   $\alpha\in(1/2,1)$. It is withdrawn ("mistake in Lemma 2"). A correct version would have
   contradicted the Vojta-conditional theorem of §4.1, which covers arcs $CR^\alpha$, $\alpha<1$.
9. **Near-axis results.** Chan (IJNT 2015, Thm 3: at most 10 points with $|b|<n^{1/4}(\log n)^{1/7}$
   for $n=N^2$) and Erdős–Rosenfeld (divisors near $\sqrt n$). Both are bounded only near a direction
   of bounded height. After a Gaussian rotation they miss general position by a factor about
   $R^{1/2}$ (`lit18/direct.md` §2.3).

**Bounds at the scale $C\sqrt R$ in this program's lineage.**

| bound on $M$ | source | status |
|---|---|---|
| $\ll_C\log R$ | Cilleruelo–Córdoba 1992; Bourgain–Rudnick | literature |
| $K\log R/\log\log R$ for every $K>8$, $R\ge R_0(C,K)$ | `paper/jarnik_sublog.tex` Thm (geometric arc form), main branch | paper; accompanied by a Lean development (`GaussianChain/`), not audited in this review |
| $(4+\varepsilon)\log R/\log\log R$ | "the previous prose leading constant", as cited in `inert_prime_cofactor_bound.md` (branch) | cited only; **not audited** by the review |
| $(1+\varepsilon)\log R/\log\log R$ | `inert_prime_cofactor_bound.md` (branch), inequalities (6) and (10) | verified by the review (`README.md` §1) |
| $(2/3+\varepsilon)\log R/\log\log R$ | `two_thirds/proof.md` Cor. 8.1 | **proved**, two referees |

### 1.4 The gap

Between $\log R/\log\log R$ and $O(1)$ **nothing is proved for every $R$**. The data and heuristics
point to boundedness:

* **Exhaustive scan, $R^2\le10^{12}$** (`growth/families.md`; **evidence**). The author's computation
  gives: at most 3 points lie on arcs of length $\frac12\sqrt R$, and at most 4 on arcs of length
  $\sqrt R$; the best constants $C$ for $k=3,\dots,8$ points are
  $0.2524,\ 0.5060,\ 1.9962,\ 4.0001,\ 5.2699,\ 5.8978$. The referee did **not** re-run the $10^{12}$
  scan ("which would take about a day"); its independent brute force reproduced the data only for
  $n\le10^{10}$, and it found that the write-up's rows labelled "1e8" and "1e10" are really $n<2^{26}$
  and $n<2^{33}$ (`growth/referee_families.md` §5 and addendum).
* **Lower bounds for $M(C)$** (proved by explicit clusters and families; `growth/families.md`,
  `growth/referee_families.md` addendum and §8). Single clusters: $M(C)\ge7$ for $C\ge5.2700$,
  $M(C)\ge8$ for $C\ge5.6740$, $M(C)\ge9$ for $C\ge7.8804$, and $M(C)\ge10$ for $C\ge8.2394$. The last
  two come from a primitive **10-point** cluster found by the referee (exact constant
  $C=8.239326866\ldots$) on $n=1176852625=5^3\,13^2\,17\,29\,113$, and from the 9-point sub-cluster
  obtained by dropping one end point ($C=7.880379573\ldots$). Single clusters bound $M(C)$, not $M_\infty(C)$. Families: the golden-unit template gives
  $M_\infty(C)\ge8$ for $C>5.89771$ and $M_\infty(C)\ge7$ for $C>5.5365$; the axis Pell family gives
  $M_\infty(C)\ge6$ (hence $M(C)\ge6$) for $C>4$.
* **Heuristics** predict $M(1/2)=M(1)=5$, and "a uniform bound is plausible" (`growth/families.md`;
  heuristic).

## 2. The main unconditional theorem: $M \le (2/3+\varepsilon)\log R/\log\log R$

Source: `two_thirds/proof.md` (Theorems 5.1, 5.2, 7.5 and Corollary 8.1). Status: **proved**. Two
independent adversarial referees re-derived every step. They found no mathematical error, only the
six minor points listed in the file's review box and in §8.2 below.

### 2.1 Statements

> **Theorem 2.1 (effective inequality; `two_thirds/proof.md` Thm 7.5).** Let $N\ge1$, $R=\sqrt N$,
> $W=\log N$ and $C>0$. Every set of $M\ge2$ distinct lattice points on an arc of length
> $\le C\sqrt R$ of $x^2+y^2=N$, in any position, satisfies
> $$3M\log M\ \le\ W'+(K_0+4\log^+C)\,M\ \le\ \log(R^2)+(14+4\log^+C)\,M,$$
> where $W'=\log N'$ is the norm after removing the common Gaussian gcd, and
> $K_0=3\beta+\xi+2(2a_3+a_1)+6\log2<13.9092$.

> **Corollary 2.2 (`two_thirds/proof.md` Cor. 8.1).** For every $C>0$ and $\varepsilon>0$ there is an
> explicit $R_0(C,\varepsilon)$ such that for all $R\ge R_0$ with $R^2\in\mathbb N$, every arc of
> $x^2+y^2=R^2$ of length at most $C\sqrt R$ contains at most
> $$\Big(\frac23+\varepsilon\Big)\frac{\log R}{\log\log R}$$
> lattice points. Qualitatively, $M\le\frac23\frac{\log R}{\log\log R}\big(1+O_C(\frac{\log\log\log R}{\log\log R})\big)$.

$R_0$ depends only on $(C,\varepsilon)$; the proof uses only the chord bound, so it is uniform in
the position of the arc and in the factorisation of $R^2$. For $C\le1$ one may take
$R_0=\exp\exp\ell_0$ with $\ell_0=24.9,\ 71.7,\ 795.2$ for $\varepsilon=\frac13,\ 0.1,\ 0.01$.

**What this does not do** (`two_thirds/proof.md` §10). It changes only the constant. The growth rate
is still $\log R/\log\log R$: the method does not exclude configurations with $M\asymp\log R/\log\log R$
(up to $(2/3+\varepsilon)\log R/\log\log R$), because its total arithmetic lower bound is about
$\frac38M^2\log M$ (inert plus split) against the diagonal cut slack $\frac M4\log R$ (the right side of
(UB) below). The
constant $2/3$ is the limit of these ingredients. The chain constant is sharp, the pigeonhole bounds
are sharp up to $O(M)$ per prime, and deeper prime powers add only $O(M^2)$ in total. With every
small split prime dividing $N$ and levels following the optimal chain profile, the arithmetic lower
bound is $\frac38M^2\log M+O(M^2)$; this needs the easy upper bound
$F^*(M)\le\frac38M^2\log M+O(M^2)$ (take $e=1$, $n_0=n_1=M/2$), which a referee supplied. Whether
actual configurations realise this worst case was not examined.

### 2.2 Proof sketch

The proof compares one archimedean upper bound with one arithmetic lower bound for the product of
the **reduced chords**.

**Step 1: normalisation (Lemma 1.1).** Divide by $g=\gcd(z_1,\dots,z_M)$. The points stay distinct,
$R'=R/|g|\le R$, and the arc length satisfies $L'\le C\sqrt{R'}$. Afterwards $N$ is odd with no inert
factor, $z_j=\varepsilon_j\prod_p\pi_p^{a_{jp}}\bar\pi_p^{e_p-a_{jp}}$, and for each $p\mid N$,
$\min_ja_{jp}=0$ and $\max_ja_{jp}=e_p$. The units are arbitrary. No common unit class is imposed.

**Step 2: pair factorisation and the cut identity (§2).** Put $g_{ij}=\gcd(z_i,z_j)$ and
$c_{ij}=(z_i-z_j)/g_{ij}\in\mathbb Z[i]\setminus\{0\}$. Then
$\log|g_{ij}|^2=W-d_{ij}$ with $d_{ij}=\sum_p|a_{ip}-a_{jp}|\log p$. With
$S_{p\tau}=\#\{i:a_{ip}<\tau\}$,
$$\sum_{i<j}|a_{ip}-a_{jp}|=\sum_{\tau=1}^{e_p}S_{p\tau}(M-S_{p\tau})=\frac{e_pM^2}4-\sum_{\tau=1}^{e_p}\Big(S_{p\tau}-\frac M2\Big)^2.$$
Put $Q=\sum_{p\mid N}\log p\sum_\tau(S_{p\tau}-M/2)^2\ge0$ (the **cut slack**). Then exactly
$$\sum_{i<j}\log|c_{ij}|+\tfrac12Q=\sum_{i<j}\log|z_i-z_j|-\frac{M(M-2)}4\log R,$$
equivalently
$\prod_{i<j}|z_i-z_j|^8=N^{M(M-2)}\prod_{p\mid N}p^{\sum_\tau(2S_{p\tau}-M)^2}\prod_{i<j}|c_{ij}|^8$.
This identity was verified as a big-integer equality on 6,176 configurations. With
$|z_i-z_j|\le C\sqrt R$:
$$\sum_{i<j}\log|c_{ij}|+\tfrac12Q\ \le\ \frac M4\log R+\binom M2\log C. \tag{UB}$$

**Step 3: arithmetic lower bound (§3).** Collisions at distinct primes multiply:
$\log|c|\ge\frac12\log2\cdot\mathbf 1[(1+i)\mid c]+\sum_{q\text{ odd}}v_q(c)\log q$. The conic count
$\#\{x^2+y^2\equiv n\bmod q^a\}=(q-\chi_4(q))q^{a-1}$ for $q\nmid n$ gives the collision criteria:

* $(1+i)\mid c_{ij}$ for every pair;
* inert $\ell$: $\ell^a\mid c_{ij}$ iff $z_i\equiv z_j\pmod{\ell^a}$, among $(\ell+1)\ell^{a-1}$ classes;
* split $p\nmid N$: the same, among $(p-1)p^{a-1}$ classes;
* split $p\mid N$: only pairs on the **same exponent level** $a_{ip}=a_{jp}$ can collide. They do
  so iff the cofactors $r_i=z_i/(\pi^{a_{ip}}\bar\pi^{e-a_{ip}})$ agree mod $p^a$, among
  $(p-1)p^{a-1}$ classes per level. Different levels never collide at $p$.

With $E(n,Q)=\min\sum_r\binom{m_r}2\ge n^2/(2Q)-n/2$ (pigeonhole) and level sizes
$n_{pt}=\#\{i:a_{ip}=t\}$, this gives a lower bound for $\sum\log|c_{ij}|$.

**Step 4: the chain lemma (§4).** This is the new ingredient. At a split $p\mid N$ the points spread
over levels $0,\dots,e$. Few points per level means few collisions, but it also forces the cut
counts $S_\tau$ away from $M/2$, and so a large cut slack.

> **Lemma 2.3 (chain lemma; Lemma 4.1).** Let $\lambda>0$ and $c(\lambda)=\frac{\sqrt{1+4\lambda}-1}2$,
> the positive root of $c^2+c=\lambda$. For any distribution $(n_0,\dots,n_e)$ of $M$ points into
> levels, with $S_\tau=n_0+\dots+n_{\tau-1}$,
> $$\sum_{\tau=1}^{e}\Big(S_\tau-\frac M2\Big)^2+\lambda\sum_{t=0}^{e}n_t^2\ \ge\ c(\lambda)\,\frac{M^2}2 .$$
> The constant is sharp ($\inf_e\min\Phi=cM^2/2$), and the bound holds for all real profiles.

*Proof idea.* Put $x_\tau=S_\tau-M/2$, so $x_0=-M/2$, $x_{e+1}=M/2$ and $n_t=x_{t+1}-x_t$. For one-sided
chains $F_n=\sum_{k\ge1}y_k^2+\lambda\sum(y_k-y_{k+1})^2+\lambda y_n^2$, induction gives
$F_n\ge\gamma_ny_0^2$ with $\gamma_0=\lambda$ and $\gamma_n=T(\gamma_{n-1})$,
$T(\gamma)=\lambda(1+\gamma)/(1+\gamma+\lambda)$. $T$ is increasing with fixed point $c$, so
$\gamma_n\ge c$. Split at the last index with $x_m\le0$; there
$(x_{m+1}-x_m)^2\ge x_m^2+x_{m+1}^2$, which yields two one-sided chains of height $M/2$ each. $\square$

For $\lambda=1/(p-1)$, $c_p=c(\lambda_p)\ge\lambda_p-\lambda_p^2$. So at each split prime, whether or
not it divides $N$, the combined contribution (collisions plus half the cut slack) is at least
$(c_pM^2/4-M/2)\log p$ (Lemma 4.3). Against the squarefree model this loses only a factor
$c/\lambda=1-O(\lambda)$, uniformly in the exponent pattern. It removes the "dilution by exponent
levels" that an earlier note had claimed was an obstruction (a correction, §8.2).

**Step 5: the exact finite inequality (Theorems 5.1, 5.2).** Combining Steps 2–4,

$$F^*(M):=\frac{K\log2}2+\sum_\ell\log\ell\sum_{a\ge1}E\big(M,(\ell+1)\ell^{a-1}\big)+\sum_pG_p(M)\log p\ \le\ \frac M4\log R+\binom M2\log C,$$

where $K=\binom M2$ and $G_p(M)=\min\{\sum_aE(M,(p-1)p^{a-1}),\Psi_p(M)\}$, with $\Psi_p$ the minimum
of the level-collision-plus-cut bracket over all level distributions. $F^*(M)$ is finite and exactly
computable: only $\ell\le M-2$ and $p\le M$ contribute (corrected from "$\ell<M$"; `two_thirds/proof.md`
review box item (v)), and the minimum in $\Psi_p$ is attained with
$e\le M-1$ and no empty level. Dropping the split-prime terms recovers inequality (6) of
`inert_prime_cofactor_bound.md`. Numerically (evidence) $F^*(200)=49482.5$, so for $C=1$ at least 200
points need $\log R\ge989.7$.

**Step 6: prime sums (§§6–7).** Keep the ramified term, inert $\ell\le M$ at depth one, and split
$p\le M$ via Lemma 4.3. Multiplying by $8/M$ gives
$M(4A(M)+2B(M))-4\theta_{\rm odd}(M)+2(M-1)\log2\le W'+4(M-1)\log C$, with
$A(x)=\sum_{\ell\le x}\frac{\log\ell}{\ell+1}$ and $B(x)=\sum_{p\le x}c_p\log p$. The weight is
$2$ on inert primes and $1$ on split primes, so the leading term is
$M(2S_3+S_1)\cdot 2\approx 3M\log M$. Only one-sided explicit bounds are needed:

* $\theta(y)<y\log4$ and $\psi(y)<y(\log4+2/e)$ (Chebyshev, proved in full);
* $\sum_{3\le q\le n}\frac{\log q}q\ge\log n-\beta+\frac1n$ with $\beta=1+\sigma+\frac12\log2$
  (Legendre's formula);
* $\sum_{q\le y}\chi_4(q)\frac{\log q}q\le\frac4\pi\psi(y)/y<\xi=2.70189$ (Leibniz series and
  an alternating-series sign argument);
* $\sigma\in[0.7553665,0.7553683]$, $a_3\le0.1743774$, $a_1\le0.0225136$. This is the only numerical
  input: partial sums to $10^7$ with explicit tails, recomputed by a referee with correctly rounded
  34-digit logarithms.

This gives Theorem 2.1. Corollary 2.2 follows by monotonicity of $h(M)=M(3\log M-K_C)$.

**Where the $\log\log R$ comes from.** In the inequality $3M\log M\le W+O(M)$, the $\log M$ is the
mass $\sum_{q\le M}\log q/q$ of primes up to the cluster size, at which pigeonhole forces
collisions. Counting and pigeonhole force collisions only at moduli below about $M$; nothing in this
proof forces them at primes much larger than $M$ (interpretation; compare `lit18/README.md` §1.3 on
the combinatorial large sieve, which has the same mechanism and the same rate). This is the source of
the barrier described in §5.

## 3. Restricted unconditional theorems

These are proved and refereed (one item in §3.4 is conditional on abc and is labelled). Each either
bounds $M$ for actual circles under a **restriction** on the exponent profile, the norm or the scale,
or is a structural statement that does not bound $M$ by itself. None extends to general circles: the
known fake circles (§5) violate each restriction, and no restriction is known to follow from being a
large cluster. They show which structures are already excluded.

### 3.1 Hadamard and Walsh profiles (`growth/walsh-extraction.md`; referee: survives)

*Setting* (§§1, 4). After removing the Gaussian gcd $g$, write
$z_x=g\,\varepsilon_x\prod_a G_a^{H(x,a)}$, $x=1,\dots,M$, where $H$ is a normalised Hadamard matrix of
order $M$ (first column all ones, $H^tH=MI$), $a$ runs over the nonconstant columns, and the blocks
$G_a$ are conjugate-primitive Gaussian integers (possibly $1$) with pairwise disjoint rational prime
supports; $G^{-1}$ means $\bar G$. This is a **pure Hadamard-core cluster**. It covers the pure
(unflipped) repeated-Walsh (Sylvester) and Paley profiles with any number of copies per label.

> **Theorem 3.1 (Theorem B).** Let $M\ge28$ be a Hadamard order and
> $s_M=\lfloor(1+\sqrt{4M-3})/2\rfloor+8$. A pure Hadamard-core cluster of order $M$ on an arc of
> length $\le C\sqrt R$ satisfies
> $$W\le\max\Big(\frac{8Ms_M\log^+C}{M-2s_M},\ 36\log^+(2C)\Big),$$
> where $W=\sum_j\log p_j$ is the log of the primitive norm.
>
> **Corollary 3.2 (B.2, B.3).** Every pure Hadamard-core cluster on an arc of length $\le\sqrt R$ has
> $M\le24$. For every $C$, pure Hadamard-core clusters have $M<\max(28,M_1(C))$ with
> $M_1(C)=C^{2+o(1)}$. Effective and elementary.

The referee added the hypothesis $M\ge28$ (i.e. $M>2s_M$), which the original summary dropped, and the
definition of $M_1(C)$ as the threshold beyond which the inequality holds for all larger orders.

*Proof idea.* The arc condition pins each block angle:
$\arg G_a\in\frac{\pi}{2M}m_a+[-\Delta/2,\Delta/2]$ with integers $m_a=\sum_xH(x,a)k_x$ built from
the winding integers. Comparing **different** labels through their residues $m_a\bmod M$ shows that
the light labels form a Sidon set in $\mathbb Z/M$: a coincidence $r_a+r_b=r_c+r_d$ would produce a
conjugate-primitive Gaussian integer of small height within $\Delta$ of $(\pi/2)\mathbb Z$, which
the axis gap $|\beta|\sin\epsilon\ge1$ forbids. So there are at most $\sim\sqrt M$ light labels and
at most 8 heavy ones. On the other hand, a second-moment count (Lemma B.1) forces at least
$MW/(2W+8M\log^+C)$ labels to carry weight. The two bounds contradict each other. The cross-label
residue comparison is new (referee). It replaces the Roth/Matveev single-target estimate of earlier
notes.

*Comparison with earlier notes, as corrected by the referee* (`growth/referee_walsh-extraction.md` §§1, 6;
the write-up understates these items). Item 276 already proves an **effective uniform bound for pure
cores of any normalised Hadamard type, Paley included**, via an integral character and Matveev, but only
in the unit regimes $C\le\sqrt2$ (or $C\le2$ with a common unit). Item 278 covers pure Walsh cores for
every $C$. What Theorem B adds is every $C$ for non-Walsh cores, with small explicit constants. For
Theorem B'' below, item 277 already allowed about $\log_2M$ deleted rows in unit regimes.

> **Theorem 3.3 (Theorem B'').** Let $M\ge4096$ be a Hadamard order, and let $F$ be a set of
> $f\le M/(40\log M)$ rows such that the sub-cluster on the other rows is a pure Hadamard-core
> cluster restricted to those rows. If $W\ge\max(8M\log^+C,32\log^+(LC))$ with $L=\lfloor M/60\rfloor$
> (automatic for $M\ge C^{6+o(1)}$, per the referee's corrected statement), the configuration does not
> exist. So arbitrary modifications confined to at most $M/(40\log M)$ rows are excluded; this extends
> item 277's $\sim\log_2M$ deletions in unit regimes.

> **Theorem 3.4 (Theorem C; growth for a model class).** In the repeated-Walsh model WF$(t,b,F)$
> (rows $\mathbb F_2^t$, $M=2^t$, $b\ge5$ physical columns per nonzero label, an **arbitrary** set $F$
> of rows each flipped at one column, injectively), on arcs of length $\le C\sqrt R$:
> $W_0/(M\log M)\to\infty$, uniformly in $b$, $F$, prime sizes and units, where $W_0\le\log R^2$ is the
profile weight of item 333. Quantitatively $W_0\ge cM\log M\sqrt T/\log_2T$ and
> $M=O_C\big(\log R\,\log_2T_R/(\log\log R\sqrt{T_R})\big)$ with $T=\log_2^*M$. Also, uniformly in
> $R$, fewer than $2S_0^{1/S_0}M^{1-1/S_0}$ rows are unflipped, with $S_0=32$ for $C\le1$.

The referee notes that the growth part of Theorem C already follows from items 278 and 333 of the
notes. What Theorem C adds is that every flip pattern is covered.

**Where it stops** (§§7, 10). *Fully flipped* capacity-two prime-order Paley profiles on actual nearby
primes pass every pair inequality, every integral and rational character inequality, the aggregate
small-prime collisions and the copy capacity, at $W=\Theta(M\log M)$ (Prop P; most of it restates items
334, 391–393). The missing statement is the **fully flipped Hadamard phase lemma**: $M-1$ simultaneous
Diophantine approximations of Gaussian blocks, with targets that move with $M$ free prime angles. The
referee notes that the claim that this is *exactly* the barrier is not proved. Round 4
(`round4/flipped.md`) reduces each label pinning to a nonzero linear form in two logarithms,
$M\,\mathrm{Log}(G_a/\bar G_a)-2\,\mathrm{Log}(P_a/\bar P_a)+b\,\mathrm{Log}\,i$, with error at most
$M\Delta$. A Lang–Waldschmidt-type two-logarithm bound over $\mathbb Q(i)$ with naive heights would
give the lemma when the flipped mass is at most $(1/2-\delta)W$ (**conditional** Theorem L). Proved
Baker-type bounds are too weak. For the canonical Paley assignment the heavy-flip regime is vacuous
($W_F\le W/2+O(M\log C)$), so the residual case is the critical band $W_F/W\to1/2$.

**Conditional consequence for the known fakes** (`round8/README.md` §1.1, `round8/lw.md` Thm A item 4;
refereed; **conditional on LW**). Theorem L is correct, and under the Lang–Waldschmidt-type bound the
genuine versions of the `round5/construct.md`, `round5/lower.md` and `round5/round6_sharp/sharp.md`
(P3-char) fake profiles do not exist at Paley orders $M=q+1\ge M_0$ ($M_0$ ineffective). At other $M$
these fakes are row-deleted sub-profiles, which Theorem L does not cover. The skew double core (§5.2
item 1) escapes single-form LW altogether (§5.4 rows 11–12).

### 3.2 The cotangent trade-off (`growth/cotangent.md`; referee: survives, one major correction)

Let $L_0$ be the **intrinsic cotangent scale** of a normalised cluster: the lcm of all reduced
half-angle-cotangent denominators. For a split prime $p$, each exponent level holds at most
$c_p=(p-1)p^{v_p(L_0)}$ points (item 151).

> **Theorem 3.5 (Theorem A).** Every cluster of $M$ points on a primitive circle with all distances
> $\le CN^{1/4}$ satisfies
> $$\sum_{p\equiv1(4)}\log p\cdot\frac{(M-4c_p)_+^3}{12\,c_p}\ \le\ \frac M4W+\binom M2\log\frac{C^2}2,\qquad c_p=(p-1)p^{v_p(L_0)}.$$

The proof combines the cut identity with a cubic **level-spreading deficit**: with at most $c$ points per
level, $\sum_\tau(S_\tau-M/2)^2\ge(M-4c)_+^3/(12c)$. An exact shortest-path DP over all level vectors
found no violation and shows that the deficit lemma is asymptotically sharp.

**Corollaries, as strengthened by the referee.**
* If $v_p(L_0)\le B$ for **one** split prime $p$, then $M\le(\sqrt{6c_p/\log p}+o(1))\sqrt{\log R}$. For
  example $5\nmid L_0$ gives $M\le(3.86+o(1))\sqrt{\log R}$.
* If $L_0$ has no split prime factor $p\le M/8$, then
  $M\le(\sqrt{24}+o(1))\sqrt{\log R/\log\log R}$. (The write-up had $\sqrt{192}$.)

**Scope.** This does **not** improve the general rate. A cluster can put high powers of the small split
primes into $L_0$ at negligible cost, since one colliding pair per level suffices to raise $v_p(L_0)$.

> **Theorem 3.6 (Theorem B; rigidity at the Cilleruelo–Córdoba scale).** Let $\theta_M=(M-1)/(4M)$ for
> odd $M$ and $(M-2)/(4(M-1))$ for even $M$. If $M$ points of a primitive circle have all mutual
> distances $\le KN^{\theta_M}$, then the product of the reduced cotangent denominators, the cut
> deficit and $L_0$ are bounded by explicit functions of $K$ and $M$. For five points,
> $\sigma\prod_{i<j}s_{ij}\le K^{10}/128$ and $6\mid L_0$.
>
> **Corollary 3.7.** Any five lattice points on a circle of radius $R$ have diameter at least
> $1.9433\,R^{2/5}$.
>
> **Corollary 3.8.** The five-point question (can five lattice points lie on arcs of length
> $O(R^{2/5})$ for infinitely many $R$?) is equivalent to an integer clique problem with a **bounded**
> cotangent scale $L\in6\mathbb Z$. The referee checked both directions. (That this is
> Cilleruelo–Granville's question is taken from the project notes and was not verified against the CG
> papers; see §1.3 item 5.)

Calibration (referee): the plain Cilleruelo–Córdoba product with $s\ge1$ already gives
$2^{7/10}R^{2/5}=1.6245\,R^{2/5}$. The residue input $\prod s\ge6$ accounts only for the step to
$1.9433$. The referee rates Theorem B "new but modest". Exhaustive check: the diameter bound holds on
all 3,138,560 cyclic five-windows with $N\le10^6$.

**Evidence** (exact searches, no proof uses them). The best five-point clique with $L_0=6$ and
$A\le10^6$ has $C_*=1.187$. Over $L\in6\mathbb Z$, $L\le240$, $A\le2\cdot10^4$ (and $L\le120$,
$A\le10^5$) there is no five-point clique of diameter $<4.94\,R^{2/5}$, **in the span-at-most-$90°$ regime
that was searched**. The write-up's unrestricted "never below 4.9475" is false ($N=25$ gives 4.2024).

> **Theorem 3.9 (Theorem C, corrected; direction pinning for even $M$).** For even $M\ge4$ at the
> Cilleruelo–Córdoba scale, $\prod_iz_i=\varepsilon\,N_{\rm bal}^{M/2}\,m\,\Phi$ with
> $m=\prod p^{\min(h,M-h)}$ and $\Phi$ a Gaussian integer of bounded norm. So the cluster lies within
> $O_K(N^{\theta_M-1/2})$ of one of finitely many directions $\theta$ with $e^{2iM\theta}\in\mathbb Q(i)$.

The write-up's identity omitted the rational factor $m$ and fails on 34,164 of 100,256 even windows;
the corrected one holds on all of them. The mechanism is in the notes for exactly balanced clusters
(item 15).

### 3.3 Large prime factors, smoothness and rank (`round8/stability.md`, `round8/weakest.md`)

> **Theorem 3.10 (`round8/stability.md` Thm 3.1, Cor 3.2; referee: survives).** For $M\ge4$, every prime
> $q\mid N$ of a normalised cluster satisfies
> $\log q<2W/M-3\log(M/2)+K_0+4\log^+C$. Hence
> * $M<2W/(\log P_+(N)-4\log^+C)$, and $M=o(\log R/\log\log R)$ along any sequence with
>   $\log P_+(N)/\log\log R\to\infty$ ($C$ fixed);
> * $M<\max(4,2/\lambda)$ if $P_+(N)\ge N^\lambda C_+^4$.

The referee attributes the crude forms ($q<N^{2/M}C_+^4$ and these corollaries) to a two-line
consequence of Cilleruelo–Córdoba (1992) after dividing by a shared Gaussian factor. What is new is the
second-order $F^*$ gain. Every fake of size $O(M\log M)$ has $\log P_+=O(\log M)$, so the restriction
fails for them.

> **Theorem 3.11 (near-extremal structure; `round8/stability.md` §§2–5).** Assume
> $M\ge\varepsilon\log R/\log\log R$ and $R\ge R_0(\varepsilon,C)$. Then
> $P_+(N)\le M^{4/\varepsilon-3+o(1)}$; $p^e\Vert N'$ implies $p^e<C_+^8(\log R)^{8/\varepsilon}$;
> $(M/4)(1-o(1))<\omega(N')\le\Omega(N')\le(2/\varepsilon+o(1))M$; and the primes below $M^{1-\delta}$ carry
> an $O(M^{-\delta}/\log M)$ fraction of $W'$. Also, for every subset $B$ of size $s\ge2$,
> $\log\mathrm{Norm}\gcd_B<W/s-3\log s+K_0+4\log^+C$.

So near-extremal clusters live on $M^{O(1)}$-smooth norms with $\Theta(M)$ prime factors. The fakes
share this.

> **Theorem 3.12 (rank reduction; `round8/weakest.md` Thm 1.3).** If a normalised configuration lies on an
> arc of length $\le C\sqrt R$, its exponent lattice has rank $k$, and
> $\log N'\ge18(2\log^+C+\log2)$, then
> $$M\le\frac{(1+2W/\rho)^k+1}2\le\frac{5.5^k+1}2,\qquad\rho=W/2-2\log C-\log2.$$

So the uniform theorem is **equivalent to bounded rank** of the exponent lattice, and rank
$\le\theta\log W'$ with $\theta<1/\log5.5$ would give $M=O((\log R)^{1-\eta})$. This is a reduction only:
no rank bound is known, and the fakes have full rank $M-1$.

### 3.4 Other proved statements about actual clusters

* **abc, restricted** (`growth/conditional.md` Cor 8.2; **conditional on abc over $\mathbb Q$**). If
  $\mathrm{rad}(N)\le N^{1/4-\eta}$, there is at most one point per arc $C\sqrt R$ for
  $R\ge R_0(C,\eta)$. abc gives no growth improvement: every abc-type inequality is vacuous on
  squarefree near-extremal block relations, and all fakes are squarefree (`round8/weakest.md` Props 4.1–4.2).
* **Bisector lemma** (`round9/L1L2.md`; refereed): $M\le\sqrt2\,C^2|g_{xy}|+2$ for any pair.
* **Pairs outside one matching** (`round10/additive.md` Prop 3.1; refereed, capped class):
  $\mathrm{Norm}\,g_{xy}\ge N^{1/4}/(\sqrt2C)$ for every pair outside one matching.
* **Additive rigidity** (`round10/additive.md` Thms 2.1–2.6; refereed). If $M\ge5$ points of a circle
  satisfy all linear relations of an actual cluster, they are that cluster up to similarity. The
  relation lattice has rank $M-3$ iff the points lie in a lattice coset. Exact identity
  $\mathrm{covol}(K)^2I^2=\sum_{\rm triples}(2\,\mathrm{Area})^2$.
* **Forced divisibility is supercritical in the growth regime** (`round10/determinant.md` Thm 2.1,
  Thm 4.1; refereed). Once the auxiliary form is chosen **after** the profile, the forced-divisibility
  LP is supercritical by a factor $\ge0.12\sqrt{(n-1)/\ln4r}$. Consequently an actual cluster with
  $n\ge64$, $n-1\ge300\ln4r$ and $W\ge16\ln n+8\log^+C+14$ satisfies an **exact** algebraic system:
  $F_c(z)=0$ for every integral $c\in V(A,m)$ with $|c|_\infty\le H$. This corrects the earlier
  record, which called forced divisibility "exactly critical" (true only for a form fixed before the
  cluster). Almost every fake violates the system (Prop 4.2). The theorem needs $n\ge2137$ points, so
  it may be vacuous, and it gives no bound by itself; it is the starting point of route 5 (§6.2).

### 3.5 The constants actually achieved (lower bounds for $M(C)$)

See §1.4. The record proves $M(C)\ge7,8,9,10$ for $C\ge5.2700,\,5.6740,\,7.8804,\,8.2394$ (single explicit
clusters), and $M_\infty(C)\ge6,7,8$ for $C>4,\,5.5365,\,5.89771$ (families). Nothing in the record proves $M_\infty(C)=\infty$
for any $C$.

## 4. Conditional theorems (Vojta)

Both theorems in this section are **conditional**. Each assumes **one instance** of Vojta's Main
Conjecture (Vojta, LNM 1239, Conj. 3.4.3; Bombieri–Gubler Ch. 14) for **rational points**, with
$D=0$, on **one fixed** smooth projective variety over $\mathbb Q$. No set of places $S$ enters, so
no uniformity in $S$ is assumed; the exceptional set depends only on the variety, the big divisor
and $\varepsilon$. These are the only routes in the record from a standard conjecture, and they give the
**uniform bound** directly (not merely a growth improvement). Deriving approximation exponents from
Vojta on blow-ups is standard (Silverman, McKinnon, McKinnon–Roth); what is new is the variety and
the uniformity mechanism (`growth/README.md`, referee novelty check).

### 4.1 Vojta on the blow-up of $X_7$ (`growth/conditional.md` Thms 5.1–5.2; referee: survives as stated)

**The variety.** $X_M=\{x_1^2+y_1^2=\dots=x_M^2+y_M^2\}\subset\mathbb P^{2M-1}$. Over $\mathbb Q(i)$, with
$u_j=x_j+iy_j$, $v_j=x_j-iy_j$, it is the projective cone over $\{u_1v_1=\dots=u_Mv_M\}$: an irreducible,
normal, Gorenstein toric variety with canonical singularities, $\omega=\mathcal O(-2)$, singular exactly
where at least two pairs $(u_j,v_j)$ vanish (Lemma 2.1). Its **diagonal line** is
$Y=\{(x,y,x,y,\dots,x,y)\}$, which lies in the smooth locus. Let $\pi:\mathrm{Bl}_YX_M\to X_M$ and let
$X^\#_M$ be a smooth projective model. Then $K_{\mathrm{Bl}}=-2H+(M-2)E$.

**Hypothesis $V(M,\varepsilon)$.** There are a proper Zariski-closed $Z\subset X^\#_M$ and a constant $c$
with $h_{K_{X^\#}}(P)\le\varepsilon\,h_H(P)+c$ for all $P\in X^\#_M(\mathbb Q)\setminus Z$. ($H$ is the
pull-back of the hyperplane class; an ample $A$ may be used instead after shrinking $\varepsilon$.)

> **Theorem 4.1 (`growth/conditional.md` Thm 5.1).** Assume $V(7,1/4)$. Then there are an absolute integer
> $M_0$ and an ineffective constant $c_0\ge1$ such that for every $C>0$, every $R$ with
> $R^2\in\mathbb Z$ and $R>c_0^2C^{42}$, every arc of length $\le C\sqrt R$ on $x^2+y^2=R^2$, in any
> position, contains at most $M_0$ lattice points. Hence $M(C)\le\max(M_0,\,c_0C^{22}+1)<\infty$ for
> every $C$.

(`growth/README.md` records that one fixed $\varepsilon<1/2$ suffices for the argument; the theorem is
stated and refereed with $\varepsilon=1/4$.)

> **Theorem 4.2 (`growth/conditional.md` Thm 5.2; general exponent).** Let $\alpha\in[0,1)$,
> $M_\alpha=3+\lfloor2/(1-\alpha)\rfloor$ and $\varepsilon_\alpha=((M_\alpha-2)(1-\alpha)-2)/2>0$. Then
> $V(M_\alpha,\varepsilon_\alpha)$ implies: there is an absolute $M_0(\alpha)$ such that for every $C$
> and $R\ge R_1(\alpha,C)$, every arc of length $\le CR^\alpha$ has at most $M_0(\alpha)$ lattice points.
> In particular Vojta's conjecture implies the Cilleruelo–Granville conjecture for arcs $R^{1-\varepsilon}$.

**Proof sketch.** A short-arc $7$-tuple $P$ lifts uniquely to $P^\#\in X^\#_7$, off $E$ and off the
exceptional divisors $F_i$ of the resolution. Write $R'$ for the radius of the tuple after removing its
content. The arc condition says that $P$ is archimedean-close to $Y$, so
$h_E(P^\#)\ge\lambda_{E,\infty}(P^\#)-O(1)\ge\frac12\log R'-\log C-O(1)$ (Lemma 4.1 with $\alpha=1/2$;
Silverman's arithmetic distance functions). With $h_H(P)=\log R'+O(1)$,
$$h_K(P^\#)=-2h(P)+5h_E(P^\#)+\sum a_ih_{F_i}(P^\#)+O(1)\ \ge\ \tfrac12\log R'-5\log C-O(1).$$
If $P^\#\notin Z$, then $V(7,1/4)$ gives $h_K\le\frac14\log R'+O(1)$, so $R'\le c_0C^{20}$. For $R>c_0^2C^{42}$
one has $R'\ge\sqrt R/C>c_0C^{20}$, so every short-arc $7$-tuple lies in the fixed proper closed set
$W=\phi\pi(Z)\subset X_7$, which does not depend on $C$. **Bezout grid lemma** (Lemma 3.1, Prop. 3.2): if
a form $F_W$ vanishes on $W$ but not on $X_7$, at most $2\deg F_W+6$ points of the circle can lie on one
arc, since otherwise $F_W$ would vanish on a whole conic $Q_N$, and by homogeneity on a dense subset
of the irreducible cone. $\square$

**Unconditional by-products** (proved, refereed).
* *Non-density equivalence* (Prop. 3.2). The uniform theorem holds iff for each $C$ the set of short-arc
  $M$-tuples is not Zariski dense in $X_M$ for some $M$. One fixed variety serves every $N$
  (isotriviality). The absolute form is Roth-type finiteness at Dirichlet exponent $1/2$ for the diagonal
  line in a fixed $X_M$; the known families force $M\ge9$ there (`round4/hardness.md`).
* *Ru–Vojta with $S=\{\infty\}$* (Props. 6.1–6.2) reproves, ineffectively, the qualitative
  Cilleruelo–Córdoba theorem below exponent $1/2$, from $\beta(H,E)\ge2n/(n+1)^{1+1/n}\to2$.
* *This route cannot be made unconditional with these tools.* `threshold/upper.md` proves
  $\beta_M\le\tau^*(M)=2-1/\lceil M/2\rceil<2$ for every $M$ (§5.4, row 8), so Ru–Vojta with $S=\{\infty\}$ never
  reaches $\alpha=1/2$.

**Further conditional statements in the same note** (not unconditional by-products).
* *Hyperbola analogue* (Cor. 5.3; **conditional** on $V(7,1/4)$ for the split model
  $\{x_1y_1=\dots=x_7y_7\}$): the Erdős–Rosenfeld question has a positive answer, with an absolute $K$; and
  $V(3+\lfloor1/\varepsilon\rfloor,\cdot)$ implies Ruzsa's conjecture on divisors in
  $\sqrt n\pm n^{1/2-\varepsilon}$. Cor. 5.3 is only sketched in the note; `growth/referee_conditional.md`
  §1 supplies the details. That referee's literature check was **WebSearch only** (direct fetches of arXiv,
  UCL and UAM were blocked) and found no prior derivation of the $C\sqrt R$ arc bound, the CG
  $R^{1-\varepsilon}$ conjecture or Erdős–Rosenfeld/Ruzsa from Vojta, Bombieri–Lang or abc
  (`growth/referee_conditional.md` §5).
* *Other conjectures* (`growth/conditional.md` §§7–8; `round8/weakest.md` Prop 6.2). For Bombieri–Lang the
  note found nothing it could use (the auxiliary varieties are rational); the referee (minor item 4)
  accepts only this weaker wording, not "Bombieri–Lang gives no implication". The vacuity argument is
  `round8/weakest.md` Prop 6.2. abc gives only the restricted statement of §3.4.

### 4.2 Vojta on $\overline M_{0,8}$ and the lemma $(L_8)$ (`growth/ptolemy.md` Thm D, §7; referee: survives)

**Ptolemy data** (`growth/ptolemy.md` §1; item 532). Retain one Gaussian-unit class of a short-arc cluster and
divide by its Gaussian gcd. Order a $k$-tuple $(z_0,\dots,z_m)$, $k=m+1$, with anchor $z_0$. The anchored
numerators $P_i=\prod_p\pi_p^{(a_i-a_0)_+}\bar\pi_p^{(a_0-a_i)_+}=X_i+iY_i$ ($P_0=1$) are
conjugate-primitive with $z_i/z_0=P_i/\bar P_i$. Threshold layers give blocks $H_T$ ($\emptyset\ne T\subset[m]$)
with $n_T=\mathrm{Norm}(H_T)$, $P_i=\prod_{T\ni i}H_T$ and $\mathrm{Norm}\gcd(P_i,P_j)=\prod_{T\ni i,j}n_T$. At
each rational prime the blocks with $p\mid n_T$ form two nested chains on disjoint row sets. The
**residues** are $t_{ij}=\det(P_i,P_j)/\prod_{T\ni i,j}n_T\in\mathbb Z\setminus\{0\}$, $t_{0j}=Y_j$; the
**residual height** is $H_\Sigma=\sum_j\log|Y_j|+\sum_{i<j}\log|t_{ij}|$, and the scale is $w=W/2^m$.
The cross-ratios of the $P_i$ give a point $\Phi\in M_{0,k}(\mathbb Q)$.

**Hypothesis $V_8$.** Vojta's Main Conjecture for $X=\overline M_{0,8}$ over $\mathbb Q$ with $D=0$ and
$A=K_X+B$ ($B$ the boundary; ample by Keel–McKernan): for every $\varepsilon>0$ there are a proper
Zariski-closed $Z_\varepsilon$ and $C_\varepsilon$ with $h_{K_X}(P)\le\varepsilon h_A(P)+C_\varepsilon$ off $Z_\varepsilon$.

> **Theorem 4.3 (`growth/ptolemy.md` Thm D).** Assume $V_8$. Then for every $C>0$ there is $M(C)<\infty$
> such that every arc of length at most $C\sqrt R$ on $x^2+y^2=R^2$ ($R^2\in\mathbb Z$, any position)
> contains at most $M(C)$ lattice points.

The referee checked every step and constant and found no hidden uniformity assumption, no discarded gcd
or denominator cost, and no drifting phase. Theorem D gives a weaker conclusion than Theorem 4.1 (no
absolute $M_0$, no $R^\alpha$ arcs). Its value is structural: the **form-blind cross-ratio data alone**
already suffice under Vojta. Neither instance is known to imply the other.

**Why $k=8$.** For balanced configurations the normalised height vector of $\Phi$ is $w\beta+o(w)$, where
$\beta$ is the balanced class ($\beta\cdot D_S=1$ for every boundary divisor). Then
$$K\cdot\beta=2^{k-3}(k-8)+k+2,$$
which is $-7<0$ at $k=7$ and $10>0$ at $k=8$. So Vojta's $D=0$ inequality is vacuous on balanced points
for $k\le7$, and $k=8$ is the first case where it predicts a height gap.

**The arithmetic content: $(L_8)$.** Under $V_8$ the proof needs only the following lemma, which $V_8$
implies; conversely $(L_8)$ alone gives the uniform theorem (ptolemy Steps 0–5 with Step 4 replaced by
$(L_8)$). It is a statement about rank-2 integer matrices and divisibility of minors, with no Gaussian
integers and no circle. It is stated precisely in §7.1.

**What is known about $(L_k)$** (the same statement on $\overline M_{0,k}$).
* $(L_5)$ is **false** (infinitely many balanced bounded-residue genus-0 curves; `growth/referee_ptolemy.md`
  §3.3; also `threshold/L8.md` Prop C).
* $(L_6)$ is **false** (`threshold/L8.md` Thm B; refereed): a rational 3-parameter family of genus-one
  curves in $\overline M_{0,6}$ with generic positive rank sweeps out $\overline M_{0,6}$ and gives, off every
  proper closed set, balanced anchored configurations with $\eta\to0$, $w\to\infty$, $Y_j=\pm1$ and bounded
  residues. (The original `growth/ptolemy.md` claim of this was over-scoped; see §8.2.)
* $(L_7)$ is **open**; heuristics predict it is false (§6.5).
* $(L_8)$ is **open**.
* **No auxiliary-polynomial (Liouville) argument proves $(L_k)$ for any $k$**, nor any $H_\Sigma\ge f(w)$ with
  $f\to\infty$ (`threshold/L8.md` Thm A, in the precise form of `threshold/referee_L8.md` §2.3): the balanced class
  lies in the interior of the movable cone of $\overline M_{0,k}$ (Keel–McKernan).

## 5. The barrier map: fake circles and capped input classes

### 5.1 Fake circles: definition

Normalise a cluster as in §2.2. Then $z_x=\varepsilon_x\prod_{j=1}^r\pi_j^{a_{xj}}\bar\pi_j^{e_j-a_{xj}}$ with
$N=\prod p_j^{e_j}$ and $W=\sum e_j\log p_j$. Fix real lifts $\phi_j$ of $\arg\pi_j$. For $v\in\mathbb Z^r$ put
$A_v=\prod_j\pi_j^{v_j^+}\bar\pi_j^{v_j^-}$, so $\mathrm{Norm}A_v=e^{w(v)}$ and $\arg A_v=\langle v,\phi\rangle$.

> **Lemma 5.1 (monomial gaps; `round4/outside.md` Lemma 0.1 = `growth/walsh-extraction.md` Lemma 1.1).** For
> every $v\ne0$,
> $\mathrm{dist}(\langle v,\phi\rangle,\frac\pi2\mathbb Z)\ge\arcsin(e^{-w(v)/2})$ and
> $\mathrm{dist}(\langle v,\phi\rangle,\frac\pi4+\frac\pi2\mathbb Z)\ge\arcsin(e^{-w(v)/2}/\sqrt2)$.

(Proof: $A_v=X+iY$ is a conjugate-primitive non-unit of odd norm, so $|X|,|Y|,|X\pm Y|\ge1$.)

The cluster lies in an arc $\theta_0+[0,\Delta]$, so there are $\delta_x\in[0,\Delta]$ and integers $k_x$
(absorbing units and windings) with
$$\langle2a_x-e,\phi\rangle=\theta_0+\delta_x+\tfrac\pi2k_x\qquad(x=1,\dots,M).\tag{5.1}$$

> **Definition 5.2 (angle model, "fake circle"; `round4/outside.md` §1).** Fix distinct primes
> $p_j\equiv1\pmod4$, exponents $e$ and rows $a_x\in\prod[0,e_j]$. A *profile* is
> $(p,e,(a_x),(\delta_x),(k_x))$. A vector $\phi\in\mathbb R^r$ *realises* it if (5.1) holds for some
> $\theta_0$, and is *gap-admissible* if Lemma 5.1 holds for every $v\ne0$. An **angle model** is a
> profile with a gap-admissible realising $\phi$. Its "points" are $Ri^u\exp(i\langle2a-e,\phi\rangle)$.

So a fake circle has **genuine primes, a genuine exponent profile and real angles**: the angles
$\phi_j$ need not be arguments of Gaussian integers. Every actual normalised cluster is an angle model
(with $\phi_j=\arg\pi_j$). An input class *contains a fake* when some angle model of size
$W=O(M\log M)$ satisfies all its inequalities. Then, by the growth dictionary of §1.2, no argument using
only that class can prove $M=o(\log R/\log\log R)$. Residue data are added as extra data on the profile
(point residues, or prime residues $r_j(q)$ with $\mathrm{Norm}\,r_j(q)\equiv p_j$), with each character's
**actual** collision modulus (`round5/README.md` framing correction).

**Two sharpness facts.** (i) For a pinned cluster, the angle model class with Lemma 5.1 imposed **at
every prime with the full lattice-gap strength** is the actual problem: $\mathrm{LCG}_n(\phi)$ holds iff
$\phi=\arg\beta$ with $\beta$ primitive of norm $\le n$ (`round10/analytic.md` Thm D). So a
fake-violated input is easy to find if one allows the problem itself; the point is to find one that
is weaker than the problem. (ii) The witnesses of the fake-capped rows of §5.4 (rows 1–4, 7, 10–12)
have $W=O(M\log M)$, so those caps are **at the rate $\log R/\log\log R$**; only the constant varies.
Not every known fake is that small: the `round4/outside.md` A.3 models have
$\log R\le(15+o(1))M^2\log M$, the fakes of `round12/fake.md` Thm 2.1 have $W/(n\log n)$ unbounded
(Theorem 5.3), the P3-all members of `round5/round6_sharp/sharp.md` Thm 6.4 reach only
$O(M\log^2M/\log\log M)$, and the `round7/p3G.md` Prop 6.1 members have $W\le(66+o(1))MX$
(`round7/README.md`). A fake at $W=cM\log M$ shows
that the class cannot prove $M\le(2/c-\eta)\log R/\log\log R$. Actual clusters satisfy
$W\ge3M\log M-O(M)$ (§2).

### 5.2 Constructions

1. **The skew double core** (`round5/exact.md` Thm 4.1, part (iii) re-derived by `round8/referee_lw.md`;
   `round8/lw.md` Thm 6.1 parts 1–2). Let $K$ be a skew-Hadamard matrix of order $M\ge8$ ($K+K^T=2I$,
   $KK^T=MI$), normalised with $K(x,0)=1$. Put $S_a=\{x:K(x,a)=-1\}$, $S'_a=S_a\triangle\{0,a\}$
   ($a=1,\dots,M-1$). The $2(M-1)$ cuts are distinct and balanced, every pair of rows is separated by
   exactly $M$ of them, and every nonzero integral zero-sum character $c$ has support measure
   $\sigma_c\ge\|c\|_1/2$, with equality at pairs. For the Paley version ($M=q+1$, $q\equiv3\bmod4$
   prime) with all primes in one short window $[P,Pe^{1/(8q)}]$, every character has
   $w(\mathrm{supp}\,v)\ge W/4+\frac12\log P-\frac1{16}$ (`round8/lw.md` Thm B), and there are angle models with
   every monomial gap, every rational character and Haar-generic angles at
   $$W\le(8+o(1))M\log M+8M\log(1/C).$$
2. **Fully flipped Paley cores with a residue design** (`round5/construct.md`; refereed). Paley-I matrices
   have a skew pairing, and the fully flipped profile with $b$ copies per label has exact, explicit
   character margins (item 392 of the notes, restated via rank-one rectangles). With a residue design from auxiliary Gaussian primes this
   gives angle models with **point-level** residue data at moduli $\le Q_0=2M$ at
   $W\le(27+o(1))M\log M+9M\log(1/C)$; constant $9(2+\max(K,1))$ at $Q_0=M^K$. `round5/lower.md`
   gives $(16+32A+o(1))M\log M$ with residue data up to $M^A$, with these scope limits from its referee
   (`round5/README.md`): `lower.md` Lemma 5.2 needs all primes of $N$ to exceed $X$; the class should
   quantify over all four units, and the witnesses meet 2-adic and twisted collisions at $X=2M$, but **the
   twisted extension at $X=M^A$, $A>1$, is open**; Lemma S is correct but not needed (item 393 of the
   notes gives $(32+o(1))M\log M$ for every capacity-two assignment, $b=16$); all statements are
   asymptotic, with thresholds around $\log M\sim85$ at $C=1$. Like `construct.md`, these fakes carry
   point-level residues only.
3. **Parity-class designs** (`round5/round6_sharp/sharp.md` Thm 5.1; refereed). Residues that **factor
   through prime residues** with the genuine norm obey a forced parity law: one bit per prime and modulus,
   $\ell_j=[(p_j/q)=-1]\bmod2$ (Prop 2.1; proved by `round5/referee_construct.md` and checked on 211,136 actual
   instances). Designs with classes inside parity classes give **P3-char** fakes at
   $W\le(b\max(A,2)+o(1))M\log M$, $b\ge33$, moduli $\le3M^A$. Imposed on **every monomial** (P3-all) the
   same designs reach only $O(M\log^2M/\log\log M)$ (Thm 6.4).
4. **Grafts** (`round10/analytic.md` Thm B). The genuine Gaussian-prime angle function with a fake's $O(M)$
   angles grafted in.

### 5.3 The fake-test theorem

The program's working principle is that a proposed growth mechanism must use an input that fakes violate.
The record makes this precise for the determinant route, where the residual hypothesis is a statement
about a box, and fakes are entitled to test it.

> **Theorem 5.3 (fake directions; `round12/fake.md` Thm 2.1, corrected constant
> $\theta'=\theta/(8\sqrt n)$).** Let the pattern data satisfy the hypothesis (Mar) of `round12/fake.md` §2, and
> let $u^*$ be a nondegenerate real point of the tangent cone $TC_H(A,m)$ of the slab box (§6.2). Then as
> the prime-size parameter $P\to\infty$ there are gap-admissible
> fakes whose direction is within $\theta e^{-W/4}$ of $u^*$, with
> $\kappa_c(\hat u)\le\kappa_c(u^*)+m\theta e^{-W/4}/\sqrt d$. **These fakes are not at barrier scale:
> $W/(n\log n)$ is unbounded along them.** At the barrier scale $W=O((a/\rho_1+1)r\log n)$ the same
> conclusion holds **under the Diophantine condition (2.2)**: $|\lambda\cdot u^*|\ge n^{-a\|\lambda\|_1}$ for
> **all** relevant $\lambda$ with $\|\lambda\|_1\le r/\rho_1+1$, not only for pairs. By
> `round16/barrier.md` Thm 2.3, (2.2) is needed only for **light** characters
> ($\sigma_\lambda<r/2+O(n/\tau)$).

The unconditional part therefore produces fakes too large to matter for growth; only the barrier-scale
part, which needs (2.2), bears on the growth rate (`round12/README.md` §2.4 item 1). The referee also notes
that the write-up's claim that this test is "neither easier nor harder than ZE-TC" is overstated: the two
directions need different hypotheses.

What this means for route 5 (`round12/README.md` §0 item 5; `round16/README.md` §0): the hypotheses
(MIN$^{\rm pair}$, or $(R_H)$; §6.2) depend only on the box data $(J(A),V_{H'},G)$. So a box-level proof
must also hold on the (Mar) fake pattern boxes, and there a nondegenerate $TC_H$ point that satisfies (2.2)
yields gap-admissible fakes at barrier scale. **Deciding the box hypothesis on the fake pattern boxes is
therefore the barrier question.** A concrete version, the Fourier-circle test, is:

> **Theorem 5.4 (circle cap test; `round16/barrier.md` Thm 3.5 with the R16 corrections, made
> unconditional by `round17/hlight.md` Thm 2.8).** Let $q\equiv3\pmod4$ be prime, $q\ge3251$, and let
> $\mathrm{sdc}_q$ be the skew double core box at $T=1$ (which is invariant under $\mathrm{AGL}(1,q)$), posed
> at the **barrier scale** $W=2q\tau$ with $\tau\ge\tau_D$ and $\tau=O(\log q)$, together with the
> thresholds of `round16/barrier.md` Thm 2.3, and with $H=H'=e^{m'W/40}$. If
> $\inf_\theta\kappa_H(\hat u_\theta)\le\varepsilon_0/(2\sqrt{m'+1}\,e^{\pi^2/(4q)})$ on the Fourier circle
> $u_\theta(x)=\cos(\theta+2\pi x/q)$, $u_\infty=0$, at one such $q$, then route 5 is capped.
>
> After round 17, $\tau_D$ is computed with $c_0=1/6$ (all $q\ge7$) or $c_0=0.2349$ (regime $q$)
> (`round17/README.md` §1).

The hypothesis H$_{\rm light}$ that barrier Thm 3.5 needed is **proved** (`round17/hlight.md` Thm 2.8; two
independent referee passes): every integral non-pair zero-sum $\lambda$ on $\mathbb F_q\cup\{\infty\}$ has
$\sigma_\lambda\ge q/3=r/6$ for $q\ge7$, and $\ge0.2349r$ at every regime $q$ (Weil's bound for 4-term
characters); the constant is sharp up to $O(q^{-1/2})$. Whether the test's hypothesis holds is
**undecided** (§6.2).

A simpler instance of the same principle (`round10/README.md` §2 item 6; proved): give any fake its genuine
labels $\log(\pi_j/\bar\pi_j)$. They are parity-respecting, good for every norm at every level, and coherent
with $\zeta=\pi$. So, for every fake **whose primes are in the regime**, the hypotheses **and** the
conclusions of L4$(K)$, L4$^{\rm forced}$ and L4$^*$ hold with genuine labels. (L4$(K)$ is in any case false,
`round10/L4K.md` Thm 2.) Hence L4$^*$ uses no fake-violated input, and a lemma about labels alone cannot
separate actual clusters from such fakes.

### 5.4 The table of capped input classes

"Capped" means that a barrier result shows that no argument using only these inputs proves
$M=o(\log R/\log\log R)$. The barriers are of different kinds and have different status, so each row
carries a **kind** and a **status**:

* **F** (fake cap): the class contains fakes at $W\le c\,M\log M$ for a fixed $c$;
* **V** (variety theorem): a theorem on the relevant variety shows the method cannot reach the needed
  exponent; no fake is involved;
* **$\Theta$** (weaker cap): the inputs cannot even reach $o(\log R)$;
* **Red** (reduction): the inputs reduce to classes capped in other rows;
* **Lim** (limitation of known tools): the known results have the wrong dependence on the prime set; this
  is an observation, not a theorem that the inputs are capped.

| # | input class | witness | $W$ of witness | why it is capped | kind; status | source |
|---|---|---|---|---|---|---|
| 1 | **Local inputs, class $A_0$**: exponent profile, pair chords, cut/Plotkin, all monomial (Liouville) gaps, rational characters (equivalent to the monomial gaps for large primes, `round4/outside.md` Thm A.1), Haar-generic angles | skew double core | $(8+o(1))M\log M+8M\log(1/C)$ | these are statements about the profile and $\phi$ only, and Lemma 5.1 holds for every gap-admissible $\phi$ | F; proved, refereed | `round8/lw.md` Thm 6.1; `round7/README.md` §2.1 row 1 |
| 2 | **$A_X$**: $A_0$ plus point-level residues (small-prime collisions, gcds, cofactors) at moduli $\le X$, each character with its actual modulus; all round-7 inequalities valid on $A_X$ (Thm A, B.2, separation) | fully flipped Paley cores | $(27+o(1))M\log M$ at $X=2M$; $9(2+\max(K,1))$ at $X=M^K$; $(16+32A)$ at $M^A$ (scope limits in §5.2 item 2: all primes $>X$; twisted extension at $M^A$, $A>1$, open) | residue designs from auxiliary Gaussian primes reproduce every pigeonhole consequence | F; proved, refereed | `round5/construct.md`, `round5/lower.md`; `round5/README.md` |
| 3 | $A_X$ plus near-extremal **structure**: sub-cluster gcds, smoothness, prime counts, Plotkin balance, collision and level stability, rank | `round5/construct.md` models (transfer) | $(27+o(1))M\log M$ | fakes have full rank $M-1$ and smooth norms | F; proved, refereed | `round8/stability.md` Prop 7.1; `round8/weakest.md` Thm 1.3 |
| 4 | **Characters, prime-level parity**: P3-char (residues factor through prime residues with the genuine norm, imposed on characters) | parity-class designs | $(b\max(A,2)+o(1))M\log M$, $b\ge33$ | designs keep pair demands $O(\log M)$ inside parity classes | F; proved, refereed | `round5/round6_sharp/sharp.md` Thm 5.1 |
| 5 | **GAP/Fourier, energy, incidences, equidistribution** | angle models (Thm A.3); skew double core | $\log R\le(15+o(1))M^2\log M$ (A.3); then $(8+o(1))M\log M$ via row 1 | these reduce to monomial gaps and characters (A.1) | F (via row 1); proved, refereed | `round4/outside.md` Thms A.1, A.3 |
| 6 | **$L^p$, $\Lambda(p)$, energy, decoupling, local smoothing, discrete restriction** | the transcendental set $\{Re^{ik/R}\}$ | — | single-circle data are one Laplace eigenspace, so the estimates are identities; an $M$-point cluster has $\lvert F\rvert\ge M/2$ only on area $\sim C^{-3}R^{-1/2}$, so $\lVert F\rVert_p\le K_p\lVert F\rVert_2$ gives only $M\lesssim R^{1/p}$ | $\Theta$ (capped at $\Theta(\log R)$, not a fake cap); judge's assessment after a steelman/skeptic panel, not a formal referee report | `local_smoothing/README.md` §2, `local_smoothing/chains.md` Prop 2 |
| 7 | **Analytic ensemble inputs** ($(c,\gamma)$-robust, prime-supported: Hecke large sieve over primes with constant $\ge c(X+KX^{1/2})$; Hecke prime sums with error $\ge cX^\gamma$, including zero-free regions, zero density, GRH; Kubilius/Ricci/Harman–Lewis sector counts; sieve bounds in sectors and short intervals; congruence versions), with $c,\gamma,A$ fixed before the target constant | grafts of row-1 fakes ($A>\max(4,1/\gamma)$) or of P3-char fakes ($A>\max(2,1/\gamma)$) | $(2A+o(1))M\log M$, resp. $(bA+o(1))M\log M$ | the input cannot see $O(M)$ grafted angles; "almost all sectors" exceptional sets contain any $O(M)$-point selection (Prop C) | F; proved, refereed | `round10/analytic.md` Thm B, Prop C |
| 8 | **Single-place auxiliary polynomials** (Liouville at $\infty$ on any $X_M$), Ru–Vojta with $S=\{\infty\}$ | — | — | $\tau^*(M)=2-1/\lceil M/2\rceil<2$: no form on $X_M$ vanishes along the diagonal to order $\ge2\deg F$; Liouville reaches exactly the Cilleruelo–Córdoba exponents $\alpha<(k-1)/(2k-1)$, $k=\lceil M/2\rceil$, never $\alpha=1/2$; $\beta_M\le\tau^*(M)<2$ | V; proved, refereed | `threshold/upper.md`, `threshold/referee_upper.md` |
| 9 | **Auxiliary polynomials on $\overline M_{0,k}$** (any $k$) for $(L_k)$ | — | — | the balanced class lies in the interior of the movable cone | V; refereed as a formal barrier in the precise form of `threshold/referee_L8.md` §2.3 | `threshold/L8.md` Thm A |
| 10 | **Baker/Matveev** (proved linear forms in logarithms) | `round5/construct.md` (supports $\ge2\sqrt{(b-1)M}-1$); double core (supports $\ge M/2$) | rows 1–2 | they exclude only characters supported on $\le(1-\eta)\log M/\log\log M$ primes | F; proved, refereed | `round8/stability.md` Prop 6.1; `round4/flipped.md` |
| 11 | **Single-form Lang–Waldschmidt** (product or max form, any $\varepsilon\ge0$, any representation over $\mathbb Q(i)^*$) | skew double core | $(8+o(1))M\log M$ | a single form sees a character only through $w(\mathrm{supp}\,v)$, and the double core has $w(\mathrm{supp}\,v)\ge W/4+\frac12\log P-\frac1{16}$ for every character | F; proved, refereed | `round8/lw.md` Thms A, B; `round8/weakest.md` Thm 3.3 |
| 12 | $\mathrm{LW}_{\max}(\varepsilon,C_0)$ at a **fixed** $\varepsilon$, $0<\varepsilon\le1$, $0<C_0\le1$, for all forms whose list contains $\mathrm{Log}\,i$, simultaneously; Paley orders $M=q+1\ge M_2(\varepsilon,C_0)$ | skew double core | proved: $(\max(8,4/\varepsilon)+o(1))M\log M+8M\log(1/C_0)$; the referee's refinement lowers $4/\varepsilon$ to $2/\varepsilon+o(1)$ | union bound over representations | F; proved, refereed (max form only); the constant degrades as $\varepsilon\to0$ | `round8/lw.md` Thm C; `round8/README.md` §1.1, §2.2 row 7 |
| 13 | **Subspace route QS** (round 8): $QS(m,\delta,f,H_0)$ with any threshold $H_0(s)$ depending on the prime set only through $s=\lvert P\rvert$ | smooth $S$-integers in small squares (Prop A); $s$ Gaussian primes in an octant square (Prop B) | — | Prop A forces $f(s)\ge\exp(s^{1/2-o(1)})$; Prop B forces $\ge s-m+1$ subspaces, so the route gives at best $M\le\omega(N')$, which can be $(2+o(1))\log R/\log\log R$ | **unrefereed**: Props A, B are proved in `lit18/check_tools.md` §2.2, recounted by a checker (A) and re-derived by the lit18 judge (B), but not refereed | `lit18/README.md` §2.2 ("unrefereed; not leads") |
| 14 | S-unit counts, Subspace with known $S$-dependence, gcd theorems for fixed $S$ | — | — | counts exponential in $s$; exceptional sets depend on $S$ | Lim; remarks in refereed write-ups and in the literature sweep | `round8/weakest.md` §5, `round8/stability.md` Prop 8.2; `lit18/README.md` §1.2 |
| 15 | abc, Bombieri–Lang | — | — | abc is vacuous on squarefree block relations (all fakes are squarefree); for Bombieri–Lang the auxiliary varieties are rational, and the record found nothing usable ("nothing found", not a proof of no implication; §4.1) | V/vacuity; abc vacuity proved, refereed | `round8/weakest.md` Props 4.1–4.2, 6.2; `growth/referee_conditional.md` minor item 4 |
| 16 | **Interpolation determinants and point-counting on one circle**: Bombieri–Pila, Heath-Brown, Salberger determinants; gap principles, larger-sieve variants, Stepanov's method; also sums, Sidon sets, the $\sigma$-cocycle, the $\mathbb F_q$ polynomial method | — | — | on one circle the interpolation determinants factor as a Vandermonde times a symmetric factor with no extra arithmetic, so these arguments reduce to pair data and give nothing beyond $O(\log R/\log\log R)$ | Red; `growth/auxiliary.md` is **unrefereed** (no referee report); the additive-tools reduction is refereed | `growth/auxiliary.md`; `growth/README.md`; `round10/additive.md` §5 and §6 table |

**Refuted hypotheses** (not caps, but closed): $\mathrm{SLW}(k,g,\kappa)$ with union heights is false
for every $k$ at every $g<1$, and $\mathrm{SLW}^{\log}$ at every $g<1/2$ (`round8/referee_weakest.md`
Props R1–R2, polynomial templates). $QS(m\le3,\delta,o(s))$ is false (`round8/weakest.md` Prop 5.3).

**Caveat on residue caps** (round 9 R1; an unrefereed README remark). Rows 2–4 and 7 hold for residue
moduli up to $M^A$ with $A$ **fixed in advance**. Every construction of rounds 5–7 has $p_j>X$, so the
constant grows at least like $8A$. A growth argument may choose $A$ after the cluster's size parameter
$c_0$, and **caps uniform in $A$ are not proved**. Near-extremal actual clusters have all primes
$\le M^{O(1)}$ (Theorem 3.11). Fakes with all primes $\le M^{O(1)}$ do exist (for `round5/construct.md`,
$\log p=(3+o(1))\log M$; `round8/README.md` §1.3 item 2), but **no fake is known in regime (1.4)**, where
$p_j\le M^{h_0}$, $X=3M^A$ and $A$ is chosen after $c_0$ (`round9/README.md` §0 item 1;
`round14/README.md` §0; §6.3).

### 5.5 Not capped (open or unclassified)

* **P3-all$^*$**: prime-level residues with the genuine norm and parity, imposed on **every** monomial
  (§6.3). This includes the **arithmetic fibre inputs**, which the referee showed `round10/analytic.md`
  had misclassified; they are P3-all inputs (`round10/README.md` §7.2).
* **The sharp $(X+K)$ large sieve**: up to constants equivalent to the max-scale cross-gap
  $|\sin(\phi_p-\arg\varpi)|\ge c/\max(p,\mathrm{Norm}\,\varpi)$, not to the capped pair-gap form; its
  status is open and unclassified (`round10/README.md` §3 item 6, §7.2).
* **Composite-supported analytic inputs** (all-ideal large sieve, Type I/II sums): **unclassified**; Thm B
  does not cover them (`round10/README.md` §7.2; `lit18/README.md` §2.2 item 3).
* **Additive (archimedean) data** such as $\mathrm{Log}((z_x+z_y)/(2z_u))=O(\Delta)$: not modelled by angle
  models. Used in full they are the problem itself: the additive gap is violated by almost every fake,
  but its universal form is strictly stronger than growth (`round10/additive.md`).
* **ZE and its algebraic-point form** (the zero estimate of route 5; §6.2 items 1–2): open, unclassified.
* **Prime-weight concentration**: every known fake puts its primes in one window $[P,Pe^{1/(4r)}]$.
  Whether near-extremal actual clusters must do so is open (`round8/README.md` §2.5 item 6;
  `round8/stability.md` §9.3). If they had to, the problem would become one about angles of Gaussian
  primes in very short intervals.
* **The low-height window** $P_+(N)<M^{2-\delta}$: no fakes are known there; every known fake needs
  $P\gtrsim r^2\log r$ (`round10/README.md` §3).
* **Full Lang–Waldschmidt with $\varepsilon\to0$** used simultaneously (`round8/README.md` §2.5).
* **QS$^{\rm prod}$** (§6.7) and **H$^{\rm gen}$** (§6.1).

## 6. Routes tried, and where each stops

Status words: **capped** (a barrier result of §5.4 shows the route's inputs cannot prove
$M=o(\log R/\log\log R)$; the status of that barrier, refereed or not, is stated with it); **stopped** (the residual is precisely stated, no method exists in either direction,
and the judges ended the route); **open** (neither proved, refuted nor capped). The table first, then the
details.

| route | what it would give | missing statement | status | judge's estimate (round 17 §7.2) |
|---|---|---|---|---|
| Route 1 (pairs with uniform margins) | growth | L1 (uniform margins) and L2 ($H_1$ for actual clusters) | open; pair-cut certificates alone capped | — |
| Route 3 (H$^{\rm gen}$ + L1) | $W\gtrsim M\log^2M/\log\log M$ | H$^{\rm gen}$ (Poisson tails for genuine pair integers) and L1 | **open**, no method; the label-only part (L4) is refuted / fake-satisfied | $<0.1\%$ |
| Route 5 (determinant method at $T=1$) | $M=O(\log\log R)$ on all norms | MIN$^{\rm pair}$ (or its $L^2$ form $(R_H)$) | **stopped**: not capped, not unblocked; CIRCLE$(q)$ undecided | $\sim0$ for growth |
| P3-all$^*$ (local class) | growth if a lower bound (T) is proved; a cap if a fake is built | (T) or a construction C1 | **open**, uncapped, no fake known in the regime that matters | — |
| $(L_7)$ | uniform bound | $(L_7)$ | **open**; predicted **false** | proof $\sim0$; refutation plausible |
| $(L_8)$ | uniform bound | $(L_8)$ | **open**; predicted true by Vojta | $<1\%$ |
| QS$^{\rm prod}$ | growth or uniform | QS$^{\rm prod}$ | **open**, non-standard, at least as strong as the target | $\sim0$ |
| QS (round 8) | growth | $QS(m,\delta,o(s),H_0)$ | **capped** by lit18 Props A, B (**unrefereed**: proved in `lit18/check_tools.md` §2.2 and re-derived by the lit18 judge); $QS(m\le3)$ refuted (refereed) | — |

### 6.1 Route 3: genuine-prime statistics (H$^{\rm gen}$ + L1)

**The mechanism** (`round7/p3lower.md` Thm C; refereed). Fix $0<\eta<1$. For a member with $X\ge M$ let
$Q_\eta=\{q\text{ prime}: M^{1-\eta}<q\le M,\ q\text{ odd},\ q\nmid N\}$, $\lambda_Q=\sum_{q\in Q_\eta}1/m_q$ with
$m_q=q-\chi_4(q)$ the order of the norm-one subgroup of $(\mathbb Z[i]/q)^*$ (`round5/round6_sharp/sharp.md` Prop 2.1),
$c_{xy}=\#\{q\in Q_\eta:\text{unit-1 collision of }(x,y)\bmod q\}$, $F_k=\sum_{x<y}\binom{c_{xy}}k$ and
$E_k=\binom M2e_k\big((1/m_q)_{q\in Q_\eta}\big)$.

* **Hypothesis $H(\eta,\delta)$** (Poisson-type collision statistics, $0<\delta\le1$): for every $k\ge1$
  with $\delta^kE_k\ge1$, one has $F_k\ge\delta^kE_k$. For actual clusters it is a statement about joint
  divisibility of the reduced chords $(z_x-z_y)/g_{xy}$ of §2 by primes in $(M^{1-\eta},M]$.
* **L1** ($K$-uniform margins): $\max_{x,y}\mu_{xy}\le K\bar\mu$ (with $\bar\mu$ the mean margin), and
  $\lambda_Q\ge\lambda_0$.

> **Theorem 6.1 (`round7/p3lower.md` Thm C; conditional).** Fix $\eta\in(0,1)$, $\delta\in(0,1]$, $\lambda_0>0$,
> $K\ge1$, $C>0$, $\varepsilon\in(0,1)$. For $M\ge M_0$, every member of $A_X(M,C)$ with $X\ge M$ (in
> particular every actual cluster) satisfying $\lambda_Q\ge\lambda_0$, $H(\eta,\delta)$ and $K$-uniform
> margins has $W\ge((8(1-\eta)-\varepsilon)/K)\,M\log^2M/\log\log M$, hence
> $M\le(K/(4(1-\eta))+o(1))\log R\,\log\log\log R/(\log\log R)^2$.

The referee notes that the proof is Prop B.1 applied to one assumed heavy pair (hypothesis $H_1$: one
pair colliding at $\ge(2-\varepsilon/4)\log M/\log\log M$ window primes): **the growth sits in the
hypothesis**. The round-6 design $D_A$ violates $H$ (Lemma C.2), while independent uniform residues satisfy
it with probability $\to1$ (Lemma C.3; that model ignores the forced parity). **H$^{\rm gen}$** is $H$
stated for the genuine pair integers $c^{\rm gen}_{xy}=\#\{q\in Q_\eta:q\mid D_{k_x-k_y}(Z^\varepsilon_{a_x-a_y})\}$
(`round14/rigtheory.md` §4.2). The round-14 transfer theorem (rigidity + H$^{\rm gen}$ + L1 $\Rightarrow$
$W\gtrsim M\log^2M/\log\log M$ for P3-all$^*$ members), valid after the associate-twist correction, routes the
class-level question to H$^{\rm gen}$. For actual circles alone the theorem is vacuous: it only routes
P3-all$^*$ to statements about genuine labels (`round14/README.md` §1.2 item 2). Its size conclusion needs
the extra hypothesis $r\ge M^\gamma$, which regime (1.4) does not supply. Round 15 (`round15/README.md`;
`round15/rigSLtheory.md`, **unrefereed**: its referee was cut off) argues that no outcome of the rigidity question
moves the growth rate by itself.

**What was tried and where it stops.**
* *The label-only inverse step* L4$(K)$ of round 9 is **false** for every $K\ge2$ (`round10/L4K.md` Thm 2):
  bounded norms $\le K$ force collisions only up to $M^{\lfloor K/2\rfloor+o(1)}$ (Thm 1). The surviving
  L4$^*(K)$ is satisfied, hypotheses and conclusion, by every fake given its genuine labels (§5.3). So any
  separating power must come from the label–angle step.
* *Unconditional pieces in the regime* (`round9/L4.md` Thms A–C; refereed, given Linnik's theorem and
  Heath-Brown's quadratic large sieve): exactly coherent certificates have zero slack, and all but
  $O(M^{0.45}\log^2M)$ columns of any parity-respecting system cannot be Hamming-certified by Gaussian
  integers of norm $\le M^h$ and sub-genuine height. These concern label systems, not actual circles.
* *Literature* (`lit18` §1.3). Kubilius/Poisson models (Ford arXiv:2408.03803; Gorodetsky
  arXiv:2608.14190) have the right shape for $H$ but impossible strength: one heavy pair needs a pair
  integer divisible by $\ge M^{(1-\eta)k}$ with $k\sim2\log M/\log\log M$, about $10^{30}$ at $M=10^6$,
  against a family of $\binom M2\approx5\cdot10^{11}$ adversarially chosen integers. The combinatorial large
  sieve (Croot–Mao–Pohoata–Sheffer–Yip, arXiv:2606.17487) is the same mechanism as §2, at the same rate.
* *The structural obstacle* (`round10/analytic.md` Prop C, Thm D). Robust ensemble statements cannot see
  an adversarial $O(M)$-prime selection; pointwise statements reconstitute the problem.

**Precise missing statement.** H$^{\rm gen}(\eta,\delta)$ for genuine pair integers of every near-extremal
actual cluster, together with L1 (or an argument that avoids L1; Theorem A of `round7/p3lower.md` shows that the
pair inputs themselves allow adapted margins at real weights, so L1 needs non-pair input). **Open; no
method.**

### 6.2 Route 5: the determinant method at box size $T=1$ (rounds 10–17)

**Idea.** Choose the auxiliary form **after** the cluster's profile. Then forced divisibility is
LP-supercritical and every large actual cluster satisfies an exact algebraic system (§3.4, last
bullet). Route-5 files write $n$ for the number of points ($n=M$) and $r$ for the number of split primes. Spend the arithmetic once (rigidity, "Step 0"), and reduce growth to an inequality about the
cluster's own slab box. The chain of results:

1. **Round 10** (`round10/determinant.md`; refereed). LP supercriticality (Thm 2.1), rigidity (Thm 4.1),
   and a **conditional** theorem: under the zero estimate ZE (no nondegenerate $\delta$ near $e$ lies in the
   base locus $Z_P$), clusters on squarefree norms have $n=O(\log\log R)$ (Thm 5.7). ZE is non-standard: an
   unlikely-intersection statement for a toric variety against a linear space of codimension $\ge|A|/2$.
   All of its difficulty is in the slab cuts (Prop R1).
2. **Round 11.** ZE is open. Every stronger or $\delta$-free reformulation proposed (ZE-braid,
   SL, SL$^*$, ZE-col) is **false**, through cheap rational subtori of Hadamard-type profiles. The
   refutations of ZE-braid, SL and ZE-col are the round-11 referees'; the refutation of SL$^*$ is the
   round-11 judge's finding **J1** (unrefereed in round 11; `round11/README.md` §0 item 1, §5.1), confirmed
   in round 12 (`round12/README.md` §2.1 item 10: "J1, J2, J3 (R11) are confirmed"). ZE-col$^*$
   (restricted to the actual point) would give $M=O_C((\log\log R)^2)$ on all norms (`round11/subcluster.md` Thm 5.4,
   corrected; conditional).
3. **Round 12** (refereed). At $T=1$ the machinery is **proved for every profile and every norm**, with no
   squarefree, rank or nondegeneracy hypothesis: the slab box
   $A=\{k\in\mathbb Z^n_0:k'\in\{-1,0,1\}^d,\ |\langle k,\beta_p\rangle|\le w\ \forall p\}$, $d=n-1$,
   $w=\sqrt{8d\ln4r}$, has $|A|\ge3^d/2$ and $m_{\max}\ge4w$ once $d\ge316\ln4r$; rigidity holds for every
   order $m'\in[\lceil4w\rceil,m_{\max}]$ once $W>W_1$; Step 0 gives $\kappa(\hat u)\le\varepsilon_0=Cde^{-W/4}/(m'+1)$;
   and the **arithmetic separation lemma** $\mathrm{dist}(\hat u,L)\ge\frac12\exp(-\binom n2\bar a/(2P(L)))$
   holds for every subspace $L$ with $P(L)$ collision pairs. The cheap-span form of ZE-TC is **false**
   (Sylvester $2^{12}$ plus one light column, $m'=2256$).
4. **The residual hypothesis.** With $\Phi_c(u)=\sum_{k\in A}c_k\langle k',u\rangle^{m'}$,
   $\kappa_c(u)=|\Phi_c(u)|/(|c|_1d^{m'})$ and $\kappa_H(u)=\sup_{c\in V_H\setminus0}\kappa_c(u)$, where $V_H$ is
   spanned by the integral $c\in V(A,m')$ with $|c|_\infty\le H$:
   $$\textbf{(MIN}^{\rm pair}\textbf{)}\qquad\min\{\kappa_H(u):|u|_\infty=1,\ |u_x-u_y|\ge e^{-\epsilon_{xy}/2}\ \forall x\ne y\}>\varepsilon_0,$$
   with $\epsilon_{xy}=d_{xy}-W/2+2\ln C'\ge0$ the pair excesses ($C'$ the arc constant of `round12/t1rigour.md`,
   $C'=C$ in all numerics). **If MIN$^{\rm pair}$ holds for the cluster's own $T=1$ box at
   $m'=\lceil4w\rceil$, for all regime profiles, with $n\ge n_1$, $d\ge316\ln4r$ and all $W>W_1$, then
   $n=O(\log\log R)$ on all norms** (`round12/nearstrata.md` Thm 4.1, Cor 4.2, as strengthened by its
   referee; proved implication; `round12/README.md` §2.3 items 5–7). The threshold $n_1$ is not explicit:
   the related MIN$^*$ is false on the regime Sylvester box $n=4096$ ($m'=2256$, $g>328.4$), so
   $n_1>4096$. No bound depending only on dimension, degree, height and the qualitative
   zero set can supply the needed size (Prop 2.7).
5. **Rounds 13–16.** Locality of the leading-form ideal: neither proved nor refuted in the regime, and a
   proof would not settle the barrier (a fake model satisfies it). Round 16 recasts MIN$^{\rm pair}$ as an
   $L^2$ statement:
   $$\textbf{(}R_H\textbf{)}\qquad\min_{u\in G_1}\Psi_{A,H'}(u)>d^{2m'}\varepsilon_0^2,\qquad H'=e^{m'W/40},$$
   where $\Psi_{A,H}(u)$ is the squared $L^2(\mu_A)$ norm of the top orthogonal component of
   $\langle k',u\rangle^{m'}$ (with $V_H$ in place of $V$) and $G_1$ is the gap sphere. $(R_H)$ implies
   MIN$^{\rm pair}$ losslessly; MIN$^{\rm pair}$ implies $(R_H)$ at level $\varepsilon_0 3^{-d/2}$. Every
   anti-concentration bound depending only on degree, variance, influence or rank is provably vacuous at
   degree $\ge2j+1$ (Prop 7.1). $(R_H)$ depends only on $(J(A),V_{H'},G)$, so a box-level proof must also hold on the
   (Mar) fake pattern boxes, where by Theorem 5.3 a nondegenerate $TC_H$ point satisfying (2.2) yields
   gap-admissible fakes at barrier scale (`round16/README.md` §0).
6. **Rounds 16–17: the CIRCLE$(q)$ test.** For $q\equiv3\bmod4$ prime, $q\ge3251$, $\mathrm{sdc}_q$ at $T=1$,
   at the barrier scale $W=2q\tau$, $\tau\ge\tau_D$, $\tau=O(\log q)$, with the thresholds of barrier Thm 2.3:
   decide whether $\inf_\theta\kappa_H(\hat u_\theta)\le T^*:=\varepsilon_0/(2\sqrt{m'+1}\,e^{\pi^2/(4q)})$ with
   $H=H'$. "Yes" caps route 5 (Theorem 5.4). "No" would be the first non-zonotopal lower-bound certificate at
   an irrational direction of an overdetermined regime box. Round 17 proved H$_{\rm light}$ (so the cap test
   is unconditional) and two structural results (`round17/survive.md` Thms 3.1, 3.4; refereed): on the
   Fourier circles $\kappa_H$ is either identically $0$ or at least $2^{-m}/(|A|Hq^m)$, and any integral
   witness gives a Liouville-free, $\theta$-uniform lower bound. So the test reduces to two exact questions:
   is $C_1\subset TC_{H'}$ (then the test caps route 5)? Does an integral $c\in V(A_q(j),m')$ with $\Phi_c\not\equiv0$
   on $C_1$ and $\ln|c|_1<110q$ exist at $q=3251$ (then the test passes)? Product-type and sub-cube witnesses
   are provably excluded (Props 4.1, 4.2).

**Status: stopped** (`round17/README.md` §3). CIRCLE$(q)$ is undecided at every regime $q$; no cap and no
lower bound at regime size. Every certificate type the record can analyse stops at order
$\le8jq/(q+1)\approx m'/2$; counting (Bombieri–Vaaler) reaches the witness only at height
$\exp(e^{1073})$; toy computation cannot reach $m'/(2L+1)=4$. Two points for anyone who reopens it: a "yes"
would cap route 5 in its MIN$^{\rm pair}$/$(R_H)$ form and leave only the actual-point variant, which needs a
measure of algebraic independence for about $n$ logarithms (none exists in the literature, `lit18` §1.2); a
"no" would not prove route 5, since other barrier points are not excluded. **Route 5 cannot give growth in
either outcome without a new input.**

### 6.3 P3-all$^*$: genuine prime-level residues on every monomial (rounds 6–15)

**The class** (`round5/round6_sharp/sharp.md` Def. 1.1; class P3-all$^*$ as fixed in `round7/README.md`).
A member is an angle model (Def. 5.2) together with, for every odd prime $q\le X$ and every $j$ with
$p_j\ne q$, a unit $r_j(q)\in\mathbb Z[i]/q^{A_q}$ with $\mathrm{Norm}\,r_j(q)\equiv p_j$, such that every
$v\in\mathbb Z^r\setminus\{0\}$ and every $s\in\mathbb Z/4$ satisfy the **slab condition**
$$|\sin(\langle v,\phi\rangle-s\pi/4)|\ \ge\ \nu_s\,m_{(v,s)}\,e^{-w(v)/2},\qquad m_{(v,s)}=\prod_{q\le X,\,q\nmid N}q^{a_{q,s}(v)},$$
with $a_{q,s}(v)=\max\{a\le A_q:\prod_j(r_j/\bar r_j)^{v_j}\equiv i^s\bmod q^a\}$, $\nu_s=1$ ($s$ even),
$2^{-1/2}$ ($s$ odd), plus the character conditions (e-char). Normalised actual clusters are members for
every $X$, with their genuine residues $r_j(q)=\pi_j\bmod q^{A_q}$ (`round5/round6_sharp/sharp.md` Prop 1.2). The residues
then obey the forced Legendre parity (§5.2 item 3). The classes are nested:
P3-G $\subset$ P3-all $\subset$ P3-char $\subset A_X$.

**Bracket** (proved): $3M\log M-O(M)\le W^{\rm all}_{\min}\le(66A+o(1))M\log^2M/\log\log M$ at $X=3M^A$
(`round5/round6_sharp/sharp.md` Thm 6.4). With all $p_j>X\ge M$ the lower end rises to $(4-o(1))M\log M$ (`round7/p3lower.md` B.3).

**Target (T)** (`round7/README.md` §3). For some fixed $A\ge1$, every member of $F^{\rm all}_{3M^A}(M,C)$
has $W\ge f(M)M\log M$ with $f\to\infty$. This gives $M\le(2/f(M)+o(1))\log R/\log\log R$ for actual
circles. Round 9 sharpened the quantifiers (proved reduction, via Theorem 3.11): it suffices to prove

> **(T$'$)** For every $c_0$ there is $A(c_0)$ such that, for large $M$, no member of $F^{\rm all}_{3M^A}(M,C)$
> with all $p_j\le M^{2c_0-3+o(1)}$ has $W\le c_0M\log M$.

The setting of (T$'$), $p_j\le M^{h_0}$ with $h_0=2c_0-3+o(1)$ and $X=3M^A$ with $A$ chosen after $c_0$, is
called **regime (1.4)** in the record.

**What is proved about it.**
* *Pair-cut certificates are capped* (`round7/p3lower.md` Thm A, Cor A.5): with real weights the pair channel is
  a cut LP, $2\lambda_{\max}(D)\le W_{\rm cut}(D)\le\pi\lambda_{\max}(D)$, and the Thm 6.4 members have
  $\lambda_{\max}(D)\le36M\log M+111M$. A proof of (T) must use the monomial conditions or the discreteness
  of prime weights.
* *Parity, heights and PSD give constant factors only* (`round7/p3lower.md` B.3, B.5).
* *Several construction methods provably cannot cap it* (`round7/p3Dconstruct.md` Thms 1–2, `round7/p3Dlower.md`):
  independent random labels on primes of reciprocal mass $(\log M)^{-C}$, additive labels, exactly
  Gaussian-coherent labels (certificate level, under GRH; unconditional in regime (1.4) by round 9), and
  G-coherent designs on the Paley profile unless their column and pair heights are $\ge M^{1/2-o(1)}$.
* *Rigidity of label systems* (round 14): Rig$(K)$ is undecided; certified lower bounds on the coherent
  joint core for $\{5,13,17,29\}$ are $61,109,193,239$ for $K=4..7$; a proved ceiling
  $q_{\rm core}(K)=O_K(r^{\lfloor K/2\rfloor})$ means bounded $K$ reaches scale $M$ only when $r\ge M^\gamma$.
  Round 15 (interrupted; **unrefereed**): no outcome of the single-angle rigidity question Rig$^{\rm SL}$
  changes the growth rate of actual circles by itself.
* *Genuine residues* (P3-G): on six verified random prime sets the class collapses to the actual problem
  (exact); greedily selected prime sets evade every computed closure. Open.

**Status: open.** It is the only local class valid for actual circles that is not known to be capped.
No fake of any kind is known in regime (1.4) ($p_j\le M^{h_0}$, $X=3M^A$, $A$ chosen after $c_0$); every
construction of rounds 5–7 has $p_j>X$. A member there (construction C1) would cap the last local class;
a forced lower bound there would be the first local growth mechanism.

### 6.4 Route 1: pairs with uniform margins (`round7/README.md` §3; `round9/L1L2.md`)

* **L1 (uniformity).** A valid input forcing near-extremal actual clusters to have
  $\max\mu_{xy}\le K\bar\mu$. Open; the pair inputs allow adapted margins (Thm A).
* **L2 ($H_1$ for actual clusters).** Some pair's reduced chord (a nonzero Gaussian integer of modulus
  $\le2e^{d_{xy}/2}$) is divisible by $\ge(2-\varepsilon/4)\log M/\log\log M$ primes $q\in(M^{1-\eta},M]$.
  Open; $D_A$ violates $H$, so L2 must use genuineness.
* Round 9 (refereed): L1 and L2 are mutually exclusive on any near-extremal cluster, not each other's
  negation; Route 1 is an exact split of the growth statement into an archimedean and an arithmetic half,
  not a reduction. The spectral form of L1 is vacuous by Thm A.

### 6.5 $(L_7)$ (`round4/L7.md`, `round7/L7.md`, `threshold/L8.md` §5)

$(L_7)$ is the Ptolemy lemma on $\overline M_{0,7}$ (§7.1 with $k=7$). It would give the uniform bound,
but no standard conjecture predicts it: $K\cdot\beta=-7<0$, so Vojta is vacuous on balanced points. The
record expects it to be **false** (expected dimension of balanced rational curves $3>0$ in the integral
chart; lattice count $e^w$).
* Refutation criterion (round 4 Thm 2(a); refereed): $(L_7)$ is false if good curves with classes within
  $\varepsilon$ of the balanced ray have Zariski-dense union, for every $\varepsilon>0$.
* Round 7 (refereed): rich good $\mathbb Q$-curves are Zariski dense (Thm B: the rulings of the web quadrics
  through six cube points form an irreducible 3-dimensional $\mathbb Q$-family $V_A$, general member of class
  $\gamma_A$); but no approximately balanced good curve avoiding $\mathrm{Sing}(F_A)$ is a moving section of a
  $V_A$-ruled surface ($\lambda_A(\beta)\le-23$, Prop 5.4, exact). Whether passages through $\mathrm{Sing}(F_A)$
  can make up the deficit is open.
* The $(L_6)$ mechanism does not extend: the relevant fibre products have genus $\ge5$ (Faltings), and
  lifting by a seventh point is overdetermined by 10 (`threshold/L8.md` §5).

**Status: open; not a growth route; proof attempts not recommended.**

### 6.6 $(L_8)$

**Status: open.** It implies the uniform theorem, is predicted by Vojta ($K\cdot\beta=10>0$), and is not
tested by fakes, since it involves no circle.

*Barriers against proving it.* The round-17 judge's claim that the known unconditional methods do not
reach $(L_8)$ rests on two items of `round17/README.md` §7.1 (§7.3 there):
* **refereed**: no auxiliary-polynomial (Liouville) argument on $\overline M_{0,k}$ proves $(L_k)$
  (`threshold/L8.md` Thm A, refereed as a formal barrier in the precise form of
  `threshold/referee_L8.md` §2.3);
* **unrefereed**: no $S$-integer argument whose threshold depends on the prime set only through $|S|$
  works (lit18 Props A, B; §5.4 row 13).

Separately, abc vacuity on squarefree relations is proved and refereed in `round8/weakest.md`
Props 4.1–4.2. The author's own discussion in `growth/ptolemy.md` §7, **informal and not a list of
refereed barrier theorems**, explains why further tools stall: elimination (the Plücker variety is rational,
so no polynomial identity gives a positive multiple of $w$); single $S$-unit equations and abc (vacuous:
$n_An_Bn_C$ essentially squarefree; the gap appears only when all $\binom84=70$ Ptolemy equations are taken
jointly); and Subspace / Ru–Vojta / Evertse ($S$ unbounded, exceptional sets depend on $S$; the $D=0$ form
involves global heights $h_{D_T}$ summing contacts over unboundedly many primes, a GCD-type statement on an
iterated blow-up of $\mathbb P^5$). The round-17 judge's estimate for a
proof in a further round: **well below 1%**. Recommended first tasks (evidence only): a census of small-height
balanced anchored 8-configurations to locate the exceptional set $Z$, and a lifting test for the
$(L_6)$-refuting families $E_t$. The record contains no result on these tasks yet.

### 6.7 The Subspace route: QS (capped) and QS$^{\rm prod}$ (open)

**QS** (`round8/weakest.md` §5.1). Over $K=\mathbb Q(i)$, for a finite set $P$ of split primes, $s=|P|$,
consider the system $\Sigma_m(P)$: at $\infty$ the forms $X_1$ and $X_i-X_1$ ($i=2..m$), at each finite place
above $P$ the forms $X_1,\dots,X_m$. $QS(m,\delta,f,H_0)$ asserts: for every $P$, the solutions of
$\prod_{v\in S_P}\prod_i|L_{i,v}(x)|_v/\|x\|_v\le H(x)^{-m-\delta}$ with $H(x)\ge H_0$ lie in at most $f(s)$ proper
subspaces. **Proved implication** (Thm 5.2): $QS(m,\delta,f,H_0)$ gives $M\le f(\omega(N'))+m-1$.
$QS(m\le3)$ is false (Prop 5.3; refereed). **Capped, by unrefereed propositions**: for every $m$ and every
threshold depending on $P$ only through $s$, the route gives at best $M\le\omega(N')$ (lit18 Props A, B,
proved in `lit18/check_tools.md` §2.2, recounted by a checker and re-derived by the lit18 judge,
**unrefereed**; `lit18/README.md` §2.2).

> **QS$^{\rm prod}$** (`lit18/check_tools.md` §2.5; open, non-standard). For some $m\ge2$,
> $0<\delta<(m-1)/2$, $0<c<1/2$ and $f(s)=o(s)$: for every finite set $P$ of split primes, the solutions of
> $\Sigma_m(P)$ with $H(x)\ge\max(H_0,(\prod_{p\in P}p)^c)$ lie in at most $f(s)$ proper subspaces of $\mathbb Q(i)^m$.

It would give $M=o(\log R/\log\log R)$ (bounded $f$: uniform), it escapes Props A, B and the fakes (fakes are
not $S$-integers), but it covers unequal-norm tuples, so it says strictly more than needed, and it belongs to
the genuine-integer class. The judge's crude random-model count (unrefereed) supports it only
asymptotically: the margin is a single $\log s$ factor, and any admissible $f$ must be $\exp(\Omega(s))$ up to
$s\sim5\cdot10^3$. **Do not spend proof effort on it**; a refutation attempt is cheap.

### 6.8 Other conditional statements on record (proved implications, open or false hypotheses)

| hypothesis | conclusion | status of hypothesis | source |
|---|---|---|---|
| ZE (zero estimate for the slab-cut base locus) | $n=O(\log\log R)$ on squarefree norms | non-standard, open | `round10/determinant.md` Thm 5.7 |
| ZE-col$^*$ (ZE restricted to the actual point) | $M=O_C((\log\log R)^2)$ on all norms | open | `round11/subcluster.md` Thm 5.4 |
| MIN$^{\rm pair}$ / $(R_H)$ for the cluster's own $T=1$ box, all regime profiles, $n\ge n_1>4096$, $d\ge316\ln4r$, $W>W_1$ | $n=O(\log\log R)$ on all norms | open; route stopped | `round12`, `round16` |
| $\mathrm{SLW}_w(k,g,\kappa)$, $k\ge6$, $g<1/4$ (coefficient-weighted heights) | uniform: $M\le(5.5^{k-1}+1)/2$ | new; false for $k\le5$; open for $k\ge6$; no heuristic support | `round8/referee_weakest.md` §3.1 |
| $\mathrm{SLW}_w^{\log}(\theta,g,\kappa)$, $k_0\ge6$, $\theta<1/\log5.5$, $g<1/4$ | $M=O((\log R)^{1-\eta})$ | new, open | same |
| rank of exponent lattice $\le\theta\log W'$, $\theta<1/\log5.5$ (resp. bounded) | $M=O((\log R)^{1-\eta})$ (resp. uniform; converse holds) | open; fakes have full rank | `round8/weakest.md` Thm 1.3 |
| one instance $V(m_0,\varepsilon)$, $m_0\ge7$, $\varepsilon<(m_0-6)/2$ | uniform | standard (Vojta), open | `round8/weakest.md` Prop 6.1, threshold corrected by referee |
| Lang–Waldschmidt two-logarithm bound with naive heights over $\mathbb Q(i)$ | fully flipped Hadamard phase lemma when flipped mass $\le(1/2-\delta)W$ | a Lang–Waldschmidt-type conjecture, open | `round4/flipped.md` Thm L |

## 7. Open problems for experts

Each problem is stated in the smallest self-contained form the record supports, with what a solution
would give. Labels refer to §§4–6.

### 7.1 $(L_8)$: a statement about integer $2\times8$ matrices

*Source:* `growth/ptolemy.md` §7, restated in `threshold/L8.md` §1 and `round17/README.md` §7.3.

**Data.** An integer matrix $P=(P_0\,|\,P_1\,|\cdots|\,P_7)\in\mathbb Z^{2\times8}$ with primitive columns,
$P_0=(1,0)^T$ and $P_i=(X_i,Y_i)^T$, together with positive integers $n_T$ for the nonempty
$T\subseteq\{1,\dots,7\}$, such that for all $0\le i<j\le7$
$$\det(P_i,P_j)=t_{ij}\prod_{T\ni i,j}n_T,\qquad t_{ij}\in\mathbb Z\setminus\{0\}$$
(the empty product is $1$, so $t_{0j}=Y_j$), and such that **at each rational prime $p$, the blocks $T$
with $p\mid n_T$ form two nested chains on disjoint row sets** (the structure that item 532 produces from
threshold layers of an actual cluster). Put
$$H_\Sigma=\sum_{j=1}^7\log|Y_j|+\sum_{1\le i<j\le7}\log|t_{ij}|,$$
let $w>0$ be a scale with **balance** $|\log n_T-w|\le\eta w$, and let $\Phi=[P_0:\dots:P_7]\in M_{0,8}(\mathbb Q)$
be the point of the moduli space given by the eight columns as points of $\mathbb P^1(\mathbb Q)$ (distinct,
since $t_{ij}\ne0$).

> **$(L_8)$.** There are $c>0$, $\eta_0>0$ and a proper Zariski-closed $Z\subset\overline M_{0,8}$ such that
> every such configuration with $\Phi\notin Z$ and $\eta\le\eta_0$ satisfies $H_\Sigma\ge c\,w-O(1)$.

**What is known.** $(L_8)$ implies the uniform theorem; Vojta ($D=0$) on $\overline M_{0,8}$ implies $(L_8)$
(§4.2). The analogues $(L_5)$, $(L_6)$ are false; $(L_7)$ is open and expected false (§6.5). No
auxiliary-polynomial argument proves any $(L_k)$. The balanced class has $K\cdot\beta=10>0$ at $k=8$, the
first positive value. A proof would be an unconditional case of Vojta's conjecture with $D=0$ on a
5-dimensional rational variety, for a structured but unbounded family of points. The visible blocks are
$2\le|T|\le6$; singletons and the top block $T=\{1,\dots,7\}$ enter no Plücker relation.

**Smallest useful experiments** (evidence only, recommended by `round17/README.md` §7.3): enumerate
balanced anchored 8-configurations exactly at increasing $w$, record $\min H_\Sigma$ and the locus where it
stays bounded (the candidate $Z$; candidates are pull-backs of the $(L_6)$ genus-one families, forgetful and
Kapranov curves); test whether the $(L_6)$-refuting families $E_t\subset\overline M_{0,6}$ lift to balanced
8-point configurations with bounded residues.

### 7.2 The route-5 box inequality: MIN$^{\rm pair}$ / $(R_H)$, and its fake test

*Source:* `round12/README.md` §§1, 2.3; `round16/README.md` §§0, 5.1; `round17/README.md` §§1–3.

**Data.** $n$ points, $d=n-1$, $r$ split primes, centred profile columns $\beta_p\in[-1,1]^n$, regime
$d\ge316\ln4r$. Put $w=\sqrt{8d\ln4r}$, $j=\lfloor w/2\rfloor$, $m'=\lceil4w\rceil$, and
$$A=\{k\in\mathbb Z^n:\textstyle\sum k=0,\ k_x\in\{-1,0,1\}\ (x<n),\ |\langle k,\beta_p\rangle|\le w\ \forall p\}.$$
$V(A,m')$ is the space of $c\in\mathbb Q^A$ annihilating every polynomial of degree $<m'$ on $A$;
$V_H$ is spanned by its integral members with $|c|_\infty\le H$. Write $k'=(k_1,\dots,k_d)$ and
$\bar a=W/(2d)+2\ln C'$. For $u\in\mathbb R^d$ (the relative positions $u_x=\delta_x-\delta_n$ of the points),
$\Phi_c(u)=\sum_kc_k\langle k',u\rangle^{m'}$ and $\kappa_H(u)=\sup_{c\in V_H\setminus0}|\Phi_c(u)|/(|c|_1d^{m'})$.

> **MIN$^{\rm pair}$.** For the box $A$ of the cluster's own profile and its pair excesses
> $\epsilon_{xy}\ge0$ (with $\sum\epsilon_{xy}\le\binom n2\bar a$),
> $$\min\{\kappa_H(u):|u|_\infty=1,\ |u_x-u_y|\ge e^{-\epsilon_{xy}/2}\ \forall x\ne y\}>\varepsilon_0=\frac{C\,d\,e^{-W/4}}{m'+1},\quad H=|A|^2.$$

> **$(R_H)$ (the $L^2$ form; `round16/README.md` §0).** Let $\mu_A$ be uniform on $A$, let
> $\Psi_{A,H}(u)$ be the squared $L^2(\mu_A)$ norm of the component of $\langle k',u\rangle^{m'}$ orthogonal to
> $V_H^\perp$ (the "top component"), and let $G_1$ be the gap sphere $\{|u|_\infty=1,\ |u_x-u_y|\ge e^{-\epsilon_{xy}/2}\}$.
> Then $\min_{u\in G_1}\Psi_{A,H'}(u)>d^{2m'}\varepsilon_0^2$ with $H'=e^{m'W/40}$.

$(R_H)$ implies MIN$^{\rm pair}$ losslessly, and MIN$^{\rm pair}$ implies $(R_H)$ at the lower level
$\varepsilon_03^{-d/2}$. $\Psi_A$ is a sum-of-squares form whose zero set is the tangent cone; the slab
polynomials enter only through the leading-form ideal $J(A)$.

MIN$^{\rm pair}$, if it holds for the cluster's own $T=1$ box at $m'=\lceil4w\rceil$ for **all regime
profiles**, with $n\ge n_1$ ($n_1>4096$, not explicit), $d\ge316\ln4r$ and all $W>W_1$, gives
$n=O(\log\log R)$ on all norms (§6.2 item 4). It is open, with no evidence either way at regime
sizes, which start at $n\approx2966$ (the smallest LP-supercritical Sylvester model already has $n\ge32$).
Since it depends only on the box data, a box-level proof must also hold on the fake pattern boxes, where by
Theorem 5.3 a nondegenerate real tangent-cone point satisfying the Diophantine condition (2.2) yields
gap-admissible fakes at barrier scale. So the cleanest sub-question is the barrier test:

> **CIRCLE$(q)$ / $(W_q)$.** For $q\equiv3\pmod4$ prime, $q\ge3251$, and the skew double core box
> $A_q$ at $T=1$: is the Fourier circle $C_1$ contained in $TC_{H'}(A_q,m')$, $H'=e^{m'W/40}$? Or is there an
> integral $c\in V(A_q(j),m')$ with $\Phi_c\not\equiv0$ on $C_1$ and $\ln|c|_1<110\,q$ (at $q=3251$)?

A "yes" to the first caps route 5; a "yes" to the second passes the test (but does not prove route 5).
Nobody knows how to construct such a $c$ or rule one out.

### 7.3 H$^{\rm gen}$: Poisson tails for genuine pair integers

*Source:* `round7/p3lower.md` §§4, 7 (H-actual); `round14/rigtheory.md` §4.2; `lit18/README.md` §1.3.

> **(H-actual).** Are there fixed $\eta,\delta>0$ such that every normalised near-extremal actual cluster
> satisfies $H(\eta,\delta)$: for every $k\ge1$ with $\delta^kE_k\ge1$,
> $\sum_{x<y}\binom{c_{xy}}k\ge\delta^kE_k$, where $c_{xy}$ counts the primes $q\in(M^{1-\eta},M]$, $q\nmid N$,
> dividing the reduced chord $(z_x-z_y)/\gcd(z_x,z_y)$ of the normalised points (unit-1 collisions), and
> $E_k=\binom M2e_k((1/m_q)_q)$ with $m_q=q-\chi_4(q)$ (§6.1)?

With $K$-uniform margins (L1) it gives
$M\le(K/(4(1-\eta))+o(1))\log R\,\log\log\log R/(\log\log R)^2$ (Theorem 6.1). Two remarks constrain any
proof. The integers involved are $\binom M2$ dense Gaussian integers chosen by an adversary from an
$O(M)$-prime profile, so averaged ("almost all") statements are capped (`round10/analytic.md` Thm B,
Prop C). And one heavy pair already needs divisibility by $\ge M^{(1-\eta)k}$, $k\sim2\log M/\log\log M$.
The weakest single statement that suffices is **L2** ($H_1$): one pair with
$\ge(2-\varepsilon/4)\log M/\log\log M$ window-prime divisors of its reduced chord (still with L1).

### 7.4 P3-all$^*$: a purely local finite-combinatorial question

> **(T$'$)** (§6.3). For every $c_0$ there is $A(c_0)$ such that, for large $M$, no member of
> $F^{\rm all}_{3M^A}(M,C)$ with all $p_j\le M^{2c_0-3+o(1)}$ has $W\le c_0M\log M$.

Equivalently in the negative direction, **C1**: construct factored residue data (genuine norms, forced
parity) on a Paley-type profile with $\log P=O(\log M)$, every pair demand $O(\log M)$ and
$m_v\le M^{\beta\|v\|_1}$ for every $v$ outside the character lattice, in the regime $p_j\le M^{h_0}$. Either
answer is new: (T$'$) would be the first local growth mechanism; C1 would cap the last local class.

### 7.5 Smaller open questions recorded by the judges

* **The fully flipped Hadamard phase lemma** (`growth/walsh-extraction.md` §10): for every $C,K$ there is
  $M_0$ such that no cluster with $M\ge M_0$ points and $W\le KM\log M$ has a fully flipped Hadamard-core
  profile ($b\ge5$ copies, one flip per row at distinct primes, bounded capacity). With Theorems B'' and C
  this would settle every Hadamard-core profile with at most one flip per row. Lang–Waldschmidt with naive
  heights would give it below the critical band (§3.1).
* **The five-point question** (Corollary 3.8): finiteness of the bounded-$L$ rigid five-point system
  (`growth/cotangent.md` (B4), (B5)).
* **Bounded rank** (Theorem 3.12): is the exponent lattice of a short-arc cluster of bounded rank? This is
  equivalent to the uniform theorem.
* **$\beta_M>2$?** Closed negatively: $\beta_M\le\tau^*(M)<2$ (§5.4 row 8).
* **QS$^{\rm prod}$** (§6.7): refute it for $m=2$, $\delta,c\to1/2$ (cheap; decisive in one direction).
* **The low-height window** $P_+(N)<M^{2-\delta}$: no fakes are known; is there one?

### 7.6 What any successful argument must avoid (`round17/README.md` §7.1)

Read together, the barrier results of §5.4 say that an argument giving $o(\log R/\log\log R)$ or a
uniform bound must satisfy **all** of the following. Items 1–4 rest on refereed barrier theorems. Item 5
rests on **unrefereed** propositions, item 6 is the status of a stopped route rather than a barrier
theorem, and item 7 is a judge's assessment.

1. **No local input.** Inputs satisfied by fake circles are capped: pairs, gcds, cut, residues and
   characters (rounds 4–7, 9). (With the caveat of §5.4 on $A$-uniformity, and with P3-all$^*$ open.)
2. **No averaged analytic statement.** Robust ensemble inputs are capped (round 10 Thm B). "Almost all"
   exceptional sets swallow any $O(M)$-point selection (Prop C). Pointwise statements about Gaussian primes
   reconstitute the problem (Thm D).
3. **No single-place auxiliary polynomial on $X_M$** ($\tau^*<2$), and **no auxiliary polynomial on
   $\overline M_{0,k}$** for any $k$ (`threshold/L8.md` Thm A: the balanced class lies in the interior of the movable cone).
4. **No single linear form in logarithms** (Baker/Matveev, single-form Lang–Waldschmidt).
5. **No $S$-integer argument whose threshold depends on the prime set only through $|S|$** (lit18 Props A, B;
   **unrefereed**: proved in `lit18/check_tools.md` §2.2 and re-derived by the lit18 judge, `lit18/README.md`
   §2.2). It must use equality of norms, or a threshold tied to $\prod_{p\in S}p$.
6. **Route 5's determinant method at $T=1$** reduces, after Step 0, to a box statement that fakes may test;
   its decisive test is undecided with no method in either direction.
7. In addition (local smoothing assessment, `local_smoothing/README.md`): no $L^p$, energy, decoupling or
   restriction estimate, since these are capped at $\Theta(\log R)$ on a single circle.
8. Also (§5.4 row 16; `growth/auxiliary.md`, **unrefereed**, with the refereed `round10/additive.md` §5):
   no interpolation-determinant (Bombieri–Pila, Heath-Brown, Salberger), gap-principle, larger-sieve or
   Stepanov argument on one circle, since these reduce to pair data.

The judge's conclusion: a genuinely new mechanism must be **global, multi-place, and tied to the actual
arithmetic of $N$**. The only known statements of that type are Vojta-type inequalities, and their sharpest
arithmetic form in the record is $(L_8)$.

## 8. Index of files and referee corrections

### 8.1 Theorem → file

Paths are relative to `docs/claude_review/`. Each write-up `X.md` has its adversarial referee report
`referee_X.md` (or the round README's §2) beside it, and its checkers in `X_checks/`. The round READMEs give
the referee-corrected statements and are the authoritative summaries.

| result in this document | status | file and theorem |
|---|---|---|
| Thm 2.1, Cor 2.2: $3M\log M\le\log(R^2)+(14+4\log^+C)M$; $M\le(2/3+\varepsilon)\log R/\log\log R$ | proved | `two_thirds/proof.md` Thm 7.5, Cor 8.1 |
| Lemma 2.3: chain lemma | proved, sharp | `two_thirds/proof.md` Lemma 4.1; `check_chain_lemma.py` |
| exact finite inequality $F^*(M)\le\frac M4\log R+\binom M2\log C$ | proved | `two_thirds/proof.md` Thms 5.1, 5.2 |
| predecessor $(1+\varepsilon)$ bound | verified | branch `docs/inert_prime_cofactor_bound.md`; `README.md` §1 |
| $K\log R/\log\log R$, $K>8$ | paper (Lean development) | `paper/jarnik_sublog.tex`, main branch |
| Thms 3.1–3.4: pure / row-deleted Hadamard cores; flipped Walsh | proved | `growth/walsh-extraction.md` Thms B, B'', C |
| Thm 3.5: cotangent trade-off | proved | `growth/cotangent.md` Thm A (constants per `growth/referee_cotangent.md` §3) |
| Thm 3.6, Cors 3.7–3.8: CC-scale rigidity; five-point diameter $\ge1.9433R^{2/5}$; five-point equivalence (attribution of the question to CG unverified) | proved | `growth/cotangent.md` Thm B |
| Thm 3.9: direction pinning (corrected) | proved | `growth/cotangent.md` Thm C; `growth/referee_cotangent.md` §5 |
| Thms 3.10–3.11: smoothness, near-extremal structure | proved | `round8/stability.md` Thm 3.1, Cor 3.2, §§2–5 |
| Thm 3.12: rank reduction | proved | `round8/weakest.md` Thm 1.3 |
| forced divisibility supercritical; rigidity | proved | `round10/determinant.md` Thms 2.1, 4.1 |
| additive rigidity; pairs outside a matching | proved | `round10/additive.md` Thms 2.1–2.6, Prop 3.1 |
| lower bounds $M(C)\ge7,8,9,10$ ($C\ge5.2700,5.6740,7.8804,8.2394$); $M_\infty(C)\ge6,7,8$ ($C>4,5.5365,5.89771$) | proved | `growth/families.md`; `growth/referee_families.md` addendum, §8 |
| scan to $10^{12}$ (best $C$ for $k=3..8$; at most 3 resp. 4 points on $\frac12\sqrt R$ resp. $\sqrt R$ arcs) | evidence (author's computation; independently reproduced only to $10^{10}$) | `growth/families.md`; `growth/referee_families.md` §5, addendum |
| Thm 4.1–4.2: Vojta on $\mathrm{Bl}_YX_7$ ⇒ uniform; $R^\alpha$ arcs | conditional | `growth/conditional.md` Thms 5.1, 5.2 |
| non-density equivalence; Ru–Vojta reproof of CC | proved | `growth/conditional.md` Prop 3.2, Props 6.1–6.2 |
| Thm 4.3: Vojta on $\overline M_{0,8}$ ⇒ uniform; $(L_8)$ | conditional | `growth/ptolemy.md` Thm D, §7 |
| $(L_6)$ false; no Liouville proof of $(L_k)$ | proved | `threshold/L8.md` Thms A, B |
| $(L_5)$ false | proved | `growth/referee_ptolemy.md` §3.3; `threshold/L8.md` Prop C |
| $\tau^*(M)=2-1/\lceil M/2\rceil$ | proved | `threshold/upper.md`; reduction in `threshold/verify.md` |
| Lemma 5.1, Def 5.2: gaps, angle models; Thm A.1, A.3 | proved | `round4/outside.md` |
| skew double core; fakes at $(8+o(1))M\log M$ | proved | `round5/exact.md` Thm 4.1(iii) (re-derived); `round8/lw.md` Thm 6.1 |
| fully flipped Paley fakes with point residues | proved | `round5/construct.md`, `round5/lower.md` |
| parity law; P3-char fakes; P3-all bracket | proved | `round5/round6_sharp/sharp.md` Prop 2.1, Thms 5.1, 6.4 |
| LW single-form barrier; fixed-$\varepsilon$ simultaneous fakes | proved | `round8/lw.md` Thms A, B, C |
| SLW false; QS framework | proved | `round8/referee_weakest.md` Props R1–R2; `round8/weakest.md` Thm 5.2, Prop 5.3 |
| analytic graft barrier; Prop C; Thm D | proved | `round10/analytic.md` |
| Thm 5.3: fake directions (unconditional fakes not at barrier scale; barrier scale needs (2.2)) | proved | `round12/fake.md` Thm 2.1; `round12/README.md` §2.4 |
| Thm 5.4: circle cap test; H$_{\rm light}$ | proved | `round16/barrier.md` Thm 3.5; `round17/hlight.md` Thm 2.8 |
| Props A, B (QS capped); QS$^{\rm prod}$ | **unrefereed** (Props A, B proved in the lit18 checker file and re-derived by the judge); QS$^{\rm prod}$ open | `lit18/tools.md`, `lit18/check_tools.md` §§2.2, 2.5; `lit18/README.md` §2.2 |
| interpolation determinants, Bombieri–Pila, gap principles, Stepanov reduce to pair data | unrefereed (no referee report) | `growth/auxiliary.md`; `growth/README.md` |
| conditional LW consequence: no genuine construct/lower/sharp fakes at Paley orders | conditional (LW), refereed | `round8/README.md` §1.1; `round4/flipped.md` Thm L |
| local smoothing assessment | assessment | `local_smoothing/README.md` |
| Thm 6.1: uniform margins + $H$ ⇒ growth | conditional | `round7/p3lower.md` Thm C |
| transfer theorem (rigidity + H$^{\rm gen}$ + L1) | conditional | `round14/rigtheory.md` Thm 4.1 (corrected) |
| L4$(K)$ refuted; forced/free boundary | proved | `round10/L4K.md` Thms 1, 2 |
| route-5 $T=1$ machinery, separation, MIN$^{\rm pair}$ reduction | proved implications | `round12/t1rigour.md`, `round12/nearstrata.md` Thm 4.1 |
| ZE ⇒ $O(\log\log R)$; ZE-col$^*$ ⇒ $O((\log\log R)^2)$ | conditional | `round10/determinant.md` Thm 5.7; `round11/subcluster.md` Thm 5.4 |
| $(R_H)$; anti-concentration vacuity | proved implications / proved | `round16/reduction.md`, `round16/antitheory.md` |
| equivariant witness bound; gap principle | proved | `round17/survive.md` Thms 3.1, 3.4 |
| $(L_7)$ partial results | proved / exact | `round4/L7.md`, `round7/L7.md` Thm B, Prop 5.4 |
| literature verdicts | sweep with checks | `lit18/README.md`, `lit18/direct.md`, `lit18/tools.md`, `lit18/primes.md` |

### 8.2 Corrections the referees made to earlier claims

Grouped by where the original claim was made. Each is recorded in the file cited.

**General bound and arithmetic notes** (`README.md` §1).
1. `inert_prime_residual_growth.md` §5: the remark that prime powers dilute the split gain to
   $m^2/(2(e+1)(p-1))$ is true for the residual product alone but **false as an obstruction** to the final
   coefficient. The chain lemma (Lemma 2.3) shows why.
2. The "1,035 configurations / 7,859 checks" sentence of `inert_prime_cofactor_bound.md` names no checker
   and could not be audited; a new sweep replaces it.
3. `four_row_one_coordinate_divisor_reduction.md` §7: "Equivalently" should read "Consequently"; data show
   the statement is strictly weaker than (17).
4. `cyclotomic_nullity_four_counterexample.md` §5: calling non-squarefree compression open is stale. The
   five-row family works for every odd $k\ge19$, not only $k\equiv0,1\bmod3$.
5. The (H) sharpening: the least admissible $k$ is $\lim_{C\to0^+}M^*(C)\ge5$; the anchor-only lcm is wrong
   in 4,362 of 19,064 sub-tuples (the all-edge lcm is right).

**The $2/3$ proof** (`two_thirds/proof.md` review box). Lemma 7.1 step 3 should be non-strict at $y=1$;
the floating-point claim for $\sigma,a_3,a_1$ was confirmed by recomputation with 34-digit logarithms; the
tail bound for $a_1$ uses a factor $1/(P-1)$ not stated in the text; "2/3 is the limit" needs
$F^*(M)\le\frac38M^2\log M+O(M^2)$; "only $\ell<M$ contribute" should read $\ell\le M-2$; $F^*(48)$ rounds to
1686.6. None affects the result.

**Growth write-ups** (`growth/`).
6. `growth/cotangent.md` Thm C: the product identity omitted the rational factor $m=\prod p^{\min(h,M-h)}$ (fails on
   34,164 of 100,256 even windows). Thm A: $\sqrt{192}$ improves to $\sqrt{24}$, and one split prime suffices.
   The five-point ratio "never below 4.9475" is false ($N=25$ gives 4.2024); it holds only in the searched
   regime. Deficit lemma condition $(j+2)c\le M/2$ should read $(j+1)c\le M/2$.
7. `growth/walsh-extraction.md`: Theorem B needs $M\ge28$; "the barrier is *exactly* non-Walsh cores with almost all
   rows flipped" is not proved; most of Prop P restates items 334/391–393; Theorem C's growth part follows from
   items 278 + 333; the comparisons with items 276/277 understate them (item 276 already gives an effective
   uniform bound for pure cores of any Hadamard type when $C\le\sqrt2$; see §3.1).
8. `growth/ptolemy.md` §7: "$(L_k)$ is false for $k=5,6$ (Theorem C)" was over-scoped, since $(L_k)$ allows a closed
   exceptional set; $(L_5)$ was repaired by the referee, and $(L_6)$ was later refuted in `threshold/L8.md`.
   "Vojta implies the uniform theorem" was already in `growth/conditional.md`.
9. `growth/families.md`: an incomplete symmetric search missed a primitive 10-point cluster ($C\approx8.239$), found
   by the referee; a proof gap in Corollary S was closed by computation.
10. `growth/conditional.md`: minor overstatements about the $\beta$ numerics and a loose converse in Prop 3.2(ii).
11. `paper/jarnik_sublog.tex` line 332 says four points on arcs $o(R^{2/5})$ are open; the notes call this a
    corrected mistake, and the Pell family at $R^{1/3}$ contradicts it (`growth/referee_cotangent.md` side note).

**Threshold round.** `threshold/tau_problem.md`'s "shell decoupling" claim was wrong (`threshold/verify.md`); the conclusions use
the full-slice invariant. The summary "forced divisibility is exactly critical" holds only for a form fixed
before the cluster (`round10/determinant.md` Thm 2.1).

**Rounds 4–6.**
12. The literal strengthening of `round4/outside.md` A.3 (factor $\mathrm{lcm}(1..2M)$ for every character) is not
    valid for actual circles; it fails on $x^2+y^2=5$. The valid class uses actual collision moduli.
13. `round5/construct.md` and `round5/lower.md` fakes carry **point-level** residues only; they violate the prime-level parity
    law (the construct $M=12$ fake does). Axiom (R) must be restricted to $v(c)\ne0$ (counterexample on
    $N=65$).
14. `round5/round6_sharp/sharp.md`: the aggregated height condition is false for actual primes (8976 of 8977 primes below
    $2\cdot10^5$ violate it); the per-unit condition repairs Prop 6.5.

**Round 7.** The `round7/p3G.md` headline "P3-G has no fakes at $O(M\log M)$" did not survive (selected prime sets
evade every computed closure); the union-bound completion sum is a construction tool, not a membership
condition; the glosses "PSD = Perron" and "parity gains at most a factor 2 per modulus" were false as
phrased; `round7/p3data.md`'s "fitted profiles are infeasible for their own residues" was an artefact of fixing
windings; `round7/p3Dlower.md`'s evidence table was withdrawn (mixed units).

**Round 8.** `round8/lw.md` called the max form "the Lang–Waldschmidt conjecture" (the standard one is the
product form, incomparable at equal $\varepsilon$); its error bound $W/(2\log M)$ was false. `round8/weakest.md`'s
SLW hypotheses are false, so its conditional theorems are vacuous; Prop 6.1's threshold is $(m_0-6)/2$; the
support-lemma minimum is $(q+3)/2$, not $(q+7)/2$. `round8/stability.md`'s crude forms are due to
Cilleruelo–Córdoba.

**Round 9.** "L4 is false up to norm 2; inverse theorems must use norm $\ge3$" withdrawn; Cor D.2's clause
wrong; the unweighted wording of L5 not known to be valid for actual circles.

**Round 10.** The L4K threshold "$K\ge2L(h_0+h)$, labelled proved" withdrawn; `round10/analytic.md` misclassified
arithmetic fibre inputs and the sharp large sieve; `round10/additive.md`'s AG height $M^{c_0/2}$ corrected to
$M^{1+c_0/4}$, and its "fake-violated input" headline rejected (violation is generic only).

**Rounds 11–13.** ZE-braid, SL, SL$^*$ (judge J1) and ZE-col are false; the cheap-span ZE-TC is false;
"ZE-TC holds in $M8$" is true only in the deg form; the round-12 fake computations were zonotopal or
optimiser artefacts; Prop 5.5 "scales compatible from $n\ge5185$" withdrawn; round 13's contribution
criterion for locality is false.

**Round 14.** Rig$(K)$ as posed used the wrong class (per-unit budget drops the single-angle coupling);
the transfer theorem needs the associate-twist correction.

**Round 16.** "AC is strictly stronger than MIN$^{\rm pair}$", "all probabilistic versions are the same
up to $e^{O(d)}$", "$2m'<q$ throughout the regime" (first false at $q=3251$) and "the circle case is one
explicit algebraic question" withdrawn; regime bound $m'\le8j+8$, not $8j+4$.

**Round 17.** `round17/cap.md` inverted necessary and sufficient ($\kappa_V\le T^*$ is only sufficient for a
cap); `round17/survive.md` Remark 4.4's $+1823$ is $\ge2848$; `round17/calib.md`'s "robust no in $\kappa$" should read
"$\kappa_V$, model-dependent".

**lit18.** `round8/README.md` §§2.3, 2.5: $QS(m,\delta,o(s))$ is false for every $m$, not open for $m\ge4$;
Gabdullin's result concerns arcs near the diagonal and corresponds to general-position arcs $R^{2\gamma-1}$;
Letendre's specialisation sits in a trivial window; Ziegler's count is $\ge s$, so even a valid application
would be three times worse than the proved bound.

---

## Appendix: Checker disagreements

An adversarial accuracy check of an earlier draft of this document reported 4 major and 23 minor
problems. All of them were re-checked against the sources and applied. **No finding was rejected.** The
re-check did turn up three places where the record's own files disagree with each other, or where the
earlier draft was following a summary sentence of the record. They are listed here so that a reader who
finds the other wording in the record knows which version this document follows and why.

1. **Growth-rate summary sentence (finding 1).** The earlier draft copied the record's top-level
   sentence. `README.md`, growth-rate status: "No growth improvement over `log R/loglog R` is proved,
   either unconditionally or from a standard conjecture." The same paragraph continues: "The only
   conditional route known to work from a standard conjecture is Vojta's conjecture for one fixed variety,
   which gives the uniform bound." The record means "no intermediate improvement". This document now says
   so explicitly (Summary item 1). This is agreement with the checker, not a disagreement.
2. **Fixed-$\varepsilon$ Lang–Waldschmidt constant (finding 5).** `round8/README.md` gives two versions.
   §1.1 Thm C and the §2.2 table, row 7: "`(max(8, 4/eps) + o(1)) M log M` ; `2/eps + o(1)` by the
   referee's refinement". The summary table near the end of the same README: "fakes at
   `(max(8, 2/eps) + o(1)) M log M` satisfy it". The earlier draft followed the second. This document
   follows the first, which is the theorem statement (§5.4 row 12), and agrees with the checker.
3. **"Cut slack" wording in §2.1 (finding 17).** The checker's replacement calls $\frac M4\log R$ "the cut
   slack". `two_thirds/proof.md` §10 reads: "The method gives a total lower bound of
   $(\frac14+\frac18)M^2\log M=\frac38M^2\log M$ (inert $+$ split) against the diagonal cut slack
   $\frac M4\log R$. Configurations with $M=o(\log R/\log\log R)$ are not excluded." In §2.2 of this
   document "cut slack" means the different quantity $Q$, so §2.1 uses the source's phrase "diagonal cut
   slack" and adds "(the right side of (UB))" to avoid a clash. The substance agrees with the checker.

Nothing else in the check was disputed.

---

## References (classical and literature, as cited in the record)

[CC92] J. Cilleruelo, A. Córdoba, *Trigonometric polynomials and lattice points*, Proc. AMS 115 (1992)
899–905. [CG07] J. Cilleruelo, A. Granville, *Lattice points on circles, squares in arithmetic progressions
and sumsets of squares*, CRM Proc. Lecture Notes 43 (2007) 241–262 (arXiv math/0608109). [CG09] J.
Cilleruelo, A. Granville, *Close lattice points on circles*, Canad. J. Math. 61 (2009) 1214–1238. V. Jarník,
*Über die Gitterpunkte auf konvexen Kurven*, Math. Z. 24 (1926) 500–518. J. Cilleruelo, *Arcs containing no
three lattice points*, Acta Arith. 59 (1991) 87–90. J. Bourgain, Z. Rudnick, Ann. Henri Poincaré 12 (2011)
1027–1053, and arXiv:1012.3843, 1105.0018, 1308.4247. P. Vojta, *Diophantine Approximations and Value
Distribution Theory*, LNM 1239 (1987). S. Keel, J. McKernan, on the effective and movable cones of
$\overline M_{0,n}$. E. Bombieri, J. Pila, Duke Math. J. 59 (1989). D. R. Heath-Brown, Ann. of Math. 155
(2002). Venues and theorem numbers marked "S" or "from memory" in `lit18/README.md` and `round8/README.md`
were not re-verified.
