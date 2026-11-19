# Round 31, Part IV -- judge report: integral fakes, the trace seed, and the orbifold-BMY bound (IV.2)

Source: Part IV of `astra_free_IV/answer.md` (heading "# Part IV", lines 1011-1136; scripts `iv_trace_check.py`,
`iv_section_check.py`; question `reply4.txt`, Astra's summary `final4.txt`). Seven referees each took one block:

| key | report | scope |
|---|---|---|
| bmy | `ref4_bmy.md` | proof of (IV.2): Langer's theorem, effectivity, normalisation at infinity, ledger, final algebra |
| examples | `ref4_examples.md` | (IV.2) tested by direct intersection theory on explicit surfaces; extremal configurations |
| trace | `ref4_trace.md` | the characteristic-2 trace seed and the symmetry exclusion |
| conseq | `ref4_conseq.md` | specialisation lemma, 48-mark/Turán corollary, 34/35, "where this stops" |
| bridge | `ref4_bridge.md` | dictionary from the record's function-field clusters to (IV.2) data |
| lit | `ref4_lit.md` | Langer Thm 0.1 and its predecessors; literature on sections; char p; arithmetic analogue |
| arith | `ref4_arith.md` | arithmetic substitutes for BMY; the F3 filter |

Judge scripts are in `judge_IV/` (`checks.py`, `dict_check.py`, `dict_iff.py`). The Part I-III judge reports are
`README.md`, `README_II.md`, `README_III.md`. **Every claim labelled [J] is the judge's own and is unrefereed.**

Labels as in Part III. "Proved" means a complete proof checked line by line by a referee, or an exact computation a
referee reproduced with independent code. "Computed" means exact data reproduced independently. "Conditional on
the citation" means the only unchecked input is the wording of a published theorem that no one could read in the
original (see §3(a)). **No referee found an error in (IV.2), its proof, or the function-field corollary.** The FALSE
and GAP verdicts in the summaries concern routes the referees posed themselves (arith: Bogomolov, Hodge index,
arithmetic Bogomolov, p-derivations; trace: odd characteristic; conseq: several parameters). None of them is a
FALSE or GAP on a claim Astra made in the form she made it. I rule on each below.

## 0. Verdict in five lines

1. **(IV.2) is proved** over $\mathbb C[t]$, with Astra's constants $3/4$, $24$, $m(m+3)$, $4$, for every $m\ge3$ and
   every real $\delta$ (its sign is never used). It rests on Langer's orbifold BMY (PLMS 86 (2003), Thm 0.1), whose
   wording was confirmed only through secondary sources (arXiv is unreachable from here). Two referees recomputed the
   Chern numbers directly on about 70 surfaces and found the ledger exact, with slack $\ge0$ throughout.
2. **Astra missed her strongest consequence, and it is a theorem of the record's model.** Applied to genuine
   clusters of the round-29/30 function-field model over $\mathbb C(t)$, $\mathbb R(t)$, $\mathbb Q(i)(t)$, (IV.2)
   gives **$M(D,0)\le24$ for every $D$ ($\le23$ for odd $D$) and $\sup_DM(D,c)<\infty$ for every $c$**. Three
   referees found this independently and the judge rechecked the dictionary. Status: $8\le\sup_DM(D,0)\le24$. The
   record's "OPEN" label and its "truncated function-field Vojta" pessimism must change.
3. **The characteristic-2 trace seed is correct**: $q-1$ marks, $\delta=E=0$, constant cofactors. It violates (IV.2)
   for $q\ge32$. So any proof of a uniform bound for integral fakes must use characteristic 0 or arithmetic
   collision gains (filter F3). The symmetry exclusion is correct only from $d\ge15$ (threshold never stated).
4. **Nothing changes for the integer problem.** The growth rate stays $(2/3+\varepsilon)\log R/\log\log R$. There
   is no uniform bound. Questions A and B and route 5 are untouched. Only one-parameter rational-branch polynomial
   templates are excluded beyond 24 points; the unit (golden/Pell) templates that hold the records sit at
   $\delta=L/8$, where (IV.2) is vacuous. The arithmetic analogue is an arithmetic orbifold BMY (Parshin /
   Moret-Bailly type), at least Szpiro strength. No unconditional substitute exists.
5. **Score for Part IV on the goal: 2/10.** It is the round's only real theorem, and it closes the model's main
   open question, but that model was already closed for internal work and its transfer to $\mathbb Z$ is negative.

## 1. Referee-corrected statements per block

### 1.1 bmy: proof of (IV.2) (referee: no error; hypotheses can be weakened)

**Proved (conditional only on the citation, §3(a))**
* **(IV.2), corrected statement.** *Let $N,n_i\in\mathbb C[t]$, $L=\deg N>0$, $n_i\mid N$, $\alpha_i=\deg n_i$, and let
  $P_1,\dots,P_m\in\mathbb C[t]^2$, $m\ge3$, have coprime entries and pairwise nonzero minors $D_{ij}=[P_i,P_j]$. Put
  one layer $h$ for each root of $N$ and each depth up to its multiplicity, with $s_h$ the number of $n_i$ vanishing to
  that depth there. Put $E_{ij}=\deg\gcd(n_i,n_j)-\deg\gcd(n_i,n_j,D_{ij})$, and suppose $\deg D_{ij}\le(\alpha_i+\alpha_j)/2
  -L/4+\delta$ for all $i\ne j$, for some real $\delta$. Then*
  $$\tfrac34(m-24)L+\sum_h(s_h-m/2)^2\le m(m+3)\delta+4\sum_{i<j}E_{ij}.$$
  Corrections: $m\ge3$ suffices (Astra: $m\ge4$); $\delta\ge0$ is not used.
* **Langer's hypotheses hold.** The pair $(Y,B)$ has simple normal crossings on its positive-weight part (checked in
  local graph coordinates). The coefficients are $2/3$, $1/2$, $0$, all $<1$, so the pair is klt, hence lc.
* **Effectivity.** $K_Y+B$ is $\mathbb Q$-effective: $K_{Y_0/\mathbb P^1}+C_i+C_j\sim(C_i\cdot C_j)F$; the discrepancy
  with the strict transforms is $1-l$; the invariant $\mathrm{mult}_qV_{ij}\ge(C_i\cdot C_j)_q$ keeps every coefficient
  $\ge0$ (machine-checked on all examples). Averaging gives horizontal coefficient $2/m\le a=2/3$ (needs $m\ge3$); the
  four half-fibres cancel $\pi^*K_{\mathbb P^1}$.
* **Normalisation at infinity, (IV.3), $2I=(m-1)S$.** Elementary transformations at a point of the infinity fibre
  carrying $k>m/2$ sections reduce $S$ by $2k-m$ and terminate because $S\ge0$; the cyclic pairing by $\lfloor m/2\rfloor$
  gives (IV.3).
* **Ledger.** $F(k,0,0)=4-8k/9$, $F(2,0,a)=8/9$, $F(2,a,a)=-4/9$, chain total $\le4e/9$ (in fact $(12-8e)/9<0$ after the
  charge), $\phi(k)=\max(0,\binom k2+2k-9)$, $\phi(u+v)-\phi(u)-\phi(v)\le2uv$ (defect $uv+9$ for $u,v\ge3$).
* **Final algebra.** The $S$ coefficient is $2(m+3)/9>0$, so (IV.3) may be substituted. The per-layer term is exactly
  $4-m/6-\frac29(s_h-m/2)^2$; times $9/2$ this is (IV.2). [J] reproduced symbolically in `judge_IV/checks.py`
  (`per-layer identity ok: True`; $\frac92(4-m/6)=-\frac34(m-24)$; $\delta$ coefficient $m^2+3m$; $E$ coefficient 4;
  $\phi$ defect $\le2uv$ for all $u,v<80$).

**Computed.** A lattice-level engine ($\mathrm{Pic}(Y)$, $K_Y$, $B$, $e_{\rm orb}$, $(K_Y+B)^2$ computed from classes, not
from the ledger) reproduces the ledger exactly on about 40 random configurations, pencils (slack exactly $4r$), nested
trees, 50-blow-up chains and genuine record clusters.

**Not checked.** That Langer's local $e_{\rm orb}$ at SNC points with arbitrary weights is $1-a$ and $(1-a)(1-b)$.
[J] With those local values the global formula $e(Y)-\sum a_ie(D_i)+\sum a_ia_j(D_i\cdot D_j)$ follows by inclusion-
exclusion. The values are classical for standard weights, and Pokora uses them for real weights.

### 1.2 examples: direct tests of (IV.2) (referee: no error; one script mislabelled)

**Proved / Computed**
* (IV.2) and the ledger agree with direct intersection theory on 28 exact configurations over $\mathbb Q$,
  $\mathbb Q(i)$, $\mathbb Q(\sqrt{-3})$. The cases include two shared centres at infinity, 7- and 9-fold points with
  nested groups, prime powers with tangency chains $e=2,3$, and groups left at infinity. Direct value = ledger
  $\le$ bound every time.
* $K+B$ is effective but **not nef** in several examples (e.g. $(K+B)^2=-5$ in the family $(c+u_iN,1)$). So the
  $\mathbb Q$-effective form of Langer's theorem is needed. Miyaoka/KNS/Megyesi (nef hypotheses) would not suffice.
* Family $(c+u_iN,1)$ with $n_i=N$: exact slack $18L$ for every $m$ ($m=6,12,30$ checked).
* **New extremal data (no norm condition):** a 12-mark configuration with $\delta=0$ and $E=0$, exact over
  $\mathbb Q$ at $L=8$: a coset of $\mathbb Z/2\times\mathbb Z/4$ on $y^2=x(x+4)(x+25)$ embedded as a $(2,2)$-curve in
  $\mathbb P^1\times\mathbb P^1$ (a cube on a quadric with its 12 four-point planes). Also a 9-mark one (dual Hesse,
  $L=12$). So **the maximum number of marks in the no-norm $\mathbb C[t]$ problem at $\delta=E=0$ lies in $[12,24]$.**
* $a=2/3$ is optimal for this scheme and exceptional weights cannot help balanced layers, so 24 is the method's limit.
  (lit and arith confirm: the per-layer term is $4-\frac m4(2a-\frac32a^2)$; $a=1$ gives 32, $a=1/2$ gives 25.)

**Corrected.** Astra's "computes both Chern numbers directly" for $y=u_i-t^2/u_i$: `iv_section_check.py` only evaluates
the closed-form ledger, on an unnormalised surface. That is valid but it is not an independent check. The referee's
direct computation on the normalised surface gives $Q=(4m^2+12m)/9$, in agreement with the ledger.

### 1.3 trace: the characteristic-2 seed and the symmetry exclusion (referee: no false claim)

**Proved**
* **Trace seed.** For $q=2^d$, $d\ge3$, marks $a\in\mathbb F_q^*$, primes $t-x$, $x\in\mathbb F_q^*$, $x$ marked for
  $a$ iff $\mathrm{Tr}(ax)=0$: $f_a(x)=x^{q/2-1}$ at marked $x$, and $f_a-f_b$ has degree exactly $q/4-1$, equal to the
  number of common marked $x$. Hence $D_{ab}=$ nonzero constant $\times$ the product of all common marked primes.
  Parent weight $(q-1)\log q=m\log(m+1)$; $\delta=0$ with slack $1/4$; $E=0$. Reproduced by Euclid gcds in
  $\mathbb F_{2^d}[t]$ for $q=8,\dots,512$. [J] By hand: $f_a=\sum_{j=0}^{d-2}a^{2^j-q/2}t^{2^j-1}$, leading
  coefficient of $f_a-f_b$ is $a^{-q/4}-b^{-q/4}\ne0$; two distinct trace hyperplanes meet in $q/4$ elements. Also
  rechecked by independent code for $q=8,16,32,64$ (`checks.py`).
* **$d_{ab}$ depends only on $a+b$** for any prime assignment; translations preserve every cut.
* **Symmetry exclusion, corrected statement.** *Under pair separation, an integer realisation of the trace incidence
  table at the III.7 normalised bound (any fixed $C_1$) would give a faithful action of $(\mathbb F_q,+)$ on
  $\mathbb P^1(\mathbb Q)$ by projectivities. This is impossible for $q>4$. The exclusion therefore holds for
  $q-2\ge16390$, i.e. $d\ge15$ (about $d\ge14$ with sharp counting); $d\le14$ is open.* Pair separation is a real
  hypothesis: it fails for some honest assignments (min separation$/W=0.48$).
* An elementary abelian 2-subgroup of $\mathrm{PGL}_2(\mathbb C)$ has order $\le4$ (Klein). More generally,
  $\mathrm{Aut}(d)$ embeds in $\mathrm{PGL}_2(\mathbb Q)$, so has order $\le12$ (finite-subgroup list quoted).

**Corrected interpretation**
* With $C_1=1$ the Part IV form of (IV.1) holds for every $q$; the literal Part III factor $1/4$ needs $q\ge256$. The
  cone needs a common shear by $t^H$, $H\ge q/4$. The constant field grows with $m$, so this is **not** a growing
  family over one function field.
* **New (referee):** the only quadratic forms $Q$ with $Q(P_a)\in\mathbb F_q^*n_a$ for all marks form the pencil
  $\lambda(tX^2+Y^2)+b(XY+t^{q/2-1}Y^2)$. $b=0$ is inseparable; every $b\ne0$ form is isotropic. With $\lambda=0$,
  $Z_a=(n_a,N/n_a)$ is a **genuine split-model (hyperbola) cluster** with $q-1$ points at $D=q/2$, $c=0$: $q-1$
  divisors of $t^{q-1}-1$ in a short interval. Verified to $q=512$. This is the function-field instance of the
  record's split calibration. No anisotropic (circle-type) norm is reached.

**Gap (not a claim of Astra's).** Odd characteristic: the $p$-linearised analogue misses the critical budget by
$q(p-2)^2/(4p^2)-1/4$ (excess 1/2, 2, 11, 6 at $q=27,81,125,49$). Neither existence nor non-existence is known.

**Novelty.** The incidence is the point-hyperplane design of $\mathrm{PG}(d-1,2)$ (the Fano plane at $d=3$; lit).
Its realisation by polynomial sections with constant cofactors is new to the record. The exclusion is new only for
the no-norm model: pure Sylvester/Walsh cores are already excluded for genuine clusters (`growth/walsh-extraction.md`
Thm B, Cor B.2-B.3).

### 1.4 conseq: "where this construction stops" (referee: correct, scope corrected)

**Proved**
* **Specialisation lemma, corrected statement.** *If fixed one-parameter polynomials in $\mathbb Z[T]$ (or
  $\mathbb Q(i)[T]$ with bounded denominators) supply a nonconstant parent, row factors and primitive columns
  satisfying (IV.1) at infinitely many integers $T\to\infty$, with $\log C_1=o(\log T)$, then the $\mathbb C[t]$ model has
  $\delta=0$ and $E=0$, so $m\le24$.* If $C_1=T^\kappa$ then $\delta=\kappa$, and 25 marks are excluded only when
  $\kappa<3L/2800$. The lemma survives subsequences and bounded discrepancies (e.g. bridge Thm 1.2's $g(t_0)$).
* **Sampling corollary.** On any 48 marks some edge has $\deg D_{ij}\ge(\alpha_i+\alpha_j)/2-L/4+L/136$. A fraction
  $\ge(m-47)/(47(m-1))$ of all pairs does (Turán). 49 is the optimal sample size: $75/10192\approx1/135.89$. [J]
  `checks.py` confirms $1/136$ at 48 and the maximum $75/10192$ at $k=49$ over $25\le k<200$. Averaging over
  37-subsets with the hypergeometric square term improves the global maximum excess to $(19/1480-o(1))L\approx L/77.9$
  (referee's, not Astra's).
* **(IV.4)** follows. It restates the degree bound for a presentation with a conversion factor $\lambda$.
* **Partial divisibility.** $\overline E\ge(1/140-1/(8(m-1)))L$; constant $(k-18)/(2k(k-1))$, maximal at $k=35,36$.
  With pair separation $\ge L/2$, at most $34/35+o(1)$ of the prescribed common degree is retained. Without pair
  separation only $139/140$. The best attainable fraction is open.

**Gap (scope; repairable only partially).** Families in several parameters: the proof uses one variable twice (size
$\sim|T|^{\deg}$ and division with remainder). $T_1\mid T_2^2+1$ on a Zariski-dense integer set is a counterexample to
the argument. Repair needs infinitely many good points on one polynomial curve along which the parent is nonconstant.
Astra wrote "one-parameter", so nothing she wrote is false.

**Not excluded** (referee's list, upheld): unit/Pell/golden templates (level $c=D/4$, $\delta=L/8$, where (IV.2) is
vacuous because $6(m-24)\le m^2+3m$ for all $m$; [J] checked: $m^2-3m+144>0$), multi-parameter families, polynomially
growing $C_1$, sporadic extra points, isolated parents, anything uniform in $R$.

### 1.5 bridge: the record's function-field clusters as (IV.2) data (referee: correct)

**Proved**
* **Dictionary.** Let $z_1,\dots,z_m\in\mathbb C[t]$ be a cluster of the round-29 model (`round29/ff.md` §1: distinct,
  degree $D$, $z_j\bar z_j=N$, $\deg(z_i-z_j)\le K=\lfloor D/2+c\rfloor$; bar acts on coefficients). Remove real
  roots/common factors (normalisation (N2); this only strengthens the condition). Choose $U$ with $U\bar U=N$,
  $\gcd(U,\bar U)=1$, and put $P_i=\gcd(z_i,U)$, columns $(\mathrm{Re}\,P_i,\mathrm{Im}\,P_i)\in\mathbb R[t]^2$,
  $n_i=P_i\bar P_i$. Then $L=\deg N=2D$, $E_{ij}=0$ for all pairs, the columns are primitive, and
  $$\deg(z_i-z_j)=D+\deg D_{ij}-(\alpha_i+\alpha_j)/2,$$
  so the cluster condition is exactly the degree bound of (IV.2) with $\delta=K-D/2$. Prime powers and non-squarefree
  $N$ need no special treatment. Exact sympy checks on the cube clusters ($D=8,16$), the vertical-line clusters
  ($M(7,1)\ge8$, $M(8,1)\ge14$) and two non-squarefree clusters (bridge, and independently conseq's `ffstar.py`, bmy,
  examples).
* [J] `dict_iff.py`: on all 496 pairs of the 32 points $z_a$ over a generic $D=5$ parent with Gaussian-rational roots,
  the identity $\deg(z_a-z_b)=D+\deg D_{ab}-(\alpha_a+\alpha_b)/2$ and $g_{ab}\mid D_{ab}$ hold with no failure.
  `dict_check.py`: the 6-point $D=4$, $c=0$ cluster (equal imaginary parts) satisfies the degree bound with
  $\delta=0$; at $D=6$ the cluster condition and the degree bound fail together.

**Corollary (bridge's exact form; see ruling §3(b)).** Upheld with one correction (below).

**Corrected.** "Template lower bounds for $M_\infty$ can never exceed 24" holds for one-parameter **polynomial**
(rational-branch) templates with bounded $C$. The record's 8-point records come from unit (golden/Pell) templates,
which are not polynomial in a parameter and sit at $\delta=L/8$. The two referees (bridge, conseq) do not disagree
once the scope is stated.

**Correct and important for the record.** The round-30 conditional "unbounded $\mathbb Q(i)(t)$-clusters would refute
$A_5(\delta)$" now has a false premise: that route is dead, and nothing about $A_5$ follows. The remark that the
uniform bound "has the shape of open truncated function-field Vojta" is withdrawn as an obstruction. That Vojta
statement is the shape of the sharp (R$_4$) only. The uniform bound comes from a different, known inequality.

**Not checked.** Literature: no prior uniform function-field bound for short arcs was found (WebSearch only; not
exhaustive).

### 1.6 lit: citations and context (referee: correct; originals unread)

* **Langer, PLMS 86 (2003) 358-396, Thm 0.1, as quoted by Pokora (arXiv 1612.05141):** $X$ normal projective surface,
  $D$ a $\mathbb Q$-divisor with coefficients in $[0,1]$, $(X,D)$ lc, $K_X+D$ $\mathbb Q$-effective $\Rightarrow$
  $(K_X+D)^2\le3e_{\rm orb}(X,D)$. Part IV uses exactly this. Predecessors (Miyaoka 1984, Kobayashi-Nakamura-Sakai
  1989, Megyesi 1999, Sakai 1980) need standard weights and/or nefness, which Part IV's pair does not have. [J] My own
  WebSearch confirms the paper (arXiv math/0012180), the abstract's "in the log canonical case we prove a BMY type
  inequality", and Pokora's quotation of $(K_X+D)^2\le3e_{\rm orb}$. It does **not** show the effectivity hypothesis
  verbatim; arxiv.org does not resolve from this session.
* **Lineage.** Hirzebruch-type inequalities via Langer (Bojanowski, Pokora; conic-line arrangements with tangency
  chains); Urzúa's arrangements of sections (Rend. Lincei 2011). No prior bound on the number of sections with
  prescribed contact was found. Yamanoi, moving-target SMT, Brownawell-Masser and the canonical-class inequality give
  nothing uniform in $D$. The technique is routine for a specialist; the application is new.
* **Characteristic $p$.** The trace incidence is the point-hyperplane design of $\mathrm{PG}(d-1,2)$, the same
  Artin-Schreier mechanism by which $\mathbb F_q$-line arrangements violate Hirzebruch. (IV.2) is homogeneous under
  $t\mapsto t^p$, so Frobenius base change alone cannot violate it.
* **New reading (referee).** At $\delta=E=0$, (IV.2) is equivalent to $\sum_hs_h(m-s_h)\ge(m^2/4+3m/4-18)L$, i.e. a
  lower bound for the finite part of $h_{K+\partial}$ on $\overline M_{0,m}$ at the primes of $N$. The arithmetic
  analogue has no exceptional set and an explicit constant, so it is not a case of Vojta. It is formally
  incomparable with Vojta $D=0$ on $\overline M_{0,8}$ and at least of abc/Szpiro strength.

### 1.7 arith: arithmetic substitutes (referee: none exists)

**Proved / correct**
* (IV.2) fails in characteristic 2 for $q\ge32$ (LHS $341/2$, $3717/2$, $19685/2$ at $q=32,64,128$ against RHS 0).
  [J] `checks.py`: LHS $=-175/2,-195/2,341/2,3717/2,19685/2$ for $q=8..128$, matching lit's values.
* **Any constant suffices for a bound.** With $(K+B)^2\le\kappa e_{\rm orb}$ the scheme gives $m\le2\kappa(\kappa+1)$
  ($\kappa=3$: 24). With chain weights for $\kappa>3$: $m\le40,60,144$ at $\kappa=4,5,8$. So the strength of BMY is
  irrelevant; its existence is what matters.
* **F3 filter (stated precisely).** Any inequality with a characteristic-$p$ analogue that the trace surface
  satisfies cannot bound integral fakes. This kills: Bogomolov for rank-2 bundles on rational surfaces whose
  semistability comes from the configuration (Shepherd-Barron 1991, secondary sources), arithmetic Hodge index
  (Faltings-Hriljac; the record's round-20 identification with the cut identity), Noether, Riemann-Roch, geometry of
  numbers. The referee's FALSE verdicts are on these referee-posed routes, and I uphold them as no-go arguments.
* The record's bound $W\ge3m\log m-O(m)$ transfers to no-norm fakes over $\mathbb F_q[t]$, with gain exactly 0 over
  $\mathbb F_{m+1}[t]$ (no residue field smaller than $m$). There is no contradiction with the trace seed. Its weight
  $m\log(m+1)$ is below the integer bound, so the trace table cannot live at its own weight over $\mathbb Z$.

**Gap (open routes, not errors).** Arithmetic Bogomolov (Miyaoka, Moriwaki, Soulé): its geometric shadow on a
$\mathbb P^1$-fibration has no characteristic-0 content; Miyaoka's semistable extension has no $\mathrm{Spec}\,\mathbb Z$
analogue (the missing Kodaira-Spencer class). p-derivations and the archimedean derivative are F1/pair-capped. Not
repairable with known tools.

### 1.8 Withdrawn or corrected in Part IV (collected)

1. (IV.2): $m\ge3$ suffices; $\delta$ may be any real (negative $\delta$ is used for odd $D$).
2. "Computes both Chern numbers directly" (`iv_section_check.py`): it evaluates the ledger; the independent direct
   computation is the referees'.
3. Symmetry exclusion: "sufficiently large $q$" means $d\ge15$ ($q-2\ge16390$); $d\le14$ open; needs pair separation.
4. Trace seed: the Part III normalisation (factor 1/4) needs $q\ge256$; the constant field grows with $m$; it is a
   genuine split-form cluster, isotropic, not a circle analogue.
5. 48 marks / $1/136$: correct but not optimal ($k=49$: $75/10192$; global $\approx L/77.9$ by 37-subset averaging).
6. Specialisation lemma: one parameter only; $\log C_1=o(\log T)$ required.
7. Omission: the function-field corollary (§3(b)), Astra's strongest result, is not stated in Part IV.

## 2. What is genuinely new relative to the record, and correct

* **Inequality (IV.2).** New. The record's only BMY computation is `round4/L7.md` (a universal-family check on
  $\overline M_{0,7}$, no obstruction); grep finds no other occurrence of Langer, Miyaoka or BMY. It is the first
  exceptional-set-free, D-uniform inequality in the record that genuine function-field clusters satisfy.
* **The uniform bound in the function-field model** (§3(b)). New, and it answers the question rounds 29-30 left open.
* **Characteristic-2 trace seed** and the identification of the seed as a genuine split-form cluster with $2D-1$
  points at $c=0$ (trace referee). New to the record (no hits for linearised or subspace polynomials). Elementary.
* **F3 filter** and the "any $\kappa$ gives $m\le2\kappa(\kappa+1)$" remark (arith). New framing; F3 is a useful
  standing caution.
* **Extremal no-norm data**: 12 marks (elliptic $\mathbb Z/2\times\mathbb Z/4$ on a quadric) and 9 marks (dual Hesse) at
  $\delta=E=0$, against the record's 8 (examples). New.
* **Reading of (IV.2) as a $K+\partial$ statement on $\overline M_{0,m}$** (lit). New to the record.
* Minor: specialisation lemma, sampling and Turán corollaries, $34/35$, the symmetry exclusion for $d\ge15$.

Not new: the Plücker/star dictionary (`growth/ptolemy.md` Thms A, B); exclusion of pure Walsh cores for genuine
clusters (`growth/walsh-extraction.md`); "the arc condition over $\mathbb C(t)$ is linear" (round 29 Lemma 1.1).

## 3. Rulings

### (a) Is (IV.2) proved, and with which constants?

**Yes**, conditional only on the citation. Statement as in §1.1: over $\mathbb C[t]$ (hence over any field of
characteristic 0, by the Lefschetz principle), for $m\ge3$ and any real $\delta$,
$\frac34(m-24)L+\sum_h(s_h-m/2)^2\le m(m+3)\delta+4\sum_{i<j}E_{ij}$. Constants exactly Astra's. Consequently
$m\le24$ when $\delta\le0$ and $E=0$, and $m\le23$ when $\delta<0$ and $E=0$ ([J]: at $m=24$ the right side is
$24\cdot27\,\delta<0$ while the left side is $\ge0$).

The basis for "proved": two referees checked every step independently (bmy by hand plus a class-level engine; examples
by direct intersection theory on 28 exact configurations). I redid the final algebra symbolically [J]. The one
unverified input is the wording of Langer's Theorem 0.1. Three referees and the judge confirmed it only through
secondary quotations (Pokora 1612.05141): lc pair, coefficients in $[0,1]$, $K+D$ $\mathbb Q$-effective. The effective
(not nef) form is necessary here (examples: $K+B$ is not nef in several cases). My own recollection [J] agrees with
the effective form, and Pokora's papers apply it to arrangements with non-nef $K+D$ in the same way. **Classification:
PROVED modulo a published theorem whose statement has been checked against secondary sources only.** The record should
carry that caveat until someone reads PLMS 86 (2003) 358-396 directly. The 24 is the method's limit ($a=2/3$ optimal);
it is not shown to be sharp. In the no-norm problem the truth lies in $[12,24]$.

### (b) Is a bound for $M(D,c)$ uniform in $D$ now a THEOREM?

**Yes**, with the same caveat. Exact statement for the record:

> **Theorem (round 31, Part IV + referees bmy, examples, conseq, bridge).** Let $k=\mathbb C$, $\mathbb R$ or
> $\mathbb Q(i)$ (any subfield of $\mathbb C$ stable under conjugation), and let $M(D,c)$ be the largest number of
> distinct $z_j\in k[t]$ of degree $D$ with $z_j\bar z_j=N$ (one fixed $N$, any multiplicities) and
> $\deg(z_i-z_j)\le D/2+c$. Then, with $\delta=\lfloor D/2+c\rfloor-D/2\le c$:
> (i) $M(D,0)\le24$ for every $D$, and $M(D,0)\le23$ for odd $D$; more generally $M(D,c)\le24$ whenever $0\le c<1/2$;
> (ii) $M(D,c)\le m-1$ whenever $3(m-24)D>2m(m+3)\delta$; in particular $M(D,c)\le48$ for $D>5096c/75$,
> $M(D,c)\le47$ for $D>68c$, and $M(D,c)\le24$ for $D>1400c/3$;
> (iii) hence $\sup_DM(D,c)\le\max\big(47,\;2^{2c+2}(\lfloor68c\rfloor+2)\big)<\infty$ for every $c\ge0$.
> Together with the cube family: $8\le\sup_DM(D,0)\le24$. The proof uses Langer's orbifold BMY (characteristic 0);
> it says nothing over $\mathbb F_{q^2}(t)$.

[J] I rechecked the thresholds in `checks.py` ($1400/3$, $68$, $5096/75$) and the dictionary on two further examples
(§1.5). The record's Theorem A gives $M(D,0)\le4(D+2)$, so (i) is better for every $D\ge5$. (R$_4$) (rank $\le4$,
$\limsup M=8$) stays **open**. The code rank is now $\le23$ at $c=0$.

**Round 30's status must change.** "Bounded rank (R$_4$) stays open" remains true. But "the problem stays open" (round
29 §0 item 5, expert note §5 item 10), "uniform-in-$D$ bound OPEN", and "(R$_4$)/the uniform bound has exactly the shape
of truncated function-field Vojta, which is itself open, so this is calibration" must be corrected. The **uniform
bound** is now a theorem; only the **sharp constant** has the Vojta shape. The round-30 conditional "unbounded
$\mathbb Q(i)(t)$-clusters would refute $A_5(\delta)$" now has a false premise. The decision "close the line for
internal work" stands, now on firmer ground.

### (c) The finite-characteristic construction and the filter it gives

**Correct** (§1.3), as a no-norm integral fake over $\mathbb F_{2^d}(t)$ with $m=q-1$, $\delta=E=0$, constant
cofactors and weight $m\log(m+1)$, for every $d\ge3$. The constant field grows with $m$. In fact it is a genuine
**split-form** cluster (hyperbola), isotropic, and no anisotropic norm fits it. **New filter (F3):** *a proof of a
uniform or growth bound for integral fakes, or for anything that does not use the exact norm, must use either
characteristic 0 (a cotangent/Kodaira-Spencer input, as BMY does) or collision gains from residue fields smaller than
$m$ (as the record's pigeonhole does). Every inequality whose characteristic-$p$ analogue the trace surface satisfies is
excluded.* This kills Bogomolov-type (rational surface), Hodge-index, Riemann-Roch and geometry-of-numbers routes for
integral fakes (arith). The symmetry exclusion of the trace table over $\mathbb Z$ is correct for $d\ge15$ and is a
five-line corollary of III.7. Its target (the pure Walsh core) is already excluded for genuine clusters.

### (d) What follows for the integer problem

Bluntly:
* **Growth rate:** unchanged, $(2/3+\varepsilon)\log R/\log\log R$. The constant 3 in $W\ge3m\log m-O(m)$ is untouched.
* **Uniform bound $M(C)$:** unchanged. Nothing is proved for any individual integer parent. Astra says so herself
  (summary item 11), and no referee found a transfer.
* **Polynomial templates:** one-parameter rational-branch templates in $\mathbb Z[i][u]$ (or $\mathbb Q(i)(t)$ with
  bounded denominators) with bounded $C$ (or $\log C=o(\log u)$) carry at most **24** points, uniformly in $D$ and $C$.
  This is the weakest mechanism (`round30/bridge.md`). **Not covered:** unit/golden/Pell templates, which hold the
  record 8-point families and sit at $\delta=L/8$ where (IV.2) is vacuous; multi-parameter families; polynomially
  growing $C$; sporadic points; isolated parents.
* **Question A ($(L_8)$):** unchanged. (IV.2) has no content for $m\le8$ (bridge: it holds automatically for $m\le7$,
  and for $m=8$ without trivial layers).
* **Question B ($(L_5^\circ)$, $A_5(\delta)$):** unchanged. One conditional route to refuting $A_5$ (unbounded
  function-field clusters) is now dead.
* **Route 5:** unchanged.
* **Integral fakes over $\mathbb Z$ at growth scale:** still open, in both directions.

### (e) The arithmetic analogue

It is an **arithmetic orbifold Bogomolov-Miyaoka-Yau inequality** for $\mathbb P^1_{\mathbb Z}$ blown up along the
contact tree of $m$ sections, with boundary $\frac23\sum C_i$ and layers weighted by $\log p$. It is of the Parshin
(1988) / Moret-Bailly (1990) type and implies Szpiro-type statements, so it is at least of abc/Szpiro strength. The
naive arithmetic BMY is false, and Miyaoka's 1988 claim was withdrawn (secondary sources). Equivalently (lit): a lower
bound for the finite part of $h_{K+\partial}$ on $\overline M_{0,m}$ with no exceptional set. That is formally
incomparable with, and not weaker than, the record's Vojta $D=0$ on $\overline M_{0,8}$. **The arith block found no
viable unconditional replacement.** Any constant $\kappa$ and an $O(W)$ error would suffice, which sharpens what is
missing. But every proved Arakelov tool fails F3. The missing piece is the Euler-number/cotangent side, since
$\mathrm{Spec}\,\mathbb Z$ has no $\Omega^1$.

## 4. Sentences of the record that should be updated (recommendations only; the record is read-only)

1. `expert_note/expert_note.md` §5 item 10, bullet *Function-field model*: "and the problem stays **open** there:
   bounded clusters are equivalent to bounded rank of the ratio lattice" $\to$ "the uniform question was open there
   until round 31, Part IV: $\sup_DM(D,c)<\infty$ for every $c$ and $M(D,0)\le24$ over $\mathbb C(t)$, $\mathbb Q(i)(t)$
   (orbifold BMY; item 12)."
2. Same item, point 5 of `round29/ff.md` §0 as quoted in the note ("Proving it is ... the function-field Lang-Vojta
   exceptional-set problem"): restrict this to (R$_4$), the sharp rank-4 statement.
3. §5 item 11, bullet *Rank*: after "(R$_4$) stays **open**. It has exactly the shape of truncated function-field
   Vojta ... so this is calibration and not evidence." add: "This applies only to the sharp constant. A uniform bound
   follows from Langer's orbifold BMY on the section surface (round 31, Part IV): $8\le\sup_DM(D,0)\le24$, with
   $M(D,0)\le23$ for odd $D$, and $\sup_DM(D,c)\le\max(47,2^{2c+2}(\lfloor68c\rfloor+2))$. It is characteristic 0 only,
   and Langer's Theorem 0.1 has been checked through secondary sources only."
4. §5 item 12 (round 31 summary): add a Part IV bullet: (IV.2) and the function-field theorem as in §3(b); the trace
   seed (characteristic 2, a genuine split-form cluster with $2D-1$ points, isotropic) and the filter F3; "nothing for
   integer parents; the arithmetic analogue is an arithmetic orbifold BMY of at least Szpiro strength; one-parameter
   polynomial templates carry $\le24$ points; unit templates are untouched."
5. `round30/README.md` §0 item 4 and §4 "For: (R$_4$) has exactly the Vojta shape": add "(uniform bound: proved in
   round 31, Part IV; the Vojta shape concerns only the sharp constant)". In the small-cases table, row $8\mid D$, upper
   bound "$4(D+2)$" $\to$ "$\min(4(D+2),24)$"; the same 24 applies to every row with $D\ge5$.
6. `round30/README.md`, the conditional "unbounded $\mathbb Q(i)(t)$-clusters would refute $A_5(\delta)$" (round 30 §6
   per bridge): mark it "premise false (round 31, Part IV); route dead".
7. Standing cautions (expert note): add F3: *"An argument that bounds integral fakes without the exact norm must use
   characteristic 0 or small-residue-field collisions; the characteristic-2 trace seed satisfies every
   characteristic-free inequality."*

## 5. Open sub-questions, ranked by honest probability of a worst-case growth-rate improvement

All probabilities are the judge's estimates [J]. All are small. Part IV mostly locates the obstruction more precisely.

1. **An arithmetic orbifold BMY for $\mathbb P^1_{\mathbb Z}$ with sections, any constant $\kappa$, $O(W)$ error, in
   the growth regime only.** It would give a uniform bound for integral fakes and hence for clusters. P ≈ 0.003 (it
   is at least Szpiro strength). *First step:* write the exact inequality needed (lit's $K+\partial$ form on
   $\overline M_{0,m}$ at the primes of $N$, with explicit error), and decide whether it implies abc or only
   Szpiro-for-a-family. Then file it next to the Vojta $D=0$ route rather than pursue it.
2. **F3': growth-scale integral fakes over $\mathbb F_p[t]$ with $p$ fixed** (e.g. $p=2$, $\deg N=O(m\log m)$). A
   construction would show that small primes alone cannot give growth over $\mathbb Z$. A non-existence proof would
   be a template for $\mathbb Z$ that F3 does not kill. P ≈ 0.003. *First step:* exhaustive or randomised search over
   $\mathbb F_2[t]$ for configurations with $m\ge2\deg N/\log_2\deg N$, full divisibility and $\delta\le0$, using the
   record's pigeonhole bound as the target.
3. **Orbifold-BMY bookkeeping at the unit-template level $\delta=L/8$ with the template involution**, i.e. the
   mechanism that holds the integer records. P ≈ 0.002. *First step:* compute $3e_{\rm orb}-(K+B)^2$ on the
   function-field model of the golden family (involution quotient surface) and see whether any sub-vacuous bound
   survives.
4. **An anisotropic (circle-type) characteristic-$p$ cluster with growing $m$** over $\mathbb F_{q^2}(t)$. This tests
   whether exact norms, rather than characteristic 0, block constructions. P ≈ 0.001. *First step:* the referee's
   linear solve for norm-compatible forms on $p$-linearised seeds in odd characteristic, $q=27,49,81,125$.
5. **The sharp constant in the function-field model** ($24\to8$, i.e. (R$_4$)), using the isotropic points $[1:\pm i]$
   that marked sections pass through to choose better weights. Any argument that forgets the norm cannot go below 12.
   P ≈ 0.001 for the integer goal; it is worth doing only as model calibration. *First step:* build or exclude a
   16-mark $\delta=0$ no-norm configuration (biplane of order 4 cut out by $(1,1)$-curves on a quadric).
6. **Read Langer PLMS 86 (2003) Theorem 0.1 in the original.** P = 0 for the goal, but required before (IV.2) is
   recorded without a caveat. *First step:* obtain the paper (library, or arXiv math/0012180 once arxiv.org is
   reachable) and check the hypothesis "$K_X+D$ effective" and the local $e_{\rm orb}$ for arbitrary weights.

## 6. Scores (0-10)

| block | score | reason |
|---|---|---|
| bmy ((IV.2) and its proof) | 4/10 | Correct, new, carefully executed transplant of Langer's BMY with a clever normalisation at infinity. Decisive in the model; nil on integers. |
| examples (direct tests, extremal data) | 3/10 | Independent confirmation by direct Chern numbers; new 12- and 9-mark configurations; shows 24 is the method's limit. |
| trace (char-2 seed, symmetry exclusion) | 2/10 | Correct and elementary; useful calibration (F3); the exclusion needs $d\ge15$ and targets a table already excluded for genuine clusters. |
| conseq (specialisation, sampling, 34/35) | 1/10 | Correct, routine; excludes only rational-branch templates, the weakest mechanism. |
| bridge (function-field theorem) | 4/10 | Settles the model's open uniform question, a result Astra did not see; three referees concur; no integer content. |
| lit (citations, context) | 2/10 | Places the method in its lineage and gives the clean $K+\partial$ reading; originals unread. |
| arith (substitutes) | 2/10 | Sharp no-go: every proved Arakelov tool fails F3; any $\kappa$ would do; nothing constructive. |
| **Part IV, on the goal** | **2/10** | The only genuine new theorem of the round, and it closes the function-field model's uniform question (over C). It moves nothing for integer parents: growth rate, uniform bound, Questions A and B and route 5 are unchanged, and its arithmetic analogue is an unproved arithmetic orbifold BMY of Szpiro strength. |

Sources consulted by the judge (WebSearch; arxiv.org unreachable):
[Langer, arXiv math/0012180](https://arxiv.org/abs/math/0012180);
[Pokora, The orbifold Langer-Miyaoka-Yau inequality and Hirzebruch-type inequalities, arXiv 1612.05141](https://arxiv.org/pdf/1612.05141).
