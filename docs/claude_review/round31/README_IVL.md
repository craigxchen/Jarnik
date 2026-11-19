# Round 31, Part IV (Question B thread) -- judge report: the norm defect at width $N^{-3/10}$

Source: Part IV of `astra_lemma_IV/answer.md` (heading "# Part IV", §§IV.1-IV.4 and the twelve-line summary; scripts
`iv_shell_certificate.py`, `iv_quartic_conditions.py`, `iv_cubic_shortrow.py`; question `reply4.txt`, Astra's summary
`final4.txt`). This is the Question B ("lemma") thread. The Part IV judge report of the other thread is
`README_IV.md` and is unrelated. Four referees or analysts each took one block:

| key | report | scope |
|---|---|---|
| quarter | `ref4L_quarter.md` | IV.1 (lattice, prime specialisation), IV.2 (quarter-power family, escape from closed sets) |
| balanced | `ref4L_balanced.md` | IV.3 (balanced family, DFI device, numerical parent, half-power pencil obstruction) |
| dichotomy | `ref4L_dichotomy.md` | IV.4 (width-defect dichotomy) and the correction of suggestion (c) |
| rankone | `ref4L_rankone.md` | analysis of the rank-one residual branch of IV.4 (not a refereeing block: its items are mostly the analyst's own) |

Judge scripts are in `judge_IVL/` (`qstar.py`). The Part I-III judge reports are `README.md`, `README_II.md` and
`README_III.md`. **Every claim labelled [J] is the judge's own and is unrefereed.**

Labels as in Part III. "Proved" means a complete proof checked line by line by a referee, or an exact computation a
referee reproduced with independent code. "Computed" means exact data reproduced independently. "Evidence" means
numerical, not certified. **No referee found a false statement or a gap in any claim Astra made.** The only FALSE
verdict (rankone (b3)) and the two GAP verdicts (rankone (b4), (d)) concern questions the analyst posed: whether
IV.2's lower-order terms can be killed, whether every polynomial pencil stops at the quarter power, and whether
anything follows for genuine tuples. Astra claimed none of them. I rule on each in §1.4.

## 0. Verdict in five lines

1. **Nothing changes for the integer problem.** The growth rate stays $(2/3+\varepsilon)\log R/\log\log R$. There is
   no uniform bound $M(C)$, and Question B (an exponent $>15/4$, i.e. $A_5(\delta)$ with $\delta<1/30$) is exactly
   where it was. Astra says so herself (IV.4 "Evidence", summary item 11).
2. **Everything Astra claimed is correct.** The corrections: the IV.2 cofactors grow like $(\log N)^{O(1)}$, so
   exponent $10/3$ is reached only as a limit; $\theta=1/4$ itself is not attained; and $e(\theta)$ depends on the
   orientation. In the endpoint orientation ($n_i\sim N^{2/5}$), $e(\theta)=10/3$ for every $\theta>1/4$, escaping
   every fixed closed set. In the balanced orientation ($n_i\sim\sqrt N$), $e(\theta)=10/3$ for $\theta\ge3/5$
   (one fixed pencil, no escape).
3. **The dichotomy IV.4 is a correct theorem with no content at defect zero.** Below defect $n^{1/4-30\eta}$
   (endpoint) or $n^{2/5-24\eta}$ (balanced), the norms are exact or $\delta_i=\lambda\ell(V_i)^2$. Genuine tuples
   satisfy it trivially, and no known model violates it.
4. **The rank-one residual is a genuine cluster for another quadratic order** $\mathbb Z[\sqrt k]$ with
   $|k|\le N^{3/10-\varepsilon}$ (analyst's Theorem R1). Such configurations were neither constructed nor excluded
   at small defect. As an analytic fact, the IV.2 models are themselves *exact* clusters, Liouville-sharp in the
   reduced metric, for positive-definite forms of discriminant $\asymp N^{3/10}$ (Theorem R3, [J] rechecked). So exact
   norms for *some* form are worthless; a proof must use the Euclidean form, or at least a bounded discriminant.
5. **Score for Part IV (Question B thread) on the goal: 1/10.** Clean, correct and new calibration of which
   information a Question B proof must use. No progress on Question B.

## 1. Referee-corrected statements per block

### 1.1 quarter: IV.1, IV.2 (referee: no error; three corrections of scope)

**Proved**
* **IV.1 lattice.** With each cut marking its two endpoints, the valuation columns $e_i+e_j$ generate $D_5$, and with
  $\mathbf 1$ they generate $\mathbb Z^5$. Reversing a side replaces a column by $\mathbf1-(e_i+e_j)$, so the index
  is unchanged. Index one in both orientations. The anchored Smith form is $(1,1,1,1)$. Liouville/width formula
  $\sin(\text{angle}_{ij})=(1+o(1))|t_{ij}|/\sqrt{n_in_j/g_{ij}^2}$ when $\delta_i=o(n_i)$. Routine.
* **IV.1 prime specialisation, corrected statement.** *For fixed distinct shifts $b$ (including $0$ if needed), a
  fixed modulus $m$, and $A$ any subset of the primes of positive relative density, there are $P\to\infty$ and
  $q=mr$, $0<|r|\le(\log P)^{L}$, with all $P+qb\in A$* (Tao-Ziegler, *Narrow progressions in the primes*, Maier
  volume 2015, arXiv:1409.1327; Shao, arXiv:1509.04955, for the linear case with explicit $L_k$). Corrections:
  $q\le(\log P)^{O(1)}$, so the diagonalisation is unnecessary; $q$ may be negative, which is harmless; in IV.2, $P$
  need only be $\equiv1\pmod4$ and coprime to the primes of $q$. Bounded $q$ would need the ten-tuple prime
  conjecture. The statement was confirmed from abstracts and search snippets only. [J] My own search found the same
  abstract (with $L$ depending on $k$, the $P_i$ and the density parameter in the arXiv version), which is all
  that is used.
* **IV.2 construction.** For $x=(-27,-270,135,120,40)$ ($s_4=0$, $-s_1s_5=97200^2$, $b_0=3s_5/2$), with
  $J_{ij}=P+qb_{ij}$, $U_i=P+qh_i$, $B=b_0q^2$, $L=lq$, $V_i=(PU_i+Bx_i,\ LU_i)$:
  $[V_i,V_j]=LB(x_i-x_j)J_{ij}$ exactly; the columns are primitive; both shell identities hold (the second in its
  pre-$s_4=0$ form, $e_2(\{b_{ij}\}_j)=h_i^2+3s_5x_i-s_1s_5-s_4x_i^2+s_1s_4x_i$); the residual is exactly
  $\delta_i=s_5x_i(s_2x_i-s_1s_2-s_3)\,Pq^3+c_{4,i}q^4$, with mixed signs. Reproduced symbolically for generic marks
  and [J] again in `judge_IVL/qstar.py`.
* **IV.2 consequences, corrected statement.** *With $q\le(\log P)^{O(1)}$: $n_i\sim N^{2/5}$,
  $|\delta_i|\ll_xPq^3+q^4=N^{1/10}(\log N)^{O(1)}$, every $\delta_i\ne0$ for large $P$, the cofactors are
  $lb_0q^3(x_i-x_j)=(\log N)^{O(1)}$ (unbounded), and $\Delta=N^{-3/10}(\log N)^{O(1)}$. Hence
  $|\delta_i|\le n_i^\theta$ eventually for every fixed $\theta>1/4$ (not for $\theta=1/4$), and
  $\log(1/\Delta)/\log N\to3/10$.* The form $\Delta\asymp N^{-3/10}$ with bounded cofactors holds only under the
  prime ten-tuple conjecture with $q=30$.
* **Escape from closed sets, corrected quantifier.** *For every fixed proper Zariski-closed $Z$ in the ratio torus
  there is a shape $x$ (chosen after $Z$) whose models avoid $Z$ for all large $P$.* Shape density was confirmed
  independently: the parametrisation $(b,c,k,m)$ has Jacobian rank 4, and its composite with the shape map has rank
  3. The leading-term argument is valid: the common scalar $\varepsilon=2ilb_0q^3P^{-3}$ is complex, but the lowest
  Taylor form is homogeneous and $q/P\to0$. This defeats every "$\exists Z\ \forall$ tuples" statement, which is the
  form of $A_5$.

**Computed.** Two new prime instances with $q=30$, $P=1642245131041$ and $P=112908008110501$: ten primes
$\equiv1\pmod4$, primitive, $g_{ij}=J_{ij}$, all minors exact, index one, mixed signs. Pre-asymptotic: there
$\log|\delta|/\log n=0.64$-$0.74$ and $\log\mathrm{diam}/\log N=-0.12$ to $-0.14$, because the defect constant is
$\approx10^{24}$ and the width constant $7.5\cdot10^{21}$. The instances certify the algebra, not the exponents.

**Scope (referee, accepted).** The model points are not points of $T_5$ near $e$ at exponent $10/3$ relative to
their own height. The content of $\bar V_5V_i$ is exactly $1$, so the pair height is $N^{4/5}$, not the profile's
$N^{3/5}$. $J_{ij}$ never divides $V_i\cdot V_j$ (0/10 pairs), and no marked block is isotropic (0/20). In this
orientation $n_i^{1/4}\asymp P$ is one block prime, so a single norm congruence $|V_i|^2\equiv0\pmod{J_e}$ at one
marked block excludes the models.

### 1.2 balanced: IV.3 (referee: no error)

**Proved**
* **Profile and identities.** Complementing the cycle $12,23,34,45,51$ gives $|S_i|=5$, $|S_i\cap S_j|=2$. With the
  table $(b_e)$, $(R_i,d_i,f_i)$, $p_*=P+329q$, $T_i=P^2+d_iqP+f_iq^2$:
  $T_iR_j-T_jR_i=(R_j-R_i)g_{ij}$ identically, and $p_*T_i^2-n_i$ has only $P^2q^3,Pq^4,q^5$ terms (the $P^2q^3$
  coefficients are $15177600,13251249,10476400,11379456,10319616$).
* **Construction.** With $q\equiv0\pmod{2\prod R_i}$, the $J_e$ and $p_*$ split primes, $|\alpha|^2=p_*$, $(\alpha,\beta)$
  unimodular and $s=\alpha\cdot\beta$ centred: $V_i=T_i\alpha+R_i\beta$ is primitive, has fixed signed cofactors
  $R_j-R_i$ (largest $1872$), and satisfies exactly
  $\delta_i=(p_*T_i^2-n_i)+2sT_iR_i+\frac{s^2+1}{p_*}R_i^2$.
* **Corrected statement of the main claim.** *There are arbitrarily large squarefree split parents with ten prime
  blocks in the balanced orientation ($n_i\sim\sqrt N$), index one, fixed nonzero signed cofactors, width
  $\Delta\asymp N^{-3/10}$, and $\max_i|\delta_i|=o(n_i^{3/5})$.* Here $|\delta_i|\le(R_i+o(1))n_i^{3/5}$ uses only
  narrow progressions ($q=o(P^{1/3})$ suffices). The little-$o$ uses Duke-Friedlander-Iwaniec (Ann. Math. 141
  (1995)): for $\varepsilon<1/2$ the primes with a root of $X^2+1$ in $(0,\varepsilon p)$ have relative density
  $\varepsilon$. Tao-Ziegler is applied inside that one subset, and a diagonal sequence is taken.
* **Pencil obstruction, corrected statement.** *For this fixed pencil (fixed minor ratios $R_j-R_i$), no integral
  matrix whose minors are proportional to these up to a factor $P^{o(1)}$ has maximum defect $o(P^{5/2})=o(n_i^{1/2})$.*
  The proof needs only $\nu^2=o(P)$ and $q=o(P^{1/6})$. It concerns this pencil only.

**Computed.** The parent $q=924$, $K=3999318414$, $p_*=K^2+1=P+329q$ (prime), $P=15994547776559171401$ (not prime;
not needed): five blocks prime, all squarefree, pairwise coprime and split; Smith form $(1,1,1,1)$; cofactors
exactly $R_j-R_i$; $0<(1000|\delta_i|)^5<n_i^3$; $\mathrm{diam}\cdot N^{3/10}=3744$. About $99.7\%$ of each defect is
the $q^3P^2$ term. The small $s=-K$ comes from $p_*=K^2+1$, so the instance does not illustrate the DFI mechanism.

**Scope.** This is one projective pencil (cross ratios depend on $P/q$ only), so it does not escape closed sets.
Balanced defect exponents in $[2/5,3/5)$ are undecided.

### 1.3 dichotomy: IV.4 and the correction of (c) (referee: no error)

**Proved.** The clean statement, as verified (`ref4L_dichotomy.md` §4), quoted exactly:

> **Theorem (IV.4).** Let $h\in\{4,5\}$ be the orientation ($h=4$: each pair cut marks its two endpoints; $h=5$:
> the five cycle edges are complemented). Let $J_e$ ($e\in E(K_5)$) be pairwise coprime positive integers with
> $J_e=P^{1+o(1)}$. Put $n_i=\prod_{e\in S_i}J_e$ and $g_{ij}=\prod_{e\in S_i\cap S_j}J_e$. Let
> $V_1,\dots,V_5\in\mathbb Z^2$ be primitive with $g_{ij}\mid m_{ij}=\det(V_i,V_j)$, $t_{ij}=m_{ij}/g_{ij}\ne0$,
> $\max|t_{ij}|\le P^{\kappa+o(1)}$, and $\delta_i=|V_i|^2-n_i$. Fix $\varepsilon>0$. As $P\to\infty$:
> 1. if $\max|\delta_i|\le P^{h-1-6\kappa-\varepsilon}$, then $\sum_ic_iA_iV_i=0$, and there is a unique
>    $Q_*\in\mathbb Q[X,Y]_2$ with $Q_*(V_i)=n_i$ and $H_tQ_*$ integral ($H_t=\prod|t_{ab}|$);
> 2. if $\max|\delta_i|\le P^{h-3-12\kappa-\varepsilon}$, then either all $\delta_i=0$, or $E=I-Q_*=\lambda\ell^2$
>    with $\lambda\in H_t^{-1}\mathbb Z\setminus\{0\}$ and $\ell$ primitive integral. In the latter case
>    $\delta_i=\lambda\ell(V_i)^2$, and every $\delta_i$ is nonzero, with the sign of $\lambda$.
> 3. If the angular width is $\Delta\le N^{-3/10+\eta}$, one may take $\kappa=10\eta$.

In terms of $n_i$, item 2 at width $N^{-3/10+\eta}$ reads: defect below $n_i^{1/4-30\eta}$ ($h=4$) or
$n_i^{2/5-24\eta}$ ($h=5$). Hypotheses actually used: pairwise coprime (private) blocks, integral and primitive
$V_i$, $J_e=P^{1+o(1)}$. Primality is not used. The interpolation identity is Part III's (8) and is not new.

* **Strengthening (referee).** Astra's last remark holds in the whole range, not only at $\kappa=0$: a single zero
  defect forces all defects to vanish, since $|\delta_j|\ge g^2t^2/H_t\ge P^{2(h-3)-10\kappa}>P^{h-3-12\kappa}$.
* **Sharpness, corrected statement.** *For $h=4$ and $\kappa\to0$ the exponent $h-3=1$ cannot be raised (IV.2:
  mixed signs, defect $P^{1+o(1)}$, $\kappa=o(1)$). Neither the $12\kappa$ dependence nor the $h=5$ exponent $2$ is
  known to be sharp.* The best balanced model is at $P^{3-o(1)}$ (IV.3).
* **Correction of the coordinator's suggestion (c).** At a private block prime the marked columns give **one**
  projective line. The lattice $\{Q\in\mathbb Z[X,Y]_2:Q(V_i)\equiv0\bmod n_i\ \forall i\}$ has index exactly $N$
  (proof, plus Smith normal form in 29 tests). The "two lines $L_p^\pm$" premise fails in the star normalisation,
  even for genuine stars (3-4 directions mod $p$); it holds only in the common-norm normalisation. The covolume-$N^2$
  class proposed in (c) is the wrong object.

**Computed.** $\det E=-(4112/243)^2$ for every $P,q$ tested on the IV.2 family. [J] Confirmed symbolically in
$P,q$ (`judge_IVL/qstar.py`): $\det E=-16908544/59049$ identically. So IV.2 is split rank two at bounded
determinant, the equality case of the Gram bound, one step from rank one. Rank-one models with large defect exist
algebraically: stretched stars ($\lambda=1-1/d^2$), $\mathbb Z[\sqrt{-2}]$-stars ($\lambda=-1$), $\mathbb Z[\sqrt2]$-stars.

### 1.4 rankone: the residual branch (analyst; rulings on its FALSE and GAP items)

**Proved (analyst's own, refereed by nobody but checked in its report; I read the proofs)**
* **Theorem R1 (normal form), corrected statement.** *In the rank-one branch, choose $e_1,e_2$ with $\ell(e_1)=0$,
  $\ell(e_2)=1$; put $A=|\ell|^2$, $B=e_1\cdot e_2$, $u_i=Ax_i+By_i$, $k=A\lambda-1$. Then $An_i=u_i^2-ky_i^2$,
  $\operatorname{disc}Q_*=4k$; $k=0$ is impossible on private profiles; $\lambda<1/A$ gives a definite form and
  $\lambda>1/A$ an indefinite one. Below the IV.4 threshold, $|k|\le N^{3/10-\varepsilon'}$. With
  $\mathfrak a=[A,B+\sqrt k]$, $Q_*$ is the norm form of $\mathfrak a$; the $\xi_i=u_i+y_i\sqrt k$ generate ideals
  $\mathfrak b_i$ of norm $n_i$ in one class; block primes not dividing $2k\cdot\mathrm{den}$ are split, never inert;
  at split blocks the Gaussian orientation argument applies. Hence a small-defect fake is a genuine
  $\mathbb Q(\sqrt k)$-cluster, a point of $T_5^K$ with the same profile (minus ramified blocks).* [J] The non-inert
  step uses $|y_i|<p$, which follows from $|\lambda|y_i^2=|\delta_i|$ and $|\lambda|\ge1/H_t$:
  $|y_i|\le(H_t\max|\delta|)^{1/2}\le P^{(h-3)/2+O(\kappa)}<P$. Correct. Non-integral $\lambda$ costs a factor
  dividing $H_t$.
* **Theorem R2.** Liouville $10/3$ for any binary form: $\sin\theta_F\ge(\sqrt{|D|}/2)\,g/\sqrt{nn'}$. Routine.
* **Theorem R3.** *Every IV.2 model is an exact cluster for its interpolating form $Q_*$, which is positive definite
  with $\operatorname{disc}Q_*$ a cubic in $P/q$ of size $\asymp P^3\asymp N^{3/10}$, fixed denominator for fixed $q$,
  no inert block primes.* [J] Reproduced: the interpolation from three columns is exact on all five, and
  $\operatorname{disc}Q_*=(-719P^3+3704109750P^2q-\dots)/(334731302496000000\,q^3)$.
  *Corrected meaning* [J]: the points are Liouville-sharp only in the reduced metric $\mathrm{dist}/\sqrt{|D|}$. In
  absolute terms they are at distance $\asymp N^{-3/20}$, far outside Liouville. What R3 refutes is any version of
  $A_5$ uniform over imaginary quadratic orders with a constant growing like $|D|^a$, $a>7/18$ [J: from
  $|D|^aN^{-4/15}>N^{-3/20}$ at $|D|=N^{3/10}$]. Nobody proposed such a version. Its value is as a barrier: exact
  norms of *some* form do not help unless the form's height or discriminant enters.

**Ruling on (b3) FALSE.** This is the analyst's own question ("can IV.2's leading residual be made rank one and the
lower terms killed?"), not a claim of Astra. The exact reduction is right: rank one iff $s_3=-s_1s_2$, and then
$\delta_i=(q^3s_5/4)(4Ps_2-4qs_1s_2^2+5qs_5)x_i^2$ exactly; going below $P$ needs $s_2=0$, impossible because
$X^5-s_1X^4-s_5$ has at most three real roots. **Upheld as a closed route**, exact. (The $s_2\ne0$ rank-one case
still needs integral points on a degree-7 surface; none with $|y|\le30$. That is not a proof of emptiness.)

**Ruling on (b4) GAP.** Again the analyst's own conjecture: quarter-power optimality among one-parameter
polynomial pencils ($h=4$, linear blocks, constant cofactors). Proved for constant $\ell$ (inert-prime argument;
exact Gröbner basis in the split case) and for $\deg\ell\ge2$. The general case $\deg\ell=1$ is reduced (validated
to 40 digits) to the $\sigma$-system, whose found solutions are all degenerate (720 real, about 200 complex per
$\beta$, 28 by a direct solve). **This is evidence, not a theorem; the gap is repairable** by a certified solve
(msolve or an interval homotopy), which was not available here. Record it as "Evidence".

**Ruling on (d) GAP.** Nothing new for genuine tuples, for the structural reason that IV.4 at $E=0$ is vacuous and
R1-R3 concern other forms. This is not a gap in a claim. It is a correct negative assessment, which I accept.

**Literature (c): not checked.** Only snippets: Ramana (Acta Arith. 130 (2007)) gives at most 2 points on arcs of
length $R^{1/6}/|d|^{1/3}$ on $X^2+dY^2=R$. No uniform-in-the-form $k\ge3$ statement was found. The attribution of
a two-point lemma for $x^2-Dy^2=n$ to Cilleruelo-Jiménez-Urroz is unverified.

### 1.5 Withdrawn or corrected in Part IV (collected)

1. IV.1: $q\le(\log P)^{O(1)}$ (not just $P^{o(1)}$); $q$ may be negative; $P$ need not be prime in IV.2.
2. IV.2: cofactors $(\log N)^{O(1)}$, not bounded; the exponent is a limit; $\theta=1/4$ is not attained; the
   statement belongs to the endpoint orientation. "Bounded cofactors" with $\Delta\asymp N^{-3/10}$ is conditional on
   the ten-tuple conjecture.
3. Summary item 9 / final4 item 9, "the quarter-power threshold is sharp for mixed-sign defects": sharp only in the
   exponent $h-3=1$ at $\kappa\to0$, for $h=4$.
4. IV.4 last remark: holds for every $\kappa$ in range, not only $\kappa=0$.
5. IV.3: the numerical parent does not illustrate DFI; the obstruction concerns one pencil.
6. Numerical instances of both families are pre-asymptotic and certify algebra only.

## 2. What is new relative to the record, and correct

Grep: the record has the III.11 model (absolute defect $n^{4/5}$, relative $N^{-1/10}$), the C1 family and the
gap-only angle model at $10/3$. It has no "norm defect" analysis outside round 31, no rank-one dichotomy, and nothing
on quadratic forms of growing discriminant.

* **IV.2 quarter-power family** (new): integral five-mark fakes with defect $N^{1/10+o(1)}$ in the endpoint
  orientation, relative defect $N^{-3/10+o(1)}$ (the order of the width), escaping every fixed closed set. It lowers
  III.11's absolute defect by a factor $N^{3/10}$.
* **IV.3 balanced family** (new): $o(n_i^{3/5})$ with fixed cofactors at $n_i\sim\sqrt N$, via a DFI centring
  device; plus a half-power obstruction for that pencil.
* **IV.4 dichotomy** (new): mixed-sign fakes cannot go below $n^{1/4-30\eta}$ (endpoint) or $n^{2/5-24\eta}$
  (balanced); below that, exact or rank one.
* **Correction of (c)** (new, corrective): one marked line per prime; index $N$.
* **R1/R3** (analyst, new): the residual is a genuine $\mathbb Q(\sqrt k)$-cluster with $|k|\le N^{3/10}$; IV.2 models
  are exact clusters for forms of discriminant $\asymp N^{3/10}$. **(b3)** closed exactly; **(b4)** evidence.
* **[J] $\det E\equiv-(4112/243)^2$** on the IV.2 family (referee observation, symbolically confirmed).

Not new: the interpolation identity (Part III (8)); Liouville for forms (R2); the Gram device (II.9/III).

## 3. Rulings

### 3(a) The function $e(\theta)$

Definition used: $e(\theta)$ is the best (smallest) approximation exponent attained by integral five-mark models
(squarefree split parent, ten private pair blocks, index one, all ten divisibilities, cofactors $N^{o(1)}$,
$|\delta_i|\le n_i^\theta$), as a limit $\log(1/\Delta)/\log N\to1/e$. Liouville gives $e(\theta)\ge10/3$ whenever
$\theta<1$. **$e(\theta)=10/3$ means that a method using only these data at defect level $\theta$ cannot beat
Liouville.** The value depends on which side of each cut is marked, because the defect is absolute.

* **Endpoint orientation ($n_i=\prod_{e\ni i}J_e\sim N^{2/5}$, $h=4$).** Proved: $e(\theta)=10/3$ for every
  $\theta>1/4$ (IV.2), with cofactors $(\log N)^{O(1)}$, and also after removing any fixed proper closed set. Open:
  $\theta\le1/4$. Proved constraint (IV.4): a model at width $N^{-3/10+\eta}$ with defect $\le n_i^{1/4-30\eta-\varepsilon}$
  is exact (hence genuine) or rank one of one sign. Evidence: no one-parameter polynomial pencil does better than
  the quarter power (rankone (b4)).
* **Balanced orientation ($n_i\sim\sqrt N$, $h=5$, cycle complemented).** Proved: $e(\theta)=10/3$ for every
  $\theta\ge3/5$, with fixed cofactors and $\Delta\asymp N^{-3/10}$, one pencil (no closed-set escape). Earlier
  (III.11): $\theta\ge4/5$ with bounded cofactors. Open: $[2/5,3/5)$; that pencil cannot reach below $1/2$. Proved
  constraint: below $n_i^{2/5-24\eta}$, exact or rank one.
* **$\theta=0$** in either orientation is the genuine problem (exact norms with all divisibilities force a genuine
  cluster, III.11). Nothing in Part IV bears on $e(0)$.

### 3(b) The dichotomy, as a theorem and its meaning

The theorem is §1.3, quoted. Its meaning:
1. **Is it the first statement that genuine five-tuples satisfy and integral models of small defect violate?** No,
   in any useful sense. Genuine tuples satisfy it trivially (the exact branch). No known integral model violates it:
   IV.2 ($P^{1+o(1)}$), IV.3 ($o(P^3)$) and III.11 ($P^4$) all lie above the threshold. It forbids only hypothetical
   mixed-sign models below the threshold, and none were ever proposed. The record already has a sharper separation
   at defect zero: exact norms plus divisibilities force genuineness (III.11).
2. **Why nothing at defect zero.** At $\delta=0$, $E=0$ and the theorem's conclusion is its hypothesis. The
   argument measures how far the Euclidean form is from the interpolating form. It uses neither the width (beyond
   $\kappa$) nor the closed set $Z$, so it can never yield an exponent for genuine tuples.
3. **What it does say.** The IV.2 barrier cannot be pushed lower by mixed-sign constructions. Any method that
   tolerates defects below the threshold need only handle exact norms and rank-one defects, and the latter are
   exact clusters for another form (3(c)).

### 3(c) The rank-one residual case

**Yes, it is a genuine cluster for another quadratic form** (Theorem R1): a $\mathbb Q(\sqrt k)$-cluster for the order
$\mathbb Z[\sqrt k]$ (possibly non-maximal, possibly with class group, possibly real quadratic), with $|k|\le
N^{3/10-\varepsilon'}$, lying in the strip $|\ell(V)|\le P^{(h-3)/2+O(\kappa)}$ and within rank one of the Euclidean
form. **The analyst neither constructed nor excluded such configurations at small defect.** What was done:
* rank-one models were exhibited at *large* defect (stretched, $\mathbb Z[\sqrt{\pm2}]$-stars);
* the IV.2 route to rank one was closed exactly ((b3));
* polynomial pencils were excluded except for one numerical step ((b4));
* exact clusters at reduced exponent $10/3$ for forms of discriminant $\asymp N^{3/10}$ were shown to exist (R3).
  Those are rank-*two* in $E$, so they are not in the residual branch.

The residual branch is therefore the "middle discriminant" regime $N^{o(1)}<|k|\le N^{3/10}$. Fixed-$K$ Question B
does not address it. With bounded $k$ it would reduce to Question B for $\mathbb Q(\sqrt k)$, a different field from
$\mathbb Q(i)$, and that is not implied by Question B as stated.

### 3(d) What changes (blunt)

* **Growth rate:** unchanged, $(2/3+\varepsilon)\log R/\log\log R$.
* **Uniform bound $M(C)$:** unchanged.
* **Question B ($A_5(\delta)$, $\delta<1/30$; exponent $>15/4$):** unchanged. Part IV proves only negative
  calibrations:
  * archimedean norm accuracy $N^{1/10+\varepsilon}$, plus all divisibilities and index one, cannot beat $10/3$ even
    off any fixed $Z$;
  * exact values of *some* rational form of small denominator cannot either (R3);
  * so a proof must use **local isotropy** (dot-product or Gaussian-gcd divisibility: $J_{ij}\mid V_i\cdot V_j$),
    or the fact that the form is $X^2+Y^2$ (or of bounded discriminant).
  None of this is a lever. It says which information the missing lever must consume, which the record essentially
  knew (III.11: "only exact norms are missing").
* **Route closures to record:** "kill IV.2's lower-order terms to reach rank one" (closed exactly); mixed-sign
  small-defect fakes (closed by IV.4); the coordinator's quadratic-form lattice (c) (wrong object).

## 4. Recommended updates to the expert note (recommendations only; the record is read-only)

1. §5 item 12, after the "Integrality with all divisibilities is not enough" bullet, add:
   > *Norm defect (round 31, Part IV, Question B thread).* Write $\delta_i=|V_i|^2-n_i$. In the endpoint
   > orientation ($n_i\sim N^{2/5}$) there are integral five-mark models with all divisibilities, index one,
   > cofactors $(\log N)^{O(1)}$, $|\delta_i|=N^{1/10+o(1)}$ (every $\theta>1/4$) and width $N^{-3/10+o(1)}$, avoiding
   > any prescribed proper closed set of ratios. In the balanced orientation ($n_i\sim\sqrt N$) one pencil gives
   > $o(n_i^{3/5})$ with fixed cofactors. The models are non-isotropic at every block prime. Below
   > $n^{1/4-30\eta}$ (endpoint) or $n^{2/5-24\eta}$ (balanced), any model is exact or has
   > $\delta_i=\lambda\ell(V_i)^2$; the latter is a genuine cluster for $\mathbb Z[\sqrt k]$, $|k|\le N^{3/10}$. Nothing
   > here bears on defect zero.
2. §3 (Question B), where admissible inputs are discussed: add "a proof must use local isotropy or the Euclidean
   form itself; exact values of a form of discriminant $\asymp N^{3/10}$ do not suffice (the IV.2 models are exact
   clusters for such forms, Liouville-sharp in the reduced metric)".
3. Closed routes: add "killing the IV.2 lower-order terms (exact; Rolle)" and "the quadratic-form congruence lattice
   of suggestion (c): index $N$, one marked line per prime".
4. Wherever $e(\theta)$ is mentioned: always state the orientation and the normalisation of $n_i$.

## 5. Open sub-questions, ranked by honest probability of progress on Question B

All probabilities are the judge's estimates [J].

1. **Isotropic integral models.** Do integral five-mark models exist with $J_e\mid V_i\cdot V_j$ for a positive
   proportion of blocks (or all blocks) at width $N^{-3/10+o(1)}$ off a fixed $Z$? An existence result would push
   the barrier to the Gaussian-gcd level; non-existence would identify the input a Question B proof can use.
   P ≈ 0.02. *First step:* impose $J_e\mid V_i\cdot V_j$ at the five cycle edges on the IV.2 pencil, as extra
   polynomial congruences in $(P,q)$, and count dimensions (with exact checks by sympy).
2. **The middle-discriminant regime.** Prove or refute an $A_5$-type statement for $T_5^K$ with explicit
   polynomial dependence on $|D_K|$, uniform for $|D_K|\le N^{3/10}$, strong enough to kill the rank-one branch. P ≈
   0.01 (this is at least as hard as Question B). *First step:* compute the reduced-metric exponent of the record's
   C1 and golden families transferred to $\mathbb Z[\sqrt{-2}]$, to see whether $10/3$ is attained there off fixed
   $Z$.
3. **Balanced window $[2/5,3/5)$.** Vary the minor ratios of the IV.3 pencil to get below the half power, or prove a
   rank-one-style obstruction in the window. P ≈ 0.003 for Question B (peripheral either way). *First step:* rerun
   the cubic-pencil interpolation with the $R_i$ as unknowns and impose vanishing of the $P^2q^3$ coefficient.
4. **Certify the $\sigma$-system.** Turn "quarter power optimal among polynomial pencils" into a theorem (msolve or
   interval homotopy). P ≈ 0.002 for Question B; P ≈ 0.6 of success on its own terms. *First step:* run a certified
   Gröbner computation over $\mathbb F_p$ for several large $p$, then lift.

## 6. Scores (0-10, on the goal)

| block | score | reason |
|---|---|---|
| quarter (IV.1-IV.2) | 2 | correct, new, escapes every $Z$; calibrates the barrier at the archimedean level only |
| balanced (IV.3) | 1 | correct and neat (DFI device); one pencil, no escape, superseded in strength by IV.2 |
| dichotomy (IV.4, (c)) | 2 | a clean correct theorem and a needed correction of (c); vacuous at defect zero |
| rankone (analysis) | 2 | R1/R3 identify that the form's discriminant must enter; (b3) closed; (b4) evidence |
| **overall Part IV** | **1** | correct and well calibrated, no movement on Question B, the growth rate or the uniform bound |

Sources consulted by the judge: [arXiv:1409.1327 (Tao-Ziegler)](https://www.arxiv.org/abs/1409.1327),
[arXiv:1509.04955 (Shao)](https://ar5iv.arxiv.org/html/1509.04955) (abstract and snippet level only).
