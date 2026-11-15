# Round 29 — adversarial review of lane ff (function fields)

## 0. Verdict

**Survives, as an obstruction report, with corrections. Novelty: mostly a repackaging of known barriers; one
genuinely new exact object.**

1. **Theorem C: confirmed exactly** (own sympy code). It is the record's first explicit 8-point function-field cluster
   and the first 8-point polynomial template. It also corrects `families.md`'s "≤6 for even $D$" impression.
2. **Theorem B: proof gap.** Sign normalisation can make columns coincide. At $L=k$, including ff's own $D=8$, $c=0$,
   the construction is invalid. It is repaired for $L\ge k+1$ by putting a triple $\{a,b,a+b\}$ into $B$. The F1
   conclusion survives.
3. **Theorem D: correct after two repairs** (complexification; the full-rank case). It is **the record's** free-curve
   argument (`gauss.md` Thm 4.3(ii), `toric.md` Prop 4.3), transplanted. ff missed that citation.
4. **The linearisation is in the record** (`families.md` Lemma P). The "rank bounded iff $M$ bounded" equivalence is
   the record's item 8 in function-field dress, and for non-squarefree $N$ it needs a packing argument instead of
   Odlyzko.
5. **Table corrections.** $M(7,0)\ge4$ is exact (record). New exact bounds $M(7,1)\ge8$ and $M(8,1)\ge14$ (ff had 6 and
   8), via roots on a vertical line.
6. **Confirmed independently:** Prop 3.3's LP values (6, 8; 4, 5) and $M(5,0)=3$, $M(6,0)=5$.
7. **Not refuted.** (R$_4$) and the prediction $\limsup M=8$ were not refuted: compositions of the cube at $D=16$ stay at
   rank 4. They remain a heuristic plus an open question.
8. **Assessments without proof.** The SMT / log-Wronskian "exactly critical" claims are unproved sketches.
9. **No growth improvement.** The arithmetic transfer is negative, consistently with filter F1.

## 1. Re-derivation of the linearisation and pair bound (Lemmas 1.1, 1.2, Thm A)

**Lemma 1.1: correct, but not new.** Monic $z_a,z_b$ of degree $D$ have $\deg(z_a-z_b)\le K$ iff $e_1,\dots,e_{D-1-K}=e_L$
agree. In characteristic 0, or $p>L$, Newton turns this into agreement of $p_1,\dots,p_L$, and
$p_m(z_a)-p_m(z_b)=2i\sum_k(a_k-b_k)\operatorname{Im}\rho_k^m$. Re-derived; correct. **Novelty correction:** ff §2 says the
linearisation is "not in the record". It is: `growth/families.md` §4 **Lemma P** is the same statement, proved the same way
(Newton: bounded $C$ iff $\operatorname{Im}p_r(S_j)$ agree for $r<D/2$). The record also uses the coupled-PTE form
(`round22/toric.md` §7: "four coupled complex ideal Prouhet–Tarry–Escott problems"). New in ff is the base-point
**linear-section** form $V'\cap\{0,1\}^D$, together with the rank bookkeeping.

**Lemma 1.2: correct (char-free).** Write $z_a=gA$, $z_b=gB$ with $\deg A=\deg B=w=\|a-b\|_1$. Then $\deg(A-B)\le w-L-1$, so
$w\le L$ forces $A=B=1$. In PTE language this is the classical fact that size exceeds degree.

**Theorem A: correct.** The unary encoding is an isometry from $\ell^1$ on $\prod[0,e_k]$ to Hamming distance. Shortening
gives $A_2(D,d)\le2^{D-2d}A_2(2d,d)$ with $D-2d=2K-D\le2c$. Plotkin gives $A_2(2d,d)\le4d$; ff writes $4(d+1)$, which is
weaker and still true. For $D<2d$, Plotkin gives $\le2d$. The final constant is generous but true. (N1) needs
$D>2c+1$; otherwise $L\le0$ and the bound is trivial.

**Minor gap in "sup $M<\infty$ iff the rank is bounded".** ff cites Odlyzko. That is valid only for squarefree $N$, where
the points are 0/1 vectors. The $\{0,1\}$ form of the lemma is confirmed by WebSearch (Costello–Tao–Vu, arXiv
math/0505156, Lemma 3.1). For general $N$ the points lie in a box, and an $r$-dimensional affine section of a box can
contain more than $2^r$ points. **Repair (proved here):** cluster points lie in an $r$-dimensional affine space, have
pairwise $\ell^1$ distance $\ge L+1$ and $\ell^1$ diameter $\le D$. A volume packing argument in the induced norm then
gives $M\le(1+2D/(L+1))^r\le(5+o(1))^r$ for fixed $c$. The converse $M\ge r+1$ is trivial. So the "iff" holds, and it is
elementary. It is also the function-field form of the record's "bounded rank is equivalent to the uniform theorem"
(`expert_note.md` §5 item 8). So it is a repackaging, not a new reduction.

## 2. Theorem B (fakes): a gap in the proof, false as proved at the boundary, repaired

The orthogonality part is correct. $\ker W\supseteq\operatorname{span}\{\mathbf 1,\chi_a:a\in A\}$, which gives $2|A|+2=D+2c+2$
exhibited 0/1 points at pairwise distance $D/2$ or $D$. This was checked exactly in `kill_ff_checks/thmB.py`.

**The gap.** ff argues that the columns are distinct "because $B$ spans". That gives distinctness **before** the sign
normalisation. The normalisation flips every column $u$ with $g_u<0$. Two columns with $w_u=-w_{u'}$ therefore become
**equal** after it, and equal columns mean a repeated root. Now $w_{u'}=-w_u$ iff $b\cdot(u+u')=1$ for all $b\in B$, i.e. iff $B$
lies in an "odd" affine hyperplane $\{b:b\cdot x=1\}$.
* **At $L=k$**, i.e. $c=D/2-k-1$, the top of ff's stated range, $B$ is forced to be a basis. Then $x=$ (the dual-basis sum)
  always works, and every column pairs off with its negative. The normalised fake has only $D/2$ distinct columns, so
  it is **not a fake**. This covers in particular **$D=8$, $c=0$**, the case ff's table uses. `thmB2.py` checks it:
  `(8,3,0,'basis',obstructions=1,distinct=False)`; it also fails at $(16,L=4)$ and $(32,L=5)$.
  For the basis choice, $\ker W\cap\{0,1\}^D$ also contains points at distance 2 (`thmB.py`: min distance 2 at $D=8$).
* **Repair (proved).** If $L\ge k+1$, choose $B$ spanning **and containing a triple $\{a,b,a+b\}$**. The values
  $b\cdot x$ on a triple sum to 0, so they are never all 1, and no odd hyperplane contains $B$. The normalised columns
  are then distinct and none is a negative of another. Checked numerically for $k=4$, $c=0,1,2$ and $k=5$, $c=0,5$
  (`thmB2.py`, all `True`).
* **Corrected Theorem B.** For $D=2^k$ and $0\le c\le D/2-k-2$ there are fakes with $M\ge D+2c+2$. At fixed $c$ this
  costs nothing asymptotically. So the qualitative (F1)-transplant conclusion, that linear structure plus pair
  inequalities allow linear growth, **survives**. The case $D=8$, $c=0$, where a fake would have 10 points, is
  **unproved**.

## 3. Theorem C (8-point cube family): confirmed exactly

`kill_ff_checks/thmC.py` (own code, sympy) checks the following.
* The identity holds: $abcd\sum(3/y-y^3)+3(a+b)(a+c)(b+c)(1-abcd)\equiv0$ for $d=-(a+b+c)$; the numerator factors as ff states.
* For $y=(-\frac94,-\frac23,\frac14,\frac83)$ there are exactly 8 of the 256 transversals within degree $\le4$ of $Q$, and
  they are the 8 listed sets. The maximum pairwise difference degree is 4, and $N=Q\bar Q$ is real and squarefree.
  Every cluster containing $Q$ lies inside these 8, so the maximal cluster through $Q$ has exactly 8 points.
* The proof re-derives: $\operatorname{Im}r^2=2xy=2\varepsilon$ and $\operatorname{Im}r^3=3x^2y-y^3=3/y-y^3$. The algebra holds in any
  characteristic $\ne2,3$, so the $\mathbb F_q$ claims follow from the identity. Composition multiplies degrees
  exactly, so the extension to $D=8q'$ needs no Newton identities and is valid over $\mathbb F_q$ too.
* Corollary C′ (`template.py`): scaling by 72 gives roots in $\mathbb Z[i]$. At $u=10^6$ and $u=10^9+7$ the 8 values are
  distinct Gaussian integers of one norm, with chord$/\sqrt R=3.1796\cdot10^9$. Confirmed.
* **Composition test** (the lead ff proposes): $Q_{\rm cube}(h(t))$ for $h=t^2$, $t^2+t$ and $t^2+\frac37t-\frac25$ at $D=16$,
  $c=0$, $L=7$. All $2^{16}$ subsets were enumerated (double precision, re-checked to 40 digits). There are exactly 8
  zero-sum sets, of rank 4 (`compose16.py`). Compositions add no rank, so they do not contradict (R$_4$).
* **Record context.** This is the first exact 8-point function-field cluster in the record. It also corrects the
  record's remark that polynomial templates top out at 6 points for even $D$ (`families.md` §4, line 284; that
  remark was a box-search observation, and Theorem C′ shows the box was too small). At bounded $C$ it is not a record
  ($M_\infty(C)\ge8$ for $C>5.898$ is known). It has $c=0$, i.e. exponent exactly $\frac12$, so it gives no A-points and
  says nothing about $(L_5^\circ)$ or Question B.

## 4. Small cases (Props 3.1–3.3, table): independent recomputation

**Prop 3.1: correct.** For $L=1$, a fibre of $\sum a_ky_k$ with $y_k>0$ is an antichain. The level-surjection argument is
valid (the theorem of de Bruijn–Tengbergen–Kruyswijk is a memory-quote, standard).

**Prop 3.2: confirmed by a different method.** I reformulated $L=2$ myself. In base-point form the columns
$w_k=v(r_k)=(y_k,2x_ky_k)$ are any vectors with $y_k\neq0$, subject to $w_k\neq\pm w_l$. So
$M=\max|\Sigma\cap\{0,1\}^D|$ over subspaces $\Sigma$ spanned by 0/1 vectors with $e_k\notin\Sigma$ and $e_k\pm e_l\notin\Sigma$:
a generic $W$ in $(\Sigma)^{\perp}$ realises it, and the span of any cluster is such a $\Sigma$. The case $L=1$ drops the
$e_k\pm e_l$ condition. `kill_ff_checks/L2closed.py` does an exhaustive breadth-first search over closed sets, with
floating SVD rank tests on 0/1 matrices of size $\le7$. Results:
* $L=1$: $6,10,20$ at $D=4,5,6$ (Sperner).
* $L=2$: **$M(5,0)=3$, $M(6,0)=5$**, agreeing with ff.

**Prop 3.3: confirmed.** `lp_rank.py` uses its own LP over all $S\subseteq\{0,1\}^r$ (HiGHS). At $\theta=\frac12$:
rank 3 gives at most 6, rank 4 at most 8 (the best $|S|=9$ has optimum $0.4$). At $\theta=\frac47$: rank 3 gives 4,
rank 4 gives 5. As $\theta_D=(L+1)/D\uparrow\frac12$ at fixed $c$, finiteness of the set of $S$ gives a gap below $\frac12$, so the
$\theta=\frac12$ value applies for $D\ge D_1(c)$. The LP is a relaxation (real weights), so the bounds are valid.

**Table corrections.**
* $M(7,0)\ge4$ is **exact**, from the record (`round22/toric.md` Thm 4.4 / `gauss.md` Thm 4.1, $D=7$, contact 4). ff
  labelled it "numerical".
* **$M(7,1)\ge8$** (ff: $\ge6$) and **$M(8,1)\ge14$** (ff: $\ge8$), both **exact** (`L2find.py`, `rank1L2.py`,
  `verify_rank1.py`; sympy over $\mathbb Q(i)$). In each case the cluster through $Q$ has exactly that many points, the
  pairwise difference degree is $\le K$, and $N$ is squarefree. The mechanism ff missed: put all roots **on one vertical
  line** $x=x_0$. Then $\operatorname{Im}\rho^2=2x_0\operatorname{Im}\rho$, the $m=2$ row is proportional to the $m=1$ row, and the cluster is
  the family of zero-sum subsets of signed reals $y_k$ with distinct $|y_k|$:
  - $y=(1,2,3,4,5,-6,-7)$ gives 8 points;
  - $y=(1,2,3,4,-5,-6,-7,8)$ gives 14 points.

  Upper bounds are $A(7,3)=16$ and $A(8,3)=20$. The exhaustive $L=2$, $D=7$ run (`L2_7.txt`) had reached 8 at the time
  of writing; see §10.
* This trick does not scale. For $L\ge3$ the vertical line keeps $\lceil L/2\rceil$ conditions on $D-1$ parameters, which is
  the same ratio as generic, so there is no excess. A purely imaginary integer search at $c=0$ found only 3 points
  ($D=8$, $|y|\le18$; $D=10$, $|y|\le15$; `vert.py`). It affects only $L\le2$, i.e. $D\le2c+6$, so it is irrelevant to
  $\sup_D$ at fixed $c$.

## 5. Theorem D (dominating families): correct after two small repairs; not new in method

**Re-derivation.**
* The cocharacter at $\bar\rho_k$ is the $k$-th coordinate functional restricted to $\Lambda$. It is onto $\mathbb Z$ because
  $\Lambda$ is spanned by vectors with entries in $\{0,\pm1\}$, so it is primitive whatever $\mathbb Z$-basis is chosen.
* Each boundary point has a ray-generator cocharacter. So $f$ meets the open orbit of one $D_\rho$ transversally, and
  $-K_X\cdot f=2\#\{k:\text{column }k\neq0\}\le2D$.
* Contact: $z_F/Q=1+O(s^{L+1})$ with $s=1/t$, so $m=L+1$. In log coordinates $df=O(s^{m-1})ds$. That makes $df$ a
  section of $f^*T_X\otimes\Omega_{\mathbb P^1}(-(m-1)\infty)=f^*T_X(-m-1)=F(-1)$, as ff states.
* Generic smoothness: the image $G$ of $H^0(F)\otimes\mathcal O$ has rank $\ge r-1$ and is globally generated, so $\deg G\ge0$.
  - If $\operatorname{rk}G=r-1$: $df\notin G$ at a general point, so $df$ gives a nonzero section of $(F/G)_{\rm tf}(-1)$.
    Hence $\deg F\ge\deg G+1\ge1$.
  - If $\operatorname{rk}G=r$: ff's "some $c_j\ge1$" is **not justified as written**, since all $c_j$ may be 0 with $G\subsetneq F$.
    The correct argument: either $F\ne G$, and then $\deg F>\deg G\ge0$; or $F=G=\mathcal O^r$, and then $F(-1)$ has no
    nonzero section, which contradicts $df\ne0$.
  - Either way $\deg F=-K_X\cdot f-rm\ge1$, so $r(L+1)+1\le2D$. **Correct.**
* **Repair 2: complexification.** $B$ is a *real* variety, since the conditions involve $\rho$ and $\bar\rho$, and generic
  smoothness is a statement about complex varieties. Fix: treat $\rho_k,\sigma_k(=\bar\rho_k)$ as independent and pass to
  the complex Zariski closure $B_{\mathbb C}$. The family $f_b$, $b\in B_{\mathbb C}$, still has $\le2D$ boundary points with the
  same cocharacters and contact $m$, and a Zariski-dense real image makes $B_{\mathbb C}\times\mathbb P^1\to X$ dominant.
  The proof then runs on $B_{\mathbb C}$, and the inequality is combinatorial. ff omits this; it is easy.
* Arithmetic of the corollaries, checked:
  - $c=0$: $r\le3$ for every $D$, even or odd. For odd $D$, $m=(D+1)/2$ and $r\le(4D-2)/(D+1)<4$.
  - $c\ge1$: $r=5$ needs $D\le10c-2$ for even $D$ and $D\le10c-7$ for odd $D$, so $D\ge10c-1$ gives $r\le4$.
  - Then Prop. 3.3 gives $M\le6$.
* It equals the naive dimension count: image dimension $\le(2D-2-rL)+1\ge r$.

**Not new in method.** This is the record's argument: `round20/gauss.md` Thm 4.3(ii) ("tangent space to maps with this
fixed contact ... generically globally generated, hence nef: $(\#\text{fixed})\mu_0\le-K\cdot C$"), re-derived in
`kill_gauss.md` §4 and restated as `round22/toric.md` Prop 4.3, "covering rational families", $(B\cdot C)/\mu_0\ge k-1$.
ff's §2 cites only `round18/theory.md` §2.3. Theorem D transplants the same free-curve inequality from
$\overline M_{0,k+2}$/$Y_k$ to a toric compactification of $T_\Lambda$ and parametrises by code rank. The "4 false / 5 true for
covering families" pattern is already the record's: `gauss.md` Thm 4.1 at $k=4$ (a covering family, $3\cdot4\le14$) and
Thm 4.3 at $k\ge5$. **That family is itself a function-field cluster:** 14 rays give $D=7$, contact 4, $c=0$, rank 3,
exact (`round22/toric.md` Thm 4.4, `k4class.py`). Hence:

* **Correction to ff's table:** $M(7,0)\ge4$ is **exact** (record), not "numerical".
* `round22/toric.md` §7 already ran the rank-4 dimension count, in the case of the $D=15$ "balanced" 5-point simplex code
  (expected dimension $30-2-28=0$, none found numerically). ff's dimension count agrees with it, and that case is
  consistent with ff's own count.

**F2 for Theorem D.** It allows dominating rank-3 (4-point) families, as the record's $k=4$ family shows. So it does not
prove the false 4-point analogue, and it passes F2. But it says nothing about single configurations. Being a family
statement is exactly why it passes.

## 6. (R$_4$), the "iff", and the tool assessments

* **(R$_4$) is correctly stated as a sufficient condition.** It forbids any 5 independent 0/1 kernel vectors of weight
  $\ge L+1$, which is stronger than needed. Together with Prop. 3.3 it implies $M\le8$ for $D\ge D_1(c)$. The prediction
  "$\limsup M=8$" is a **heuristic**: it rests on the generic dimension count, which is exactly what (R$_4$) would make
  rigorous. It is not a result. My tests did not contradict it: three quadratic compositions of the cube at $D=16$
  give rank 4, and two vertical-line / purely imaginary loci at $c=0$ show no excess.
* **The "iff" is elementary**, with the box repair of §1. It is the function-field twin of the record's "bounded rank ⇔
  uniform theorem" (`expert_note.md` §5 item 8). ff presents it as a new reduction; it is a repackaging.
* **The single-relation tools.**
  - Mason–Stothers on a reduced pair $A-B=c$ with $\deg A=\deg B=w$ and $\deg c\le w-L-1$: for squarefree $N$ the radical has
    degree $\ge2w$, so Mason gives nothing. Confirmed: vacuous, the twin of item 6.
  - The $k$-term log-Wronskian bound $m\le(k+1)(|S|-1)/2$ and the claim that the truncated SMT error is $\asymp|S|$ are
    stated **without proof** in ff.md, so I treat them as unverified sketches. ff's notes contain a sharper
    heuristic: a truncated function-field Vojta of the form $(r-1)m\le|S|-2$ would give $r\le4$, i.e. (R$_4$) itself.
    That is a conditional route; it does not show that every SMT route fails. The write-up's "exactly critical"
    verdict for SMT is therefore **plausible, not proved**.

## 7. Filters F1 and F2 (applied myself)

**F1.** The lane yields no input for the arithmetic problem, so F1 applies only to the function-field model. Inside
the model:
* The corrected Theorem B shows that fakes, meaning arbitrary columns in $\mathbb R^L$, satisfy the linear structure and
  the pair bound with $M\ge D+2c+2$.
* Theorem D does not touch fakes. It concerns families, and a fake is a single configuration.

So **F1 is not passed, and nothing is claimed that would need to pass it.** ff states this correctly. In the
*arithmetic* problem every ingredient ff names as a possible transfer (Fermat quotients, a $p$-adic congruence, the
pair bound) is a point-level residue or a pair inequality. Those are inputs the record's fakes satisfy (item 1). So the
negative transfer verdict is consistent with F1.

**F2.** The single-configuration tools are vacuous (Mason) or plausibly critical (SMT; unproved, §6). Theorem D is
non-critical **only because it is a family statement**: it allows the dominating 4-point family, which exists
(record, $k=4$), and forbids dominating rank $\ge4$. That is the record's covering-family dichotomy, not a proof of the
false 4-point analogue. **F2 is passed, trivially.** Theorem D cannot give a pointwise bound, and the record already
lists covering-family exclusions under item 9 as insufficient ("what remains are dense supplies with finitely many
points per curve, curves of unbounded degree, and sporadic points"). The cube family is exactly a non-covering
rational family of the kind that item excludes from coverage.

## 8. Arithmetic transfer, literature, novelty

**Transfer: the negative verdict is confirmed.**
* $U_1/U_k$ for $\mathbb Z_p$ is cyclic, so a $p$-adic arc gives one congruence. The function-field $\asymp D/2$ independent
  equations come from the infinite residue field. This is correct, and it is the real reason the model is easier.
* Over $\mathbb F_q(t)$, Lemma 1.1 fails once $L\ge p$ (Newton). ff notes the counting heuristic there; Theorem C
  survives over $\mathbb F_q$ only by composition.
* The Buium point is partly literature-checked. WebSearch confirms that the $p$-derivation is the Fermat quotient
  $\delta x=(\phi(x)-x^p)/p$ (arXiv 2212.02137). I did not find a source for the explicit $\delta$-character of $\mathbb G_m$, so that
  stays a memory-quote. The conclusion, a residue mod $p^2$ and hence point-level and faked, is consistent with
  item 1 under the record's scope limits (moduli $\le M^A$, primes above the modulus bound).

**Literature.**
* WebSearch found no function-field version of the Cilleruelo–Granville short-arc problem. Results returned:
  Cilleruelo–Granville, *Close lattice points on circles*, CJM 2009; Sedunova's function-field Bombieri–Pila (Acta
  Arith. 2017); Cilleruelo–Shparlinski short polynomial intervals over $\mathbb F_q[T]$.
* Odlyzko's 0/1 lemma: confirmed (Costello–Tao–Vu, math/0505156, Lemma 3.1).
* Plotkin, Sperner and Kollár II.3 were not re-checked. They are standard memory-quotes.

**Novelty, item by item.**
| item | status |
|---|---|
| Lemma 1.1 (linearisation) | **in record** (`families.md` Lemma P; `toric.md` §7 coupled PTE) |
| Lemma 1.2, Thm A | elementary; the function-field Cilleruelo–Córdoba; new as written, routine |
| Thm B (fakes) | new in this form; **proof gap**, repaired for $L\ge k+1$ |
| Thm C (8-point cube), C′ | **genuinely new** to the record; exact; refutes the "polynomial templates ≤6" impression |
| Thm D | **record's argument** (`gauss.md` Thm 4.3(ii), `toric.md` Prop 4.3) transplanted to $T_\Lambda$ |
| rank reduction / iff, (R$_4$) | repackaging of record item 8 ("bounded rank ⇔ uniform"); (R$_4$) is a clean new open question |
| small-case table | partly new; two entries improved here; $M(7,0)\ge4$ was already exact in the record |
| arithmetic transfer | negative, consistent with items 1, 6, 8 |

Overall the lane is a **repackaging of known barriers** (exceptional set / rank / covering families) in the
function-field model. It adds one genuinely new exact object, the cube family.

## 9. Corrected statement

Over $\mathbb C(t)$ (twisted by complex conjugation), the following hold.
1. **(proved; in record as Lemma P)** The arc condition is $L=D-1-\lfloor D/2+c\rfloor$ real-linear equations
   $\sum_k(a_k-b_k)\operatorname{Im}\rho_k^m=0$. Distinct cluster points have $\ell^1$ distance $\ge L+1$ in every
   characteristic.
2. **(proved)** $M(D,c)\le2^{2c+2}(D+2)$. For $D=2^k$ and $0\le c\le D/2-k-2$ there are fakes (arbitrary columns) with
   $\ge D+2c+2$ points. The case $c=D/2-k-1$, including $D=8$, $c=0$, is not covered: ff's proof yields repeated
   columns there.
3. **(proved, exact)** The cube family gives 8-point clusters with $c=0$ for all $D\equiv0\ (8)$, over $\mathbb Q(i)$ and
   $\mathbb F_{q^2}$. Hence $\sup_DM(D,c)\ge8$ for every $c\ge0$. It also gives an 8-point $\mathbb Z[i][u]$ template with
   $C\approx3.18\cdot10^9$.
4. **(proved; the record's free-curve argument)** A dominating family of rank-$r$ clusters satisfies $r(L+1)+1\le2D$.
   The proof needs complexification of the real parameter variety, and the case $\operatorname{rk}G=r$ needs the corrected
   argument of §5.
5. **(elementary)** $\sup_DM(D,c)<\infty$ iff the code rank is bounded; for non-squarefree $N$ use volume packing, not
   Odlyzko. Rank $\le4$ for large $D$ would give $M\le8$; this is (R$_4$), **open**. "$\limsup=8$" is a heuristic
   prediction.
6. **(exact)** Squarefree table: $M=4,3,6,3,5$ for $c=0$, $D=2..6$. $M(7,0)\ge4$ is exact (record). $M(7,1)\ge8$ and
   $M(8,1)\ge14$ are exact, from vertical-line configurations.
7. **No arithmetic consequence.** The growth rate $(\frac23+\varepsilon)\log R/\log\log R$ is unchanged.

## 10. Files (`scratchpad/round29/kill_ff_checks/`)
| file | content |
|---|---|
| `thmC.py` | Thm C identity and exact 8-point cluster over $\mathbb Q(i)$ |
| `template.py` | C′ template: 8 Gaussian integers of one norm, chord$/\sqrt R\to3.1796\cdot10^9$ |
| `thmB.py`, `thmB2.py` | Thm B fibres; column coincidence after normalisation; triple repair |
| `lp_rank.py` | Prop 3.3 LP (own implementation) |
| `L2closed.py`, `L2_7.txt` | exhaustive closed-set enumeration, $L=1,2$; $D=7$, $L=2$ run (best 8 at the last level written) |
| `L2find.py`, `rank1L2.py`, `verify_rank1.py` | vertical-line clusters: exact $M(7,1)\ge8$, $M(8,1)\ge14$ |
| `compose16.py` | $Q_{\rm cube}\circ h$ at $D=16$: 8 zero-sum sets, rank 4 |
| `vert.py` | purely imaginary locus at $c=0$: only 3 points |
Notes: `scratchpad/round29/kill_ff_notes.md`.
