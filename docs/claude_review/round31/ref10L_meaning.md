# Referee block "meaning" on Part X (signed residues on one fixed modulus)

Source refereed for meaning: `astra_lemma_X/answer.md` lines 2151-2262 (Part X), `final10.txt`, `x_checks.py` (not imported).
Context read: `round31/README_IXL.md` §§0, 1.6, 2, 3, 5; the question `astra_lemma_X/reply10.txt`; the parallel circle question
`astra/free/reply10.txt`. I did not re-referee the proofs of X.1-X.4 (other blocks do); where I rely on a step I say so.
My scripts are in `round31/ref10L_meaning/` (own code, exact integers or exact rank tests; outputs `*.out`).

**Definitions used below.** $q$ squarefree, $s=\omega(q)$, $W=\log q$, all primes in $[P,2P]$ ("dyadic parent"),
$L=C\sqrt q$. A *cut pattern* $S$ is a $B\times s$ sign matrix, column $S_p$ giving the two classes at $p$; $b_p=(S_p+\mathbf 1)/2$.
A *genuine configuration* is $B$ distinct integers $y_i$ in an interval of length $L$ with $y_i$ in class $r_p$ or $r'_p$ mod $p$
according to $S_{ip}$, $r_p\ne r'_p$, $\delta_p=r'_p-r_p$. A *sign-free table* is $(d_{ij})$ with $0<|d_{ij}|\le L$,
$\gcd(d_{ij},q)=h_{ij}$, and $d_{ij}\equiv\pm\delta_p\pmod p$ at each separating $p$, the sign free per pair and per prime.
$M_S=\{z\in\mathbb Z^B:\ z \bmod p \text{ is constant on each class of } p,\ \forall p\mid q\}$ (classes may get equal residues).

---

## 0. Results (summary)

1. **Delivered: option (iv), plus a reduction (iii) ⇒ (ii). Not delivered: (i), (ii).** Precisely: no sign-free table (hence no
   genuine configuration) on the full Paley core for every modulus and prime assignment of the IX.4 type, with margins
   $L=q^{5/8-\epsilon}$ (full core) and every row subset of size $\ge(3/5+\epsilon)m$ (two primes per column) or
   $\ge(2/3+\epsilon)m$ (one prime per column), $m\ge m_0(\epsilon,C)$. And: a sign-free table with $B\ge cs$ gives a genuine
   configuration on the same $q$ and the same $C$ with $\ge B/2$ marks by the *short* half of X.4 (one prime $p>C^2$, which every
   dyadic parent with $P>C^2$ has), and with $B-O_{c,C}(1)$ marks by the long anchor reconstruction (refereed elsewhere).
   Thus (iii) would give (ii); the factor-2 version does not depend on the delicate part of X.4. (Proved, relying on X.1, X.2.)
2. **For genuine configurations the Paley news is only the row-subset statement.** The full Paley core was already
   unrealisable by genuine configurations up to length $q^{1-O(1/m)}$ (IX.4, referee's strengthening, README_IXL §1.6). What
   is new for *genuine* configurations is the exclusion on every $\ge(3/5+\epsilon)m$ rows; what is new for the *relaxation* is
   everything in X.3.
3. **Map of the twin (fair, with a qualification).** The sign-free relaxation is the same as asking for a short square root of a
   fixed residue for each pair: $|d_{ij}|\le L$, $d_{ij}^2\equiv\Delta_{ij}\pmod q$, with $\Delta_{ij}$ fixed by $S$ and $\delta$
   (Proved, §2). Part X shows this relaxation equals the genuine problem up to a factor 2 in $B$ (or up to $O_C(1+s/B)$ marks).
   So among the relaxations in the record, nothing sits strictly between the valuation data of IX.4 and the genuine
   problem. Relaxations that keep residues only on a proper sub-product $q_E$ of $q$ have not been studied. X.1 forces
   degeneracy only when $q_E\,g>9C^2q$, so for them the question is open (§2).
4. **Column orthogonality cannot give the rate for arbitrary cuts (Proved, §6).** No argument that uses only the linear
   consequences of class-constancy mod $p$ and the interval bound can show that $y$ is constant for every cut pattern with
   $B\ge cs$. For every dyadic $q$ there are cut patterns with $2^{(1/2-o(1))s}$ distinct rows whose lattice $M_S$ contains a
   nonconstant two-valued vector of spread $\le L$ (a *shadow*), and every such linear certificate also annihilates the shadow.
   The Siegel count asked for in the block needs $B>\tfrac12(1+\tfrac{\log2}{\log P})s+O(\cdot)$. At the smooth scale
   $P=s^{1+o(1)}$, where the valuation barrier lives, that regime is already empty by the larger sieve. **No gain beyond
   the constant $\tfrac12$ comes from this route.**
5. **A structural fact (Proved, §6, probably folklore).** For $s\ge s_0(C)$, all vectors of $M_S$ with spread $\le L$, for every
   residue choice, lie in $\mathbb Z\mathbf 1+V_S$ for one rational subspace $V_S$ of dimension $\le2$ that depends only on
   $(q,S,L)$. So linear relations can never pin down more than $B-3$ independent combinations of a genuine configuration.
6. **Half-row Paley (the circle fakes' row class).** Open asymptotically. Computed: for the class
   $\{i:S_{i0}=+1\}$ at $m=19$ and $m=23$, kernel vectors of support $8<9.5$ and $11<11.5$ span the whole kernel. Her
   deletion argument therefore excludes these half-row patterns (9 and 11 rows) once $P\ge P_0(m,C)$. At $m=11$ it does
   not apply (support 6 > 5.5). Excluding the pattern for all large $m$ needs spanning kernel vectors of support
   $\le m/2-cm/\log P$. The data (spark/$m$ = 0.57, 0.55, 0.42, 0.48 at $m$ = 7, 11, 19, 23) do not decide this.
7. **Circles.** The circle analogue of the sign-free table is formulated in §4. On circles the triangle divisibility is only
   $g^2q$ against a size of order $C^4q^2$, so the twin mechanism has no analogue there (Proved for this argument;
   the conclusion that pair-level data do not force genuineness on circles is Heuristic). The explanation in item (5)
   of the brief is **correct with one correction**: the obstructing character is $\alpha/\bar\alpha$ with
   $\log|\alpha|=\tfrac{m+1}{4}\log(Q_vQ_w)\approx W/2$, and the threshold is $W/4$ (not $W/2$) (§5).
8. **Recommendation: stop the two-class thread.** The answer is a restricted-class exclusion plus a collapse of the last
   relaxation, which is the judge's stop condition in substance. If one more exchange is wanted for completeness, ask only the
   half-row Paley question (P(rate) ≈ 0, P(clear answer) ≈ 0.35) (§7).

---

## 1. What she delivered (option by option)

| Option | Delivered? | Precise form | Label |
|---|---|---|---|
| (i) $B=o(s)$ for every dyadic $q$ or one explicit dense family | **No** | stated by her (X.5 "OPEN") | — |
| (ii) $B\ge cs$ genuine on one dense $q$ | **No** | — | — |
| (iii) sign-free table with $B\gg s$ | **No**; instead a reduction: (iii) ⇒ (ii) | (a) $B_{\rm int}\le B_{\rm free}\le 2\max(4,B_{\rm int})$ whenever $q$ has a prime $p>C^2$; (b) on dyadic squarefree $q$, delete $O_C(1+s/B)$ rows to get a simultaneous signing (X.4, long proof) | (a) Proved (short; I checked X.1's four-factor divisibility and the unit-constant step); (b) Proved per her, refereed elsewhere |
| (iv) no sign-free table for Paley cuts | **Yes** | For the IX.4 Paley core ($m$ rows, $m$ columns, a bounded positive number of primes per column, primes in $[P,2P]$): none at $L=C\sqrt q$ for $m\ge m_0(C)$, any prime assignment; none at $L=q^{5/8-\epsilon}$ with equal numbers of primes per column, $m\ge m_0(\epsilon)$; none on any $\ge(3/5+\epsilon)m$ rows (two primes per column) or $\ge(2/3+\epsilon)m$ rows (one prime per column), $m\ge m_0(\epsilon,C)$ | Proved (her X.3; not re-refereed here) |

Quantifiers: X.1 holds for every squarefree $q$ and every $C$; the forcing "$H=0$" needs $g>C^2$ (unit constant). X.2 needs $b\ge5$.
X.3 depends on $C$ and $\epsilon$; $A$ (in $P\le s^A$) enters only through $\log m\le\log P+O(1)$, which holds for any $2m$ distinct
primes in $[P,2P]$, so the conclusions are uniform in $A$. X.4(b): constants depend on $C$ only; $s\ge s_0(C)$.

**Reading of "(iii) would give (ii)".** The factor-two reduction suffices for this sentence. Partition a sign-free table by its
bit at a prime $p>C^2$. Every triple inside one side has $p\mid g$, so by X.1 it degenerates. By X.2 a side with $\ge5$ marks
embeds isometrically in $\mathbb Z$. The embedded integers keep the gcds, and the class structure follows up to the unordered
pair of residues; her CRT translation step handles that. So $B_{\rm free}\ge cs$ gives $B_{\rm int}\ge cs/2$ on the same $q$ and $C$.
The anchor reconstruction of X.4(b) is needed only for the sharper "up to $O(1)$ marks".

## 2. What Part X means for the map of the twin

**Reformulation (Proved, elementary).** On a squarefree $q$, $p\mid d_{ij}$ iff $p\mid d_{ij}^2$. The sign-free conditions
(gcd $=h_{ij}$; $d_{ij}\equiv\pm\delta_p$ at separating $p$) are therefore exactly
$$ d_{ij}^2\equiv\Delta_{ij}\pmod q,\qquad \Delta_{ij}:=\mathrm{CRT}\big(\delta_p^2\,[p \text{ separates } i,j];\ 0\ \text{otherwise}\big). $$
So a sign-free table is a choice, for each pair separately, of a square root of a fixed residue $\Delta_{ij}$ in $[-L,L]$.
It is the purely pair-level datum "$d_{ij}^2\bmod q$ plus magnitudes".

**Is "no known relaxation strictly between valuation data and the genuine problem" fair?** Yes, read as "none among the
relaxations studied in the record". The relaxations are IX.4's table (magnitudes, gcds, valuation patterns, Sidon
signatures, integral cocycle, approximate real positions; it reaches the larger-sieve constant $\tfrac12$ and fails
$\Delta$), then the sign-free table, then the genuine problem. Part X collapses the second onto the third (factor 2 in $B$,
or $O_C(1+s/B)$ marks). It is **not** fair as "there is no intermediate relaxation". An obvious unstudied family keeps the
square congruences only at a set $E$ of primes. X.1's argument then gives $q_E^2g^2\mid H$ against $|H|\le81C^4q^2$, which
forces degeneracy only when $q_E\,g>9C^2q$, essentially when $E$ is everything. So "valuation data + pair squares on a
sub-product $q_E\le q^{1-\delta}$" is a genuine intermediate whose status is unknown (Heuristic: it is probably as weak as IX.4).

**Consequence for proof strategies (Proved as a logical statement).** Any proof of a bound $B\le f(s)$ for genuine twin configurations on
dyadic $q$ also bounds sign-free tables by $2\max(4,f(s))$. A proof may therefore use only pair-level data, the
$\Delta_{ij}$ and the magnitudes. This frees nothing, because those data already determine the configuration. In particular,
no proof can gain from "sign coherence" across pairs: the coherence is a consequence. The bottleneck of the rate problem
is $\binom B2$ simultaneous short-square-root conditions with one common $\delta$, which is Conjecture 9 restated pairwise.
Heuristic size count: per pair the expected number of roots is $\approx2^{|D_{ij}|}\cdot2L/(h_{ij}\cdot q/h_{ij})$, which is
tiny for balanced pairs. The constraints are individually rare, and Part X shows they interlock completely.

## 3. The Paley exclusion and the half-row pattern

* The extremal fake profile of the circle problem (Legendre/Paley) is impossible in the twin on every modulus of IX.4 type,
  for the full core and for every subset of $\ge(3/5+\epsilon)m$ rows. (Proved, hers.)
* The circle fakes use the row class $\{i:K_{i0}=+1\}$ (density $\tfrac12$), or a density-$\rho<1$ part of it. Her deletion bound
  needs $k\le(2/5-\epsilon)m$ deleted rows. The half-row class deletes $(m+1)/2$, so X.3 does not cover it; she says so.
* **What exclusion would take (Proved reduction plus computation).** Let $T$ be the ones row together with the missing rows of
  $S$. Her argument goes through, for every residue choice, as soon as $\ker T$ is spanned by integer vectors $x$ with
  $\sum_v\log Q_v$ over $\operatorname{supp}x$ below $\tfrac12W-\log(C\,\|Sx\|_1)$. For two primes per column this means
  support $<\tfrac m2\big(1-\tfrac{\log2}{\log P}\big)-\tfrac{\log(C\|Sx\|_1)}{2\log P}$. By X.1/X.2 any table on $\ge5$ such rows is
  genuine, so the line lemma step is free. The only issue is sparse spanning circuits of $T$.
  **Computed** (`paley_half.out`, `paley_half2.out`; $S=Q-I$, $Q_{ij}=\chi(j-i)$, $m\equiv3$ mod 4; exact Gram identity checked):

  | $m$ | retained rows | $T$ | spark$(T)$ | generic | span of minimum circuits | max coefficient |
  |---|---|---|---|---|---|---|
  | 7 | 3 | 5×7 | 4 | 6 | — | — |
  | 11 | 5 | 7×11 | 6 | 8 | 4 = dim ker | 2 |
  | 19 | 9 | 11×19 | 8 | 12 | 8 = dim ker | 3 (7 at support 9) |
  | 23 | 11 | 13×23 | 11 | 14 | 10 = dim ker | 13 |

  The spark is well below the generic value, which reflects the quadratic-residue symmetry. At $m=19$ (support 8) and $m=23$
  (support 11) spanning circuits lie below $m/2$. Her inequality then holds for $P\ge P_0(m,C)$ (for example $\log P>24+\log C$
  at $m=23$ with two primes per column). So the half-row Paley pattern at $m=19,23$ admits no sign-free table and no genuine
  configuration once $P$ is large. **Computed + Proved (her argument)**, but only at fixed $m$ and large $P$: 9 or 11 marks, not
  an asymptotic statement. At $m=11$ the minimum support 6 exceeds 5.5 and the argument fails. Asymptotically one needs spanning
  circuits of support $\le m/2-cm/\log P$. **Gap/open**: the data do not show the trend.
* **Realising it.** A sign-free (equivalently genuine) half-row table would be a genuine configuration with $\approx m/2=s/4$
  marks (two primes per column). That contradicts the rate form of Conjecture 9 on that $q$. There is no reason to expect one
  (Heuristic). The random-parent theorem gives $O(\log s)$ on typical parents. The expected number of $B$-root clusters for fixed
  $q$, summed over $a$, is about $q\,(2^sL/q)^{B-1}=q\,2^{s(B-1)}(C/\sqrt q)^{B-1}$. Its logarithm is
$s\big[(B-1)\log2-\tfrac{B-3}2\log P\big]+O_C(B)$, which is negative for $B\ge4$ once $\log P>2(B-1)\log2/(B-3)$ ($P>64$ at $B=4$). And §6 shows that all short
  vectors of $M_S$ lie in a 3-dimensional space. Worst-case arithmetic structure could beat these counts, but no mechanism for
  $B\gg s$ is visible.

## 4. Circles: the analogue of the question just answered

**Pair-level data (Proved identities; `circle_checks.out`, 900 pairs, 0 failures).** Take a unit class, a squarefree split norm $q$,
and $z_i=\eta\prod_p\pi_p^{b_{ip}}\bar\pi_p^{1-b_{ip}}$. For a pair, $G_{ij}:=z_i\bar z_j=u\,h_{ij}X_{ij}^2$ with $u$ a unit and $X_{ij}=x_{ij}+it_{ij}$ of
norm $Q_{ij}=q/h_{ij}$. Hence $\mathrm{Re}\,G_{ij}=h(x^2-t^2)$, $|\mathrm{Im}\,G_{ij}|=2h|xt|$, and squared chord $=4h_{ij}t_{ij}^2$ (up to the
unit's swap of $x,t$).

**The circle "sign-free table".** For each pair, the difference set $D_{ij}$ and a representation $Q_{ij}=x_{ij}^2+t_{ij}^2$ with
$|t_{ij}|\le\tfrac C2\sqrt{Q_{ij}}\,q^{-1/4}$. Equivalently, an orientation at each $p\in D_{ij}$ chosen independently per pair and per prime
(the $2^{|D_{ij}|-1}$ representations modulo conjugation). This is the circle counterpart of the twin's independent signs.
Genuine means that these orientations come from one bit vector per vertex, i.e. $G$ has rank one: $G_{ij}G_{jk}=qG_{ik}$.
On a triangle this reads $X_{ij}X_{jk}=\pm A_jX_{ik}$, $A_j=\prod_{p\text{ separating }j}p$, which forces conjugate primes at $p\mid A_j$.
So for triangles, rank one over $\mathbb Z[i]$ already gives orientation consistency (Proved).

**Does pair-level data force genuineness? The twin mechanism does not transfer (Proved for this mechanism).** The real triangle
condition is $\det$ of the $3\times3$ cosine Gram $=0$. One has $q^3\det=q^3+2R_{12}R_{23}R_{13}-q\sum R^2$ with $R=\mathrm{Re}\,G$. At a constant
prime $p^3$ divides it; at a separating prime only $p^1$ is forced. The divisor is therefore $g^2q$. The size is
$|q^3\det|\le\tfrac34C^4q^2(1+o(1))$, from $\det=4\sin\sigma\prod\sin(\sigma-\varphi)$ and $\varphi\le Cq^{-1/4}(1+o(1))$. Since
$g\le h_{ij}\le C^2\sqrt q/4$, we get $g^2q\le C^4q^2/16$, so this comparison never forces $\det=0$, for any $C$. In the twin the
divisor is $q^2g^2$ against $81C^4q^2$. The circle loses a full factor $q$ because separated squared chords are not constant
mod $p$ (her norm-65 example). **Heuristic:** pair-level circle data probably do not force genuineness, and the twin result says
nothing about the parallel session's option (A2) ("integral-cofactor fake"). The twin result neither predicts nor excludes
such a fake. Per pair, a short-$t$ representation exists with "probability" about $2^{|D_{ij}|}Cq^{-1/4}$, so a fake must
satisfy $\binom M2$ rare near-square conditions. The difficulty is in existence, not in an obstruction. Do not cite Part X as
evidence either way for circles.

## 5. Why the twin Paley argument does not port to circles (item (5) of the brief): confirmed with a correction

For $a=S_v-S_w$: $\sum_ia_i=0$, $a\cdot b_h=0$ for $h\ne v,w$, $a\cdot b_v=(m+1)/2$, $a\cdot b_w=-(m+1)/2$, and $\|a\|_1=m+1$.
(Proved from $S^tS=(m+1)I-J$; Computed at $m=11,19,43$, `circle_checks.out`.) Hence
$$\prod_i z_i^{a_i}=\Big(\tfrac{\Pi_v}{\bar\Pi_v}\Big)^{\frac{m+1}2}\Big(\tfrac{\bar\Pi_w}{\Pi_w}\Big)^{\frac{m+1}2}=\frac{\alpha}{\bar\alpha},\qquad
\alpha=\Pi_v^{\frac{m+1}2}\bar\Pi_w^{\frac{m+1}2},\quad \Pi_v=\prod_{p\in Q_v}\pi_p .$$
The exponent vector is supported on the primes of columns $v,w$, with exponent $\pm(m+1)/2$ on each $\pi_p$ (and the opposite on
$\bar\pi_p$). The arc gives $|\arg(\alpha/\bar\alpha)|\le\tfrac{m+1}2Cq^{-1/4}$. The monomial gap gives $|\alpha/\bar\alpha-1|\ge2/|\alpha|$
unless $\alpha/\bar\alpha=1$. Triviality is therefore forced iff $\log|\alpha|<W/4-\log((m+1)C/4)$. But
$$\log|\alpha|=\tfrac{m+1}4\log(Q_vQ_w)=\tfrac W2\big(1+\tfrac1m\big)\big(1+O(1/\log P)\big)$$
for balanced columns, with any fixed number of primes per column. So $|\alpha|\approx q^{1/2}$ against the needed $q^{1/4}$.
**Correction to the brief's wording.** "Weight about $(m+1)\log P>W/2$" is right only with $\log|\alpha|$ as the weight and two
primes per column. The correct comparison is $\log|\alpha|\approx W/2$ against the threshold $W/4$, or equivalently
$\log N(\alpha)\approx W$ against $W/2$. It fails by a factor 2 on the log scale, not marginally.
**Mechanism (confirmed).** The twin uses the additive form $a^ty$, whose divisibility comes from the *other* columns
($q/(Q_vQ_w)$). Its cost is linear in $\|a\|_1$. On circles the other columns give exponent 0 (no information), and the two
columns $v,w$ carry exponent $(m+1)/2$, whose cost is exponential. No additive substitute exists on circles. Mod $\pi_p$ the
class not containing $\pi_p$ has nonconstant residues, so $\sum a_iz_i$ gets no divisibility from $a\cdot b_p=0$.

## 6. Column orthogonality and the Siegel count (the main new analysis of this block)

**Basic lemma (Proved).** If $a\in\mathbb Z^B$, $\sum a_i=0$ and $\sum_{i\in\text{class 2 at }p}a_i\equiv0\pmod p$ for $p\in G$, then for
every $y\in M_S$, $\prod_{p\in G}p\mid a^ty$. If moreover $\|a\|_1L<2\prod_{p\in G}p$ and the spread of $y$ is $\le L$, then $a^ty=0$.
The exact version $a^tS_p=0$ is the special case.

**Siegel count (Proved; this is the computation the brief asks for).** Exact orthogonality on a set $G$ of $g$ columns is $g+1$
equations with $\pm1$ coefficients. A nonzero solution exists when $B\ge g+2$, with
$\|a\|_\infty\le B^{(g+1)/(B-g-1)}$ (Siegel / Bombieri–Vaaler). The divisibility beats the size iff
$$ g\log P>\tfrac s2\log(2P)+\log\big(CB\|a\|_\infty\big),\quad\text{i.e.}\quad g>\tfrac s2\Big(1+\tfrac{\log2}{\log P}\Big)+O\Big(\tfrac{\log B}{\epsilon\log P}\Big)$$
when $B-g\ge\epsilon s$. So column orthogonality yields any relation only when $B>\tfrac12\big(1+\tfrac{\log 2}{\log P}\big)s+O_\epsilon(1)$.
For $B<s/2$ and $S$ in general position there is no kernel at all.

**Comparison with the larger sieve (Proved).** The Vandermonde form gives $B\le(\tfrac12+o(1))W/\log B$. At $P=s^{1+o(1)}$, which
includes the IX.4 scale $P\asymp s(\log s)^2$ where the valuation barrier is attained, this is $B\le(\tfrac12+o(1))s$. The regime
where column orthogonality yields any relation is therefore **already empty** there. At $P=s^A$ with $A>1$ the larger sieve allows
up to $(\tfrac A2+o(1))s$, and column orthogonality produces relations for $B>(\tfrac12+\epsilon)s$. But (next two results) these
relations never force $y$ constant for general $S$.

**Proposition M (Proved; probably folklore).** Let $k$ be the dimension of the span of the projections to $\mathbf1^\perp$ of all
$z\in M_S$ with spread $\le L$. If both classes are nonempty at every $p$, then with $n=B-1$ and $V_n$ the volume of the unit ball,
$$ \big(\tfrac k2-1\big)\log q\ \le\ k\log\!\big(\tfrac{C\sqrt B}2\big)+\log\!\big(\tfrac{n!V_n}{2^n}\big)+\tfrac12\log B .$$
*Proof.* $[\mathbb Z^B:M_S]=q^{B-2}$. Since $\mathbf1\in M_S$ is primitive, the projection $N$ has covolume $q^{B-2}/\sqrt B$.
Next, $\lambda_i(N)\le q$ (from $qe_i$) and $\lambda_1,\dots,\lambda_k\le\sqrt BL/2$. Apply the lower bound in Minkowski's second theorem.
*Consequence.* Since $\log(n!V_n/2^n)\le\tfrac n2\log n$, we get $k\le2+(B\log B+O(B))/W$. For a genuine configuration the larger
sieve gives $B\log B\le(\tfrac12+o(1))W$ (primes $>B$), hence $k\le2$ for $s\ge s_0(C)$. So every genuine configuration with cut
pattern $S$, for every residue choice, lies in $\mathbb Z\mathbf1+V_S$ with $\dim V_S\le2$ determined by $(q,S,L)$.
Linear relations can therefore give at most $B-3$ independent annihilators of a genuine $y$.
*Computed sanity check:* on her five-mark example ($q=59\cdot67\cdot71\cdot83$, $C=5.91$) LLL gives Gram–Schmidt norms
$40661,129092,179067,13449391$ against $R=57084$. So $k=1$, and the short vector is $y$ itself (`minima2.out`).

**Shadow obstruction (Proved; the key negative result).** For every dyadic $q$ ($P\ge2$) let $D$ be a set of
$\lfloor\tfrac s2\cdot\tfrac{\log P}{\log 2P}\rfloor$ primes, and $d=\prod_Dp\le\sqrt q\le L/C$. Take any cut pattern whose rows carry
arbitrary bits at the primes of $D$ and, at every $p\notin D$, the bit $[z_i=d]$ for a fixed partition of the rows into $\{z=0\},\{z=d\}$.
Then $z\in M_S$ is nonconstant with spread $d\le L$ ($C\ge1$), and $S$ can have $2^{|D|+1}=2^{(1/2-o(1))s}$ distinct rows.
Every linear certificate (any $a$ as in the basic lemma, any $G$, any heights) satisfies $a^tz=0$, so the certificates never
span $\mathbf1^\perp$. Exactly: for $p\notin D$, $a\cdot S_p=0$ and $a\cdot\mathbf1=0$ give $\sum_{z_i=d}a_i=0$. Hence **no theorem of the
form "$B=o(s)$ unless $S$ has no short dependencies" can be proved by column orthogonality for arbitrary $S$**. The cut patterns
it cannot handle carry no information about $B$: they are excluded, if at all, only by nonlinear inputs (distinct residues,
the Vandermonde product). By Banaszczyk transference ($\lambda_i(N)\lambda_{n+1-i}(N^*)\le n$, $\Lambda^{\rm mod}=qN^*$), the linear
method succeeds on $S$ **iff**, up to a factor $B^{O(1)}$ in $L$, $M_S$ has no nonconstant short vector. That is the problem
with nondegeneracy dropped (Proved).

**Is there any gain beyond the valuation barrier $\tfrac12$? No, not from this route (Proved for the linear method; Heuristic
in general).** The first use of additive structure beyond the sieve constant would have to combine Proposition M (a 3-dimensional
family) with a nonlinear input, i.e. the distinct-residue condition or the quadratic relation $y_i^2+by_i+c\equiv0\pmod q$. I see no
such combination: in $y=c+\alpha u+\beta w$ the vectors $u,w$ have spread $\asymp L\gg P$, so their reductions mod $p$ carry no
smallness. Label: Gap.

## 7. Recommendation

**Stop the Question B two-class thread.** Reasons:
1. The judge's condition: Part X delivered (iv), an exclusion for a restricted class (Paley), plus a structural collapse. The
   rate did not move. Her X.5 says so plainly.
2. The last relaxation is gone. Sign-free residues are the genuine problem (factor 2). There is no barrier left to map
   between IX.4 and the truth, except unstudied sub-product variants that do not bear on the rate.
3. The natural successor route, column orthogonality / sparse near-kernel vectors for all cut patterns, is closed by the
   shadow obstruction (§6). It cannot give $o(s)$ for arbitrary $S$, and in the only regime where it yields relations it is
   dominated by the larger sieve at the smooth scale.

If one closing exchange is wanted anyway, ask **one** question:
> On the Paley core of IX.4, let $T$ be the matrix of the all-ones row and the rows with $S_{i0}\ne+1$. Prove or disprove that
> for all large $m\equiv3\pmod 8$ the kernel of $T$ is spanned by integer vectors of support $\le(\tfrac12-c)m$ and height
> $m^{O(1)}$. Equivalently, decide whether the half-row Paley pattern (the row class of the circle fakes) is realisable in the
> twin. Known: spark$(T)$ = 6, 8, 11 at $m$ = 11, 19, 23, and the minimum circuits span at 19 and 23.

Honest probabilities: P(rate for worst-case $q$) ≈ 0 (<0.005); P(clear answer to the Paley question) ≈ 0.35. Its only value is
to close the comparison with the circle fakes. The twin's rate problem is Cilleruelo–Granville's Conjecture 9 (rate form).
Remaining effort is better spent on the parallel circle session's option (A), where the twin mechanism does not apply (§4).

## 8. Proposed wording for the expert note (§5 item 12, after the Part IX bullet)

> * ***Part X of the Question B thread: signed residues on one fixed modulus*** (`round31/README_XL.md`; referee blocks incl.
>   `ref10L_meaning.md`). In the critical two-class problem ($q$ squarefree, interval $C\sqrt q$), a *sign-free table*
>   prescribes for each pair an integer $d_{ij}$, $0<|d_{ij}|\le C\sqrt q$, with $d_{ij}^2\equiv(y_j-y_i)^2\pmod q$ (sign free per pair
>   and per prime). **Proved:** on every triple, $q^2g^2$ divides the Heron product, $g$ the product of the primes constant on
>   the triple. A nondegenerate triple therefore has a side $\ge\sqrt{qg}$. Every triple inside one class of a prime $p>C^2$
>   degenerates, and $\ge5$ marks with degenerate triples embed in $\mathbb Z$. Hence the maximal sign-free and genuine counts
>   satisfy $B_{\rm int}\le B_{\rm free}\le2\max(4,B_{\rm int})$, and on dyadic parents (all primes in $[P,2P]$) a sign-free
>   table of $B$ rows becomes genuine after deleting $O_C(1+s/B)$ rows. No sign-free table (so no genuine configuration)
>   exists on the Paley cut pattern of IX.4 for any modulus and prime assignment, at length $q^{5/8-\epsilon}$ on the full core,
>   and at critical length on every $\ge(3/5+\epsilon)m$ of its rows. For genuine configurations only the row-subset statement is
>   new. **Open:** a worst-case $o(\omega(q))$ bound, a genuine example of order $\omega(q)$, the half-row Paley pattern (the
>   circle fakes' row class; excluded computationally at $m=19,23$ for large $P$). **Meaning:** the pair-level square
>   congruences already force genuineness, so in the twin no relaxation studied in the record lies strictly between the
>   valuation data of IX.4 (larger-sieve constant $\tfrac12$) and the genuine problem. Linear (column-orthogonality) relations
>   cannot give a rate for arbitrary cut patterns: patterns with $2^{(1/2-o(1))s}$ rows admit nonconstant short "shadow"
>   vectors. The mechanism does not transfer to circles: separated squared chords are not constant mod $p$, the circle
>   triangle divisibility is only $g^2q$, and the Paley character has $\log|\alpha|\approx W/2$ against the gap threshold $W/4$.
>   **The growth rate, the uniform bound and Questions A and B are unchanged.**

**Sentences to avoid:** "Part X proves Conjecture 9 for Paley moduli" (it excludes one cut pattern, not all residue
configurations on a modulus); "the sign-free relaxation is a new barrier" (it collapses onto the genuine problem); "Part X shows
there is no intermediate relaxation" (only none among those studied); "Paley/Legendre profiles are impossible in the twin"
(without "full core or $\ge3/5$ of the rows"); "the circle fakes are impossible in the twin" (they use half the rows: open);
"(iii) and (ii) are equivalent" (loss of a factor 2 or $O_C(1+s/B)$ marks); "column orthogonality gives a rate beyond the sieve
constant"; "the twin result implies pair-level circle data force genuineness" (or the converse); "Part X gives
$B=o(s)$ on some dense modulus".

## 9. Files

* `ref10L_meaning/minima.py`, `minima2.py`, `minima2.out`: Proposition M sanity checks (LLL; her five-mark example gives $k=1$).
* `ref10L_meaning/paley_half.py`, `paley_half.out`, `paley_half2.py`, `paley_half2.out`: spark and circuit span for the half-row
  Paley pattern, $m=7,11,19,23$.
* `ref10L_meaning/circle_checks.py`, `circle_checks.out`: Paley character exponents ($m=11,19,43$); the Gram identity
  $z_i\bar z_j=u\,hX^2$ on 900 pairs; $q^3\det=0$ on 300 genuine triples.
