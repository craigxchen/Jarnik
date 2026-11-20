# Judge report, round 31: Part X of the Question B thread (Astra, "Signed residues on one fixed modulus")

Source: `astra_lemma_X/answer.md` from line 2151 (`## Part X.`) to the end (line 2262), with `final10.txt`,
`x_checks.py` and the question `reply10.txt`.
Referee reports in this directory, **all five present**: `ref10L_heron.md` (X.1, X.2), `ref10L_paley.md` (X.3),
`ref10L_recon.md` (X.4 and the two limiting examples), `ref10L_port.md` (X.5: circles, regressions, limiting examples),
`ref10L_meaning.md` (novelty, meaning, the column-orthogonality route, next steps).
Record consulted: `round31/README_IXL.md` (format; rulings on Part IX, in particular §1.6 on the IX.4 relaxation and §5 on
the question that produced Part X); `expert_note/expert_note.md` §5 item 12 (lines 929-1410).
Judge scripts (my own code; nothing imported from Astra or any referee), in `README_XL/`:
* `judge10L_heron.py` (outputs `judge10L_heron.out`, `judge10L_heron_b.out`): the divisibility $q^2g^2\mid\mathcal H$, the
  quotient modulo 4 and the unit-constant bound on every locally admissible triple, all squarefree $q\le2500$ with sides
  $\le3\sqrt q$ (10,307 triples) and all squarefree $q\le12000$ with sides $\le1.3\sqrt q$ (4,521 triples).
* `judge10L_eisen.py` (output `judge10L_eisen.out`): a family showing the unit constant of X.1 is asymptotically sharp
  already at $g=1$.
* `judge10L_recon.py` (output `judge10L_recon.out`): the five-mark example; $B_{\rm free}$, $B_{\rm int}$ and the 8-clique at
  $q=59\cdot67\cdot71\cdot83$; statement 1 of X.4, uniqueness of magnitudes and the complementary-pair count on 180 random
  moduli; an attempt to break the final step of X.4.
* `judge10L_lattice.py` (output `judge10L_lattice.out`): LLL on the certificate lattices of the meaning referee's item (6),
  for random cut patterns, shadow patterns and the five-mark example.
* `judge10L_halfrow.py` (outputs `judge10L_halfrow.out`, `judge10L_halfrow23.out`): spanning circuits of the half-row Paley
  matrix at $m=7,11,19,23$.

Claims of my own are labelled **[J]** and are unrefereed.

**Notation and definitions used throughout.**
* **Twin (critical two-class problem), squarefree case.** $q$ squarefree, $s=\omega(q)$, $W=\log q$, $L=C\sqrt q$. At each
  $p\mid q$ two residues $r_p\ne r'_p$ are prescribed, $\delta_p=r'_p-r_p\not\equiv0$. A *genuine configuration* is a set
  of $B$ distinct integers in an interval of length $L$, each $\equiv r_p$ or $r'_p\pmod p$ for every $p$. Its *cut
  pattern* is the $B\times s$ bit matrix (bit = class at $p$); its rows are the *signatures*. *Dyadic parent*: all primes
  of $q$ in $[P,2P]$.
* **Agreement product** $h_{ij}$: product of the primes at which rows $i,j$ have the same bit; $D_{ij}=\log(q/h_{ij})$. A pair
  is *complementary* if $h_{ij}=1$. For a triple, $g$ is the product of the primes at which all three bits agree.
* **Sign-free table** on a cut pattern (with $\delta$ fixed): integers $d_{ij}=d_{ji}$ with (T0) $0<|d_{ij}|\le L$;
  (T1) $p\mid d_{ij}$ if $i,j$ have the same bit at $p$; (T2) $d_{ij}\equiv\pm\delta_p\pmod p$ if they are separated at $p$,
  the sign free for every pair and every prime. Equivalently (T1)+(T2): $\gcd(d_{ij},q)=h_{ij}$ and
  $d_{ij}^2\equiv\Delta_{ij}\pmod q$ with $\Delta_{ij}$ the CRT residue that is $\delta_p^2$ at separating $p$ and $0$
  elsewhere. Only the magnitudes $a_{ij}=|d_{ij}|$ matter.
* $B_{\rm free}(q,C,\delta)$, $B_{\rm int}(q,C,\delta)$: the maximal size of a sign-free table, resp. of a genuine
  configuration (optimising over the choice of the two classes with difference $\pm\delta_p$ and over the interval).
* **Heron product** $\mathcal H(a,b,c)=(a+b+c)(-a+b+c)(a-b+c)(a+b-c)$; a triple is *degenerate* if $\mathcal H=0$ (the
  largest magnitude is the sum of the other two).
* **Paley core** (IX.4): $m$ prime, $m\equiv3\pmod 4$, $K_{iv}=\chi(v-i)$ (Legendre symbol mod $m$), $S=K-I$; column $v$
  carries $Q_v$, a product of $t_v\in[1,t]$ distinct primes, all primes distinct and in $[P,2P]$, $q=\prod_vQ_v$; the
  cut at $p\mid Q_v$ is the sign of $S_{iv}$. $S^tS=SS^t=(m+1)I-J$.
* $M_S=\{z\in\mathbb Z^B:\ z\bmod p$ constant on each class of $p$, for all $p\mid q\}$ (the two classes may get equal
  residues). *Spread* of $z$: $\max z-\min z$.
* **Conjecture 9** (Cilleruelo–Granville): uniformly in $a$, intervals of length $\sqrt b$ contain $O(1)$ roots of
  $x^2\equiv a\pmod b$; its *rate form* asks for $o(\omega(b))$. The critical twin is equivalent to it (README_IXL §1.5).

## 0. Verdict (five lines)

1. Nothing in Part X is false and nothing of Astra's mathematics is withdrawn. Five referees returned no FALSE verdict on a claim of hers (the one FALSE in the paley block concerns a coverage claim she did not make). I checked myself [J] the divisibility $q^2g^2\mid\mathcal H$ and the unit-constant corollary of X.1, the X.4 reconstruction line by line (with a computational attempt to break it), and the meaning referee's item (6); the last I partly overturn.
2. **Option (iv) was delivered, and a reduction (iii) $\Rightarrow$ (ii); not (i), not (ii).** No sign-free table exists on the full Paley core of IX.4 at length $C\sqrt q$ for $m\ge m_0(C,t)$ (any modulus, any prime assignment), nor at $q^{5/8-\epsilon}$ with equal numbers of primes per column, nor at critical length on any $\ge(3/5+\epsilon)m$ rows (two primes per column) or $\ge(2/3+\epsilon)m$ rows (one prime). For *genuine* configurations only the row-subset statement is new.
3. **The sign-free relaxation is the genuine problem**: $B_{\rm int}\le B_{\rm free}\le2\max(4,B_{\rm int})$ whenever $q$ has a prime $p>C^2$; on dyadic parents any sign-free table becomes genuine after deleting $O_C(1+s/B)$ marks (at most $1$ of them for complementary pairs if the signatures are Sidon [J]). Pair-level square congruences force sign coherence; nothing lies strictly between the IX.4 valuation data and the twin among the relaxations studied.
4. Corrections: the line lemma X.2 is Menger's theorem (1928); X.1's unit constant is sharp even at $g=1$, so nondegenerate $g=1$ triangles occur below $C\sqrt q$ for **every** $C>1$ [J] (not only $C\ge\sqrt{15/7}$); the meaning referee's "column orthogonality is closed as a route" is **not** established: his shadow patterns carry at most $2\lfloor C\rfloor+2$ genuine marks, and congruence certificates exist for every pattern with $B\ge4$ [J]. What survives is Proposition M and the reformulation "the linear method works on $S$ iff $M_S$ has no nonconstant short vector".
5. **For the integer problem nothing moves**: the worst-case growth rate in the twin and on circles, the uniform bound and Questions A and B are unchanged; the mechanism does not transfer to circles. The answer is a restricted-class exclusion plus a collapse, which meets the previous judge's stop rule: **stop the two-class thread**. Score on the goal: 2/10.

## 1. Per block: referee-corrected statements

### 1.1 X.1, the triangle obstruction (referee: heron; rechecked by the judge)

**Proved (statement as corrected).** Let $q>1$ be squarefree and let $a,b,c$ be the magnitudes of a sign-free table on
three distinct marks, $g$ the product of the primes at which the three bits agree. Then
$$q^2g^2\mid\mathcal H(a,b,c).$$
Moreover, if $\mathcal H(a,b,c)\ne0$ then $\max(a,b,c)\ge\sqrt{qg}$, strictly unless $g=q$ (equality at $(q,q,q)$).
Consequently: (a) at $C=1$ every triple of every sign-free table degenerates; (b) for every $C$, every triple whose bits
agree at some prime $p>C^2$ degenerates. Only (T1) and the "$\pm\delta_p$" part of (T2) are used.

*Judge's line-by-line check [J].*
* At a prime of $g$ all three magnitudes are $0\bmod p$, so $p^4\mid\mathcal H$. At a split prime, say $p\mid a$, the other
  two satisfy $b\equiv\pm\delta_p$, $c\equiv\pm\delta_p$, hence $b\equiv c$ or $b\equiv-c$; in the first case
  $a-b+c\equiv a+b-c\equiv0$, in the second $a+b+c\equiv-a+b+c\equiv0$; at $p=2$ all four factors are even. So $p^2\mid\mathcal H$.
  Squarefreeness gives $\prod_{p\nmid g}p^2\prod_{p\mid g}p^4=q^2g^2$. Correct.
* $g\mid a,b,c$ (each prime of $g$ divides all three, $g$ squarefree). With $u,v,w=a/g,b/g,c/g$ and $Q=q/g$:
  $\mathcal H(u,v,w)=\mathcal H(a,b,c)/g^4$ and $q^2g^2=Q^2g^4$, so $T=\mathcal H(u,v,w)/Q^2\in\mathbb Z$.
  $\mathcal H(u,v,w)=4u^2v^2-(u^2+v^2-w^2)^2\equiv-X^2\equiv0,3\pmod 4$; for odd $Q$, $Q^2\equiv1\pmod8$ gives
  $T\equiv\mathcal H\equiv0,3$. For even $Q$ the prime 2 is split on the triple, all four factors are even, $16\mid\mathcal H(u,v,w)$,
  $Q=2Q'$, $Q'$ odd, so $4\mid T$. Correct.
* With $l=w=\max$, $x=u^2,y=v^2\in(0,l^2]$: $\mathcal H=2l^2(x+y)+2xy-x^2-y^2-l^4$, $\partial_x\mathcal H=2(l^2+y-x)>0$; so
  $-l^4<\mathcal H\le3l^4$. If $l<\sqrt Q$ then $-1<T<3$, i.e. $T\in\{0,1,2\}$, and the congruences force $T=0$. Correct,
  including the negative (non-triangle) case.
* At $C=1$: if $g>1$, $\max\le\sqrt q<\sqrt{qg}$; if $g=1$, an integer $\le\sqrt q$ is $<\sqrt q$ since $q$ is not a square.
  Correct.

**Computed [J]** (`judge10L_heron.out`, `judge10L_heron_b.out`). All 1,522 squarefree $q\le2500$, all triples of positive
lengths $\le3\sqrt q$ admissible at every prime: 10,307 triples, 374 nondegenerate; 0 failures of $q^2g^2\mid\mathcal H$, 0 of the
mod-4 statement, 0 nondegenerate triples with $\max^2<qg$. All 7,292 squarefree $q\le12000$ at lengths $\le1.3\sqrt q$: 4,521
triples, 0 failures. Astra's census (182 moduli, 252 triples, 13 nondegenerate, $q\le300$, lengths $\le2\sqrt q$) is reproduced
exactly. (The referee's 914 moduli / 5,914 triples at $q\le1500$, cap $3\sqrt q$, is consistent.)

**Correction [J] to the referee's precision "$g=1$ nondegenerate triples need $C\ge\sqrt{15/7}$".** That ratio is an artefact of
the range $q\le1500$. My census finds $\max^2/(qg)=55/39$ at $q=2145$, sides $(39,49,55)$, and $119/95$ at $q=11305$, sides
$(95,109,119)$. These are Eisenstein triples: if $a^2-ac+c^2=b^2$ and $q=ac$ is squarefree, then $(a,b,c)$ is locally
admissible (at $p\mid a$: $b^2\equiv c^2$; at $p\mid c$: $b^2\equiv a^2$), $g=1$, and $\mathcal H=3a^2c^2=3q^2$ (Proved [J], one line).
With $a=m^2-n^2$, $c=2mn-n^2$ and $m/n\to2$ one gets $\max^2/q\to1$; computed (`judge10L_eisen.out`): 9,024 squarefree cases with
$m<400$, the best $\max^2/q=1.00676$ at $q=13{,}832{,}190{,}505$, sides $(117215,117613,118007)$. **So the unit constant of X.1 is
asymptotically sharp even at $g=1$, and for every fixed $C>1$ there are nondegenerate $g=1$ sign-free triangles of length
$\le C\sqrt q$** (for suitable $q$). This strengthens, and does not contradict, the reason why X.4 needs a prime $p>C^2$.

**Verdict.** Correct as stated. Gap: none. Withdrawn: none.

### 1.2 X.2, the line lemma (referee: heron)

**Proved; known.** Let $b\ge5$ points carry positive symmetric integer distances with every triple degenerate. Then there are
distinct integers $x_i$ with $d_{ij}=|x_i-x_j|$ (take $x_i=d_{ui}$ for one end $u$ of a diametral pair). For four points the
only other case is the *pseudo-linear quadruple*: a 4-cycle with sides $x,l-x,x,l-x$ and both diagonals $l$.
**This is a special case of Menger's theorem** (K. Menger, Math. Ann. 100 (1928); L. M. Blumenthal, *Theory and
Applications of Distance Geometry*, 1953): a semimetric space with more than four points embeds isometrically in $\mathbb R$ if
every triple does; the pseudo-linear quadruples are the four-point exceptions. Astra's proof is a correct self-contained version
and must cite Menger. [J: I did not consult the full texts; the attribution is standard and the statement matches.]

**Computed (referees).** heron: every all-degenerate table on 5 points (entries $\le20$), 6 ($\le14$), 7 ($\le6$) is a line metric;
on 4 points (entries $\le24$) the 828 non-line tables of 25,116 are exactly the pseudo-linear ones. paley: independent
enumeration, same result.

**Proved (Astra; checked by port).** If four signatures are Sidon, the four-point exception cannot occur: opposite edges are
equal, so every prime has an even number of zero edges, so no 3|1 cut, so the four signatures sum to $0$.

**Verdict.** Correct; citation missing. Withdrawn: none.

### 1.3 X.3, Paley cuts on a fixed modulus (referee: paley)

**Proved (statement as corrected; all quantifiers).** Let $m$ be prime, $m\equiv3\pmod4$ ($m\equiv3\pmod 8$ is **not** needed; it
served only the Sidon property of IX.4). Let $q=\prod_vQ_v$, $Q_v$ a product of $t_v\in[1,t]$ distinct primes, all distinct and in
$[P,2P]$, $P$ arbitrary, any assignment, any nonzero $\delta_p$.
1. *(Counting.)* Two distinct rows are separated at exactly $(m+1)/2$ columns, and every triple has exactly $(m-3)/4$ constant
   columns, so $g_{\min}\ge P^{(m-3)/4}$ on every row subset.
2. *(Reduction, any row set $R$ with $|R|\ge5$.)* If $L<\sqrt{q\,g_{\min}(R)}$ (in particular $L=C\sqrt q$ with
   $C^2<P^{(m-3)/4}$), every sign-free table on $R$ equals $|y_i-y_j|$ for a genuine two-class configuration $y$ on $R$ with
   the Paley cut, of diameter $\le L$; the two prescribed classes may be interchanged at some primes.
3. *(Full core.)* For $C\ge1$ there is $m_0(C,t)$ such that for $m\ge m_0$ no sign-free table on all $m$ rows exists with
   $L=C\sqrt q$. Sufficient: $P^{(m-3)/4}>C^2$ and $C(m+1)(2P)^{2t}<2P^{m/2}$. **No dependence on $A$** ($P\le s^A$ is never used).
4. *(Margin.)* If $t_v=t$ for all $v$: for every $\epsilon>0$, no sign-free table at $L=q^{5/8-\epsilon}$ for $m\ge m_1(\epsilon,t)$.
   If only $t_v\in[1,t]$: none at $L=q^{1/2+1/(8t)-\epsilon}$.
5. *(Missing rows.)* If $t_v=t$ for all $v$: for every $\epsilon>0$, $C\ge1$ and $m\ge m_2(\epsilon,t,C)$, no sign-free table at
   $L=C\sqrt q$ on any row set $R$ with $|R|\ge\big(\frac{t+1}{2t+1}+\epsilon\big)m$, i.e. $\frac35$ for $t=2$, $\frac23$ for $t=1$;
   if $\log m\le\theta\log P$ the threshold improves to $\big(1-\frac t{2t+\theta}+\epsilon\big)m$.

*Proof structure (checked by the referee; I agree).* Item 1 by counting separations ($2\cdot\#\text{nonconstant}=3(m+1)/2$).
Item 2 by X.1 and X.2, then: same-class pairs have $p\mid d_{ij}$, cross pairs have $p\nmid d_{ij}$ (gcd condition), so each class
has one residue and the two differ by $\pm\delta_p$. Item 3: $a=S_v-S_w$ satisfies $a^t\mathbf 1=0$ and $a^tS_h=0$ for $h\ne v,w$, so
$q/(Q_vQ_w)\mid a^ty$; $a$ has exactly $(m+1)/2$ nonzero entries, each $\pm2$, so $|a^ty|\le\frac{m+1}2L$; hence $S^ty$ is
constant and $y$ is constant. (Astra's bound $2mL$ is valid but crude; **the round-9 relax referee's "$m+1$ nonzero entries"
(`ref9L_relax.md` line 94) should read $(m+1)/2$**; the error was conservative and README_IXL §1.6's length $q/(16(m+1)P^4)$
remains valid.) Item 5: kernel vectors of $T=[\mathbf 1^t;S_M]$ with support $\le r+1$ and height $\le r^{r/2}$ (Cramer,
Hadamard); $S$ maps $\ker T$ onto the zero-sum vectors supported on $R$.

**Computed (referee; Astra's count reproduced).** Separations and constant-column counts at $m=7,11,19,23,43$ over all triples;
$m=11$, up to two missing rows: 67 row sets, 549 annihilators; $m=19$: 3,077; $m=11$ up to four missing rows: 3,684; sup-norm,
support and span checked (her script checks none of these).

**Comparison with IX.4 (upheld).** Stronger in hypothesis, weaker in length: IX.4 excluded *genuine* Paley configurations up to
$q^{1-O(1/m)}$; X.3 excludes *sign-free* tables up to $q^{5/8+o(1)}$ only. Once item 2 is applied, the rest is IX.4's column
argument verbatim. **New: item 2 (sign-free $\Rightarrow$ genuine on Paley cuts, any row subset) and item 5 for genuine
configurations.** Open: sign-free Paley tables at $q^{5/8+o(1)}\le L<q^{1-O(1/m)}$.

**The circle fakes' row class is not covered (paley: FALSE for the coverage question; meaning: GAP; upheld).** Those fakes use
$\{i:K_{i0}=1\}$ or a density-$\rho<1$ part of it, about $\rho m/2<\frac35m$ rows. Astra does not claim coverage (her X.5
"OPEN"). Heuristic (paley): no version of this divisor method reaches a retained fraction $\le\frac12$ for all large $m$, since
the divisor $q/\prod_{\rm supp}Q_v$ needs support $<m/2$ while generic kernel vectors of $T$ have support $\ge k+2$.
**Computed (meaning; reproduced [J], `judge10L_halfrow.out`, `judge10L_halfrow23.out`).** For the half-row pattern
($R=\{i:S_{i0}=+1\}$, $(m-1)/2$ rows): no kernel vector of $T$ of support $\le3$ at $m=7$ or $\le5$ at $m=11$ (so the argument
fails there); at $m=19$ the 36 kernel vectors of support 8 span $\ker T$ ($\dim 8$), so for two primes per column the half-row
pattern at $m=19$ admits no sign-free table and no genuine configuration once $P\ge P_0(C)$. At $m=23$ ($T$ is $13\times23$,
$\dim\ker T=10$) there is no kernel vector of support $\le10$, and the 341 of support 11 span $\ker T$; since $11<23/2$ the
divisor $q/\prod_{\rm supp}Q_v\ge P^{23}/(2P)^{22}$ beats the size for $P\ge P_0(C)$, so the half-row pattern at $m=23$ is excluded
as well (11 marks). Both agree with the referee. The margin is one column at $m=19$ and half a column at $m=23$; these are
fixed-$m$ statements and do not show a trend. The asymptotic question is open.

**Verdict.** Correct, with corrected quantifiers. Gap: none in the proofs. Withdrawn: none.

### 1.4 X.4, arbitrary cuts at the critical length (referee: recon; rechecked line by line by the judge)

**Proved, statement 1.** If $q$ (squarefree) has a prime $p>C^2$, then
$$B_{\rm int}(q,C,\delta)\le B_{\rm free}(q,C,\delta)\le2\max\{4,B_{\rm int}(q,C,\delta)\},$$
and at $C=1$, $B_{\rm free}\le\max\{4,B_{\rm int}\}$. Each side of the bit-partition at $p$ has $p\mid g$ on all its triples, so by
X.1 and X.2 a side of $\ge5$ marks is $|x_i-x_j|$ for integers $x$, which are then a genuine configuration (classes possibly
interchanged at some primes; with $B_{\rm int}$ optimised over the classes no alignment is needed).

**Proved, statement 2 (recon's explicit form, which I endorse).** Let $C>0$, $C_0=\max(1,C)$; let $q$ be squarefree with $s$ prime
factors, all in $[P,2P]$, $P>\max(16,C_0^8)$ (true for $s\ge s_0(C)$). Fix nonzero $\delta_p$ and any sign-free table on $B$ vertices.
One can delete at most
$$N\le2\Delta N_u+\Delta+1+\lfloor C^2\rfloor\le320\,s/B+97+\lfloor C^2\rfloor\qquad(\Delta=32,\ N_u<5s/B+1)$$
vertices (if $B\le2277$, delete everything) so that the survivors carry integers $y_i$ in an interval of length $C\sqrt q$, with
$y_i\equiv\rho_{p,\sigma_{ip}}\pmod p$, $\rho_{p,1}-\rho_{p,0}\equiv\pm\delta_p$ at every $p$, and $|y_i-y_j|=|d_{ij}|$ for **every**
surviving pair. ($\Delta=16$ is also valid, with threshold $B>629$.) Constants depend on $C$ only; nothing depends on $\epsilon$,
$A$, or the number of primes per column; dyadicity is used only in counting unbalanced primes (with unequal primes:
$\sum_{\rm unbalanced}\log p\le4W/B+8\log C_0$). At the unbalanced primes the surviving cut becomes constant.

*Judge's line-by-line check [J] of the three places Astra asked to have prioritised, and of the rest.*
* (a) $h_{ij}\mid d_{ij}\ne0$ gives $D_{ij}\ge W/2-\log C_0$; the unit vectors $v_i=(\pm\sqrt{\log p/W})_p$ have
  $v_i\cdot v_j=1-2D_{ij}/W\le\eta=2\log C_0/W$. Correct.
* (b) Far means $D_{ij}>5W/8$, i.e. $v_i\cdot v_j<-1/4$; with $S$ the sum of the $k$ far neighbours of $i$,
  $k^2/16<(v_i\cdot S)^2\le|S|^2\le k+\eta k(k-1)$, so $k\le16$ once $\eta\le1/256$. Correct.
* (c) In a triple $D_{12}+D_{13}+D_{23}=2(W-\log g)$, so three good edges give $\log g\ge W/16>2\log C_0$ and X.1 applies. Correct.
* (d) A colour class of a greedy $(\Delta+1)$-colouring of the far graph is a good clique $R$, $|R|\ge2\Delta+5$ for $B$ large;
  it embeds by X.2. Correct.
* (e) **Alignment on common anchors.** For $i\notin R$ the set $G_i$ of anchors good to $i$ has $\ge\Delta+5$ elements, $\{i\}\cup G_i$
  has only good edges, so it embeds; two isometric embeddings of a finite subset of $\mathbb Z$ with $\ge2$ points differ by a
  unique isometry of $\mathbb Z$ (a point of a line is fixed by its distances to two distinct points), so $y_i$ is the unique
  integer at distance $|d_{ia}|$ from every $a\in G_i$. No choice is made. Correct.
* (f) For a good pair $i,j\notin R$, $|A_{ij}|\ge|R|-2\Delta\ge5$ common good anchors; $\{i,j\}\cup A_{ij}$ embeds, and its
  coordinates for $i,j$ are forced to be $y_i,y_j$ by the distances to $A_{ij}\subset G_i\cap G_j$. Hence $|y_i-y_j|=|d_{ij}|$ on
  every good edge (including edges meeting $R$). Correct.
* (g) **Small bit classes.** A class of size $>2\Delta$ has a connected good graph (a component $X$ with $|X|\le|K|/2$ would give a
  member of $X$ more than $\Delta$ far neighbours); good same-class edges give $p\mid y_i-y_j$; so one residue per large class, and
  a good cross edge (which exists) gives $\rho_1-\rho_0\equiv\pm\delta_p$. Small classes are never used for residues; they are
  deleted in (h). Correct.
* (h) A prime with minority class $\le2\Delta$ has $|\sum_i\sigma'_{ip}|\ge B-4\Delta$; then
  $(B-4\Delta)^2\sum_{\rm unbal}\log p\le W|\sum v_i|^2\le BW+2B(B-1)\log C_0$, giving $N_u<5s/B+1$. Deleting minorities only
  restricts established properties. Correct.
* (i) Survivors have distinct $y$ (equal $y$ gives equal signature, hence $q\mid d\ne0$); deleting the $\le\Delta$ far neighbours of
  the smallest survivor puts all survivors in an interval of length $L$. Correct.
* (j) **Complementary pairs.** For survivors, $u=y_i-y_j$ and $d=d_{ij}$ satisfy $p^2\mid u^2-d^2$ at agreement primes and
  $p\mid u^2-d^2$ at separating (balanced) primes, so $qh_{ij}\mid u^2-d^2$ and $|u^2-d^2|\le C^2q$: $|u|=|d|$ whenever $h_{ij}>C^2$,
  i.e. whenever $h_{ij}>1$ (all primes exceed $C^2$). If an unbalanced prime exists there is no complementary pair. Otherwise every
  complementary pair has $y_i+y_j\equiv\sum$ of the two residues at every $p$, a single residue mod $q$ in a window of length
  $2L<q$, hence a single value; and $(y_i-y_j)^2$ in one residue class mod $q$ inside $[0,C^2q]$, hence $\le1+\lfloor C^2\rfloor$
  values; sum and $|$difference$|$ fix the unordered pair. So at most $1+\lfloor C^2\rfloor$ complementary pairs; delete one
  endpoint of each. Correct.

**Addition [J] (Proved, one line).** If the surviving signatures are Sidon, there is **at most one** complementary pair: two
complementary pairs $\{i,\bar i\}\ne\{j,\bar j\}$ give $\sigma_i+\sigma_{\bar i}=\mathbf 1=\sigma_j+\sigma_{\bar j}$. So for Sidon
patterns the last deletion step costs at most one mark.

**Referee additions (upheld).** (i) *Uniqueness* (Proved): if every prime of $q$ exceeds $C^2$, a non-complementary signature pair has
at most one admissible magnitude ($qh\mid d^2-d'^2$); sign freedom lives only on complementary pairs. (ii) *Necessity of deletion at
finite size* (Computed): at $q=59\cdot67\cdot71\cdot83$, $C=6$, Astra's $\delta$: $B_{\rm free}=8$ (the even-weight subgroup of
$\{0,1\}^4$) against $B_{\rm int}=5$.

**Computed [J]** (`judge10L_recon.out`). The five-mark example: $\mathcal H(1,4,5)=-1152q^2$, not line-embeddable, embeddable after
deleting mark 1 or mark 4. At $q=59\cdot67\cdot71\cdot83$: $B_{\rm free}=8$, $B_{\rm int}=5$ ($C^2=36$), $6$ ($C^2=49,58$); the only
multi-magnitude pattern is the complementary one, $\{73,23645\}$ ($C^2=36$), $\{73,23645,33439\}$ ($C^2=49,58$). On 180 random
moduli (4-6 primes, all $>C^2$, $C\in\{1,2,3\}$, best of 30 random $\delta$): 0 violations of statement 1, 0 violations of uniqueness,
and in every genuine configuration found the complementary pairs have one common sum and pairwise distinct squared differences
($\le1+\lfloor C^2\rfloor$ of them); the largest $B_{\rm free}-B_{\rm int}$ was 1 ($B_{\rm free}\le4$).
*Attempt to break the final step.* On the 8-clique, for 1,200 random choices of the complementary magnitudes, deleting one endpoint
of every complementary pair left 4 marks that were **never** line-embeddable: the four survivors (e.g. $0000,0011,0101,0110$) have
unique magnitudes and form a pseudo-linear quadruple. This is **not** a counterexample: X.4 uses $\ge2\Delta+5$ anchors to define
global coordinates before step (j), and the four-point exception is excluded by X.2 once five marks are present. It shows that the
anchor step cannot be dropped at small size and that random moduli with random $\delta$ are a weak test ($B_{\rm free}\le4$); the
thresholds of statement 2 cannot be reached by brute force. (Astra's and the referee's examples are all non-Sidon.)

**Verdict.** Correct. Gap: none (anchors, small classes and complementary pairs all checked). Withdrawn: none.

### 1.5 X.5, portability and exact tests (referee: port)

**Proved.** Let $q>1$ be squarefree with every prime factor $\equiv1\pmod4$ (implicit in X.5 through IX.0), and $z_1,z_2,z_3$ three
distinct marks of norm $q$ in one unit class, chords $a,b,c$. Then $a^2=4h_{12}t_{12}^2$ (etc.) and
$$\mathcal H(a,b,c)=16\,\mathrm{area}^2=(abc)^2/q=64g^2(t_{12}t_{23}t_{13})^2,\qquad \mathrm{area}=2g|t_{12}t_{23}t_{13}|.$$
In particular $\mathcal H>0$: circle triangles never degenerate.

**Proved (referee addition; new to the record).** At a split $p=\pi\bar\pi$, for $\pi\mid z_i$, $\bar\pi\mid z_j$:
$|z_i-z_j|^2\equiv-\psi_p(z_i)\varphi_p(z_j)\pmod p$, with $\varphi_p,\psi_p$ the reductions mod $\pi,\bar\pi$. The separated squared
chords mod $p$ form a rank-one matrix with free nonzero vertex labels; on a single triangle they impose nothing. The twin's
congruence is the case of equal labels. Hence X.1 has no circle counterpart beyond $\mathcal H=64g^2T^2$, which gives
$\max$ chord $\ge2g^{1/3}q^{1/6}$ (the Cilleruelo–Córdoba area argument sharpened by the Gaussian gcd), not the critical $q^{1/4}$.

**Computed (reproduced by port).** The norm-65 example (squared chords $80,52,4$; $52\equiv2$, $4\equiv4\pmod5$); the census (182, 252,
13; the 13 nondegenerate and the 13 capped triples are **disjoint** sets, the equal counts a coincidence); the $m=43$ relaxation
(12,341/12,341 triples fail; since $|\mathcal H|<1.2\cdot10^{-5}q^2$ "fails $q^2g^2\mid\mathcal H$" is the same as $\mathcal H\ne0$);
the weighted PTE families (55,026 triples; **tautological** for genuine configurations at any length, and none at critical length:
$\log q/\log L\in[1.07,1.88]$); the order-11 erasure count 549; the four-point and five-mark examples.

**Corrections (upheld).** "The fixed-modulus argument gives no circle rate" is an assessment of this argument, not a PROVED
statement. The five-mark example shows deletion can be necessary only for **non-Sidon** patterns (its signatures satisfy
$s_1+s_4=s_2+s_3$ and marks 1-4 are the four-point exception); whether deletion is ever needed for Sidon patterns is open
(with [J]'s addition in §1.4, at most $O_C(1+s/B)$ deletions of which at most one is complementary). Reproducibility nit:
`round31/astra_lemma_X/x_checks.py` imports `ix_checks`, present only in `scratchpad/astra/lemma/`.

**Verdict.** Correct; two overstatements of scope. Withdrawn: none.

### 1.6 Meaning (referee: meaning; item (6) rechecked by the judge)

**Upheld.**
* *Options:* (iv) delivered; (iii) $\Rightarrow$ (ii) delivered as a reduction (already with loss 2 from statement 1 of X.4, which
  needs only one prime $p>C^2$, X.1 and X.2; the anchor reconstruction is needed only for loss $O_C(1+s/B)$); (i), (ii) not delivered.
* *Reformulation (Proved, elementary):* a sign-free table is a choice, for each pair separately, of an integer $d$ with $|d|\le L$
  and $d^2\equiv\Delta_{ij}\pmod q$. Part X shows these pairwise short square roots interlock completely: they determine a genuine
  configuration up to a factor 2 in $B$, or up to $O_C(1+s/B)$ marks on dyadic parents.
* *"No known relaxation strictly between valuation data and the genuine problem":* fair as "none among the relaxations studied in
  the record"; not fair as "none exists". Keeping the squares only on a sub-product $q_E$ (valuations elsewhere) gives only
  $q_E^2g^2\mid\mathcal H$ against $|\mathcal H|\le81C^4q^2$ and is unstudied (Heuristic: probably as weak as IX.4).
* *Half-row Paley:* open asymptotically; excluded at fixed $m=19,23$ for large $P$ (reproduced [J], §1.3); the argument does not apply at $m=7,11$.
* *Circles (items 4, 5):* the twin mechanism does not transfer (circle triangle divisibility $g^2q$ against size $\asymp C^4q^2$; the
  Paley character $\alpha/\bar\alpha$, $\alpha=\Pi_v^{(m+1)/2}\bar\Pi_w^{(m+1)/2}$, has $\log|\alpha|=\frac{m+1}4\log(Q_vQ_w)\approx W/2$
  against the gap threshold $W/4$, failing by a factor 2 on the log scale). "Pair-level circle data do not force genuineness" is
  Heuristic. Part X says nothing either way about an integral-cofactor fake on circles.

**Item (6), the column-orthogonality route: partly overturned [J].** The referee's claims, one by one.
* *Basic lemma* (if $\sum a_i=0$ and the class sums of $a$ vanish mod $p$ for $p\in G$, then $\prod_Gp\mid a^ty$ for all $y\in M_S$;
  with $\|a\|_1L<2\prod_Gp$, $a^ty=0$ on spread $\le L$). **Proved**, trivial.
* *Siegel count* ("column orthogonality yields any relation only when $B>\frac12(1+\frac{\log2}{\log P})s+O_\epsilon(1)$; at
  $P=s^{1+o(1)}$ that regime is already empty by the larger sieve"). **Correct for EXACT orthogonality ($a^tS_p=0$) only; false for
  the congruence certificates of his own basic lemma.** [J, Proved] Let $\Lambda_S=\{a:\sum a_i=0,\ w\cdot a\equiv0\pmod q\}$ with
  $w_i\equiv$ (bit of row $i$ at $p$) mod $p$; it has index $\le q$ in the rank-$(B-1)$ zero-sum lattice, so Minkowski gives $0\ne a\in\Lambda_S$
  with $\|a\|_\infty\le(\sqrt Bq)^{1/(B-1)}$, and $\|a\|_1L/2<q$ as soon as $q^{1/2-1/(B-1)}\gg_{B} C$, i.e. for every pattern with
  $B\ge4$ and $q$ large. So a nontrivial linear relation exists for **every** cut pattern with $B\ge4$, far below $s/2$.
  Computed [J] (`judge10L_lattice.out`, LLL): random patterns with $B<s/2$ ($s=14,16,20$; $B=5,\dots,9$; 24 trials) always have
  $B-1$ independent certificates, so the linear method *succeeds* on them (they carry no nonconstant short vector).
* *Proposition M* ($\dim V_S\le2$: all vectors of $M_S$ of spread $\le L$ lie in $\mathbb Z\mathbf 1+V_S$, for $s\ge s_0(C)$, whenever
  $S$ carries a genuine configuration). **Proved [J, rechecked]**: $[\mathbb Z^B:M_S]=q^{B-2}$, the projection $N$ to
  $\mathbf 1^\perp$ has covolume $q^{B-2}/\sqrt B$, $\lambda_i(N)\le q$, $\lambda_1..\lambda_k\le\sqrt BL/2$, and Minkowski's second theorem
  gives $(k-2)W\le k\log B+B\log B+O(B\log C_0)$; with the record's sieve bound $B\log B\le(\frac12+o(1))W$ this forces $k\le2$
  (with constant 1 in place of $\frac12$ one would get $k\le3$). Simplest form: $q^{k/2}q^{n-k}\gtrsim q^{n-1}$. Computed: $k=1$ on the
  five-mark example (meaning; and [J]: certificate rank $3=B-2$). **Correction:** the referee's gloss "linear relations can give at most
  $B-3$ independent annihilators" has the inequality the wrong way round; the certificates span a space of dimension $B-1-k\ge B-3$
  (up to the factor $B^{O(1)}$ lost in transference), and at most $B-2$ for a pattern carrying a genuine configuration.
* *Shadow obstruction* (patterns with $2^{(1/2-o(1))s}$ rows whose $M_S$ contains a nonconstant two-valued $z$ of spread $\le\sqrt q$, so
  certificates never span). **Proved, but irrelevant to the rate [J].** In a shadow pattern the rows with equal $z$ agree at every prime
  outside $D$, so their agreement product is $\ge q/d\ge\sqrt q$ and genuine marks on one side are spaced by $\ge\sqrt q$: **a shadow
  pattern carries at most $2\lfloor C\rfloor+2$ genuine marks** and violates the necessary pair-valuation condition $h_{ij}\le C\sqrt q$
  for almost all its pairs. So the obstruction kills only "column orthogonality applied to arbitrary $S$, without the pair valuations",
  not the route "pair valuations, then column orthogonality". Computed [J]: shadow patterns ($s=14,16$) have certificate rank exactly
  $B-2$, as the referee says.
* *Transference statement* ("the linear method succeeds on $S$ iff, up to $B^{O(1)}$ in $L$, $M_S$ has no nonconstant short vector;
  that is the problem with nondegeneracy dropped"). **Upheld**, and made precise [J, Proved]: if $z\in M_S$ is nonconstant of spread
  $\le L$, let $E$ be the set of primes at which the two class residues of $z$ differ and $q'=q/q_E$; then $z=c\mathbf 1+q'k$ with $k$
  a two-class configuration (repetitions allowed only between rows with equal bits on $E$) for the modulus $q_E$ with the cut
  $S|_E$, in an interval of length $C\sqrt{q_E}/\sqrt{q'}$. So the linear method fails on $S$ exactly when some sub-pattern $S|_E$
  carries such a (possibly repeated) configuration at relative constant $C/\sqrt{q'}$; for $q'=1$ this is a genuine configuration.
  Excluding these is a problem of the same type as the twin itself.
* *"No gain beyond the constant $\frac12$ from this route" (FALSE verdict on the route).* **Not established; downgraded to
  Heuristic.** What is proved: Proposition M; the shadow obstruction for patterns carrying $O(C)$ genuine marks; the reformulation above.
  What is not proved: that on patterns with $B\ge cs$ rows satisfying the pair-valuation conditions $M_S$ has nonconstant short
  vectors. My reading: the route is not closed, but it is no easier than the problem (by the reformulation), and nothing in the record
  suggests a way to exploit it. Whether any use of the additive structure on a fixed modulus beats the larger-sieve constant $\frac12$
  for arbitrary cuts is exactly the open rate question; nothing in Part X or in the reports decides it.

**Recommendation (stop the thread):** upheld (§5).

## 2. What is new relative to the record and correct

New to the record and correct:
1. **X.1, the sign-free triangle obstruction**: $q^2g^2\mid\mathcal H$ for any three magnitudes satisfying the separated-square
   congruences, with the unit constant $\max\ge\sqrt{qg}$ for nondegenerate triples (asymptotically sharp even at $g=1$ [J]).
   Elementary; new in this context.
2. **The collapse of the sign-free relaxation** (X.4): $B_{\rm int}\le B_{\rm free}\le2\max(4,B_{\rm int})$ when some $p>C^2$ divides $q$
   (no loss at $C=1$); on dyadic parents, $B_{\rm int}\ge B_{\rm free}-O(1+C^2+s/B_{\rm free})$. Pair-level squares mod $q$ plus magnitudes
   determine a genuine configuration; sign coherence is a consequence, not an input. With the referee's uniqueness (magnitudes forced
   off complementary pairs) and [J]'s "at most one complementary pair for Sidon signatures".
3. **Sign-free $\Rightarrow$ genuine on Paley cuts on any $\ge5$ rows** below $\sqrt{qg_{\min}}=q^{5/8+o(1)}$, and the exclusion of
   genuine (hence sign-free) Paley configurations on every $\ge\frac{t+1}{2t+1}m$ rows at critical length (X.3, items 2 and 5). The
   full-core exclusion for genuine configurations is IX.4's.
4. **The circle rank-one label congruence** (port) and the exact explanation of non-transfer (port, meaning).
5. **Proposition M** (meaning; probably folklore) and the reformulation of the linear method (meaning, made precise [J]).

Not new: the line lemma X.2 (Menger 1928); the column argument of X.3 (IX.4); the circle area identity as a bound (Cilleruelo–Córdoba,
refined by IX.0's Gaussian gcd). Literature: no full text was read in this round.

## 3. Rulings

**(a) Which option was delivered, and what is now proved.** **(iv)**, plus the reduction **(iii) $\Rightarrow$ (ii)**; **not (i)** (no
$o(s)$ bound for any explicit dense family or for every $q$) and **not (ii)** (no genuine configuration of order $s$). The pair
congruences that give the Paley contradiction are the separated-square congruences $d_{ij}^2\equiv\delta_p^2\pmod p$ at the
$(m-3)/4$ constant columns of each triple (forcing every triangle to degenerate, hence genuine line coordinates), followed by IX.4's
column orthogonality. For an arbitrary Sidon cut pattern with $B\gg s$ they constrain nothing beyond genuineness: they force the table
to be genuine after $O_C(1+s/B)$ deletions, and give no bound on $B$. Precisely, on a fixed squarefree modulus it is now proved that:
* every triple of a sign-free table satisfies $q^2g^2\mid\mathcal H$, and degenerates if its bits agree at some $p>C^2$ (all triples at $C=1$);
* $B_{\rm int}\le B_{\rm free}\le2\max(4,B_{\rm int})$ if some $p>C^2$ divides $q$, $B_{\rm free}\le\max(4,B_{\rm int})$ at $C=1$;
* on dyadic parents with $s\ge s_0(C)$, every sign-free table on $B$ rows is genuine after at most $320s/B+97+\lfloor C^2\rfloor$
  deletions, the surviving cut possibly constant at $<5s/B+1$ primes;
* on Paley cuts (any $t$, any assignment, any $P$): sign-free tables on $\ge5$ rows below $\sqrt{qg_{\min}}$ are genuine; none on the
  full core at $C\sqrt q$ ($m\ge m_0(C,t)$) or at $q^{5/8-\epsilon}$ (equal $t$); none at $C\sqrt q$ on $\ge(\frac{t+1}{2t+1}+\epsilon)m$
  rows (equal $t$); the half-row pattern is excluded for large $P$ at $m=19,23$ (computed) and open in general.

**(b) What it means.** Among the relaxations studied in the record **there is now nothing between the valuation data and the genuine
twin problem**: the IX.4 relaxation (magnitudes, gcds, valuation patterns, Sidon signatures, integral cocycle) reaches the larger-sieve
constant $\frac12$, and the next natural relaxation, which adds the pair squares mod $q$, is already the genuine problem. So a proof of
the rate may use only pairwise short-square-root data ($|d_{ij}|\le C\sqrt q$, $d_{ij}^2\equiv\Delta_{ij}\pmod q$) and gains nothing
from sign coherence; conversely no relaxation-type barrier can be built from pair-level residues. The bottleneck is $\binom B2$
interlocking short-square-root conditions with one common $\delta$: Conjecture 9 (rate form) restated pairwise. Unstudied
intermediates exist (squares on a sub-product $q_E$); they do not bear on the rate.

**(c) Transfer to circles: none.** On circles the separated squared chords mod $p$ are a rank-one product of free vertex labels, so
the twin's step "two separated squares agree, so two Heron factors vanish" has no counterpart; circle triangles are never degenerate;
the only Heron identity is $\mathcal H=64g^2T^2$ (area $=2g|T|$), which controls arcs of length $\asymp q^{1/6}$, not $q^{1/4}$; the
circle triangle divisibility is $g^2q$, a full factor $q$ short; and the Paley character has $\log|\alpha|\approx W/2$ against the
threshold $W/4$. The twin and circle edges also live on different scales relative to the modulus ($d^2\le C^2q$ against
chord$^2\le C^2\sqrt q$). Part X is evidence neither for nor against an integral-cofactor fake on circles.

**(d) Consequences. Blunt: nothing moves.**
* Worst-case rate in the twin ($o(\omega(q))$, the rate form of Conjecture 9): **unchanged**; no bound for any explicit dense family, no
  order-$s$ example. What changed: the problem is shown to be insensitive to the sign data (sign-free $=$ genuine), and the Paley
  profile is excluded on a larger row range.
* Worst-case growth rate on circles, $(2/3+\epsilon)\log R/\log\log R$ (Theorem 1.1): **unchanged**; the mechanism does not transfer.
* Uniform bound, Questions A and B: **unchanged**.
* Column orthogonality as a route to the twin rate is **not closed** (contrary to the meaning referee) but is equivalent to excluding
  short degenerate configurations on sub-patterns, a problem of the same kind; Heuristic: no gain beyond $\frac12$ from it.

## 4. Recommended edits to the expert note

1. **§5 item 12, header (lines 929-934).** Add `README_XL.md` to the list of judge reports; add five referee blocks to the count and the
   exchange to the number of exchanges (coordinate with the parallel thread's judge, who is adding the round-10 circle report); add
   2/10 to the list of scores.
2. **§5 item 12, new bullet** after the bullet "***Part IX of the Question B thread: the critical two-class problem***" (ends line 1395),
   before "*Withdrawn or corrected in the round:*":
   > * ***Part X of the Question B thread: signed residues on one fixed modulus*** (`round31/README_XL.md`; five referee blocks;
   >   nothing false; judge's score 2/10). In the critical two-class problem ($q$ squarefree, interval $C\sqrt q$, classes at each
   >   $p\mid q$ differing by $\delta_p$), a *sign-free table* prescribes for each pair of marks an integer $d_{ij}$ with
   >   $0<|d_{ij}|\le C\sqrt q$, $\gcd(d_{ij},q)=h_{ij}$ and $d_{ij}\equiv\pm\delta_p\pmod p$ at every separating $p$, the sign free per pair
   >   and per prime. **Proved:** on every triple $q^2g^2$ divides the Heron product ($g$ the product of the primes constant on the
   >   triple), so a nondegenerate triple has a side $\ge\sqrt{qg}$ (sharp, also at $g=1$); with Menger's line-embedding theorem this
   >   gives $B_{\rm int}\le B_{\rm free}\le2\max(4,B_{\rm int})$ whenever $q$ has a prime $p>C^2$ (no loss at $C=1$), and on parents with
   >   all primes in $[P,2P]$ every sign-free table of $B$ rows becomes a genuine configuration on the same modulus and interval after
   >   deleting $O_C(1+s/B)$ rows. So the pair-level square congruences already force genuineness: among the relaxations studied in the
   >   record nothing lies strictly between the valuation data of Part IX (larger-sieve constant $\frac12$) and the genuine problem.
   >   No sign-free table exists on the Paley cut pattern of Part IX for any modulus and prime assignment: on the full core at length
   >   $C\sqrt q$ and, with equal numbers of primes per column, at $q^{5/8-\epsilon}$; and at critical length on every set of at least
   >   $(3/5+\epsilon)m$ rows (two primes per column) or $(2/3+\epsilon)m$ rows (one prime). For genuine configurations only the row-subset
   >   statement is new (the full core was excluded up to $q^{1-O(1/m)}$ in Part IX). **Open:** a worst-case $o(\omega(q))$ bound, a genuine
   >   example of order $\omega(q)$, the half-row Paley pattern (the row class of the circle fakes; excluded computationally at $m=19,23$
   >   for large primes). The mechanism does not transfer to circles: separated squared chords modulo a split prime are a rank-one product of
   >   free vertex labels, circle triangles never degenerate, and the only Heron identity is $\mathrm{area}=2g|t_{12}t_{23}t_{13}|$.
   >   Linear (column-orthogonality) certificates succeed on a cut pattern exactly when its lattice has no nonconstant short vector,
   >   which is a problem of the same type as the twin itself. **The growth rate for general $R$, the uniform bound and Questions A and B
   >   are unchanged**; the two-class thread is closed. Do not write "Part X proves Conjecture 9 for Paley moduli", "the sign-free
   >   relaxation is a new barrier", "there is no relaxation between valuation data and the genuine problem" (say "none among those
   >   studied"), "Paley profiles are impossible in the twin" (without "full core, or at least $3/5$ of the rows"), "the circle fakes are
   >   impossible in the twin", "(iii) and (ii) are equivalent" (loss of a factor 2 or $O_C(1+s/B)$ marks), "column orthogonality is
   >   closed as a route" or "column orthogonality gives a rate beyond the sieve constant", "the twin result implies that pair-level
   >   circle data force genuineness" (or the converse), or "Part X gives $B=o(s)$ on some dense modulus".
3. **§5 item 12, "*Withdrawn or corrected in the round:*" (line 1396).** Append: "in Part X, the line lemma is Menger's theorem (1928;
   citation added); the coefficient '$m+1$ nonzero entries' of a round-9 referee for $S_v-S_w$ (it has $(m+1)/2$; the derived length
   bound stays valid); and a referee's claim that column orthogonality is closed as a route (only shown for patterns carrying
   $O(C)$ genuine marks)."
4. **§6 Pointers (line 1437).** Add `README_XL.md` to the round-31 row and the phrase "sign-free residue tables collapse onto the
   genuine twin; Paley cuts excluded on $\ge3/5$ of the rows".

Sentences the note must not contain: those listed in edit 2.

## 5. Should the thread stop?

**Yes. Stop the Question B two-class thread.** The previous judge's rule was: stop if the answer is a further random-parent variant or a
further restricted class. The answer is a further restricted class (the Paley pattern on larger row sets) and a structural collapse
(sign-free $=$ genuine). The collapse is genuinely informative, because it closes the relaxation ladder from below, but it does not touch
the worst-case rate, and every remaining question on this thread is either the rate itself (Conjecture 9, rate form) or a fixed-pattern
exclusion.

If, against this advice, one closing exchange is wanted, the one question is the meaning referee's: *on the Paley core of IX.4, let $T$
be the all-ones row together with the rows $i$ with $S_{i0}\ne+1$; prove or disprove that for all large $m\equiv3\pmod4$ the kernel of
$T$ is spanned by integer vectors of support $\le(\frac12-c)m$ and height $m^{O(1)}$* (equivalently, decide the half-row Paley pattern,
the row class of the circle fakes, in the twin). Honest probabilities: of moving the worst-case rate, $<0.005$; of a clear answer, about
$0.3$. Its only value is to close the comparison with the circle fakes. Effort is better spent on the parallel circle thread.

## 6. Scores (0-10, on the goal)

| block | score | reason |
|---|---|---|
| heron (X.1, X.2) | 2 | correct and checked [J]; sharp elementary triangle obstruction; the line lemma is Menger's |
| paley (X.3) | 1.5 | correct with corrected quantifiers; the new content is sign-free $\Rightarrow$ genuine and the row-subset range; the full core was IX.4's |
| recon (X.4) | 3 | correct, checked line by line [J]; the collapse of the last relaxation onto the genuine problem at scale $s$; no bound on $B$ |
| port (X.5) | 1 | correct; explains the non-transfer to circles precisely; regressions are consistency checks |
| meaning (referee) | 1.5 | sound reading and good wording; Proposition M correct; its "route closed" verdict overturned in the range that matters |
| **Part X overall** | **2** | **nothing false; option (iv) and the collapse sign-free $=$ genuine; no movement on the worst-case rate in either problem, the uniform bound or Questions A and B; thread closed** |

My own checks: X.1 (divisibility, quotient mod 4, unit constant, sharpness at $g=1$; §1.1); X.4 line by line, with the
uniqueness, complementary-pair and statement-1 tests and an attempt to break the final step (§1.4); the meaning referee's item (6)
(congruence certificates, Proposition M, the shadow patterns, the reformulation; §1.6); the half-row Paley circuits (§1.3).
