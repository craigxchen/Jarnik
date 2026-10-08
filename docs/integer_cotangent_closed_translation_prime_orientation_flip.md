# Closed cotangent translations can flip an actual cut-prime orientation

The integer cotangent condition has a local symmetry that a sign-only
reciprocity test must respect.  It acts on **actual integer cliques**, not
merely on formal Gaussian allocations.  The symmetry has an uncontrolled
global height cost; it does not give an endpoint counterexample or improve
the unproved all-edge height target in
[the integer formulation](integer_cotangent_lcm_height_target.md).

Let $L>0$, let the $X_i$ be distinct integers, and assume every

\[
 Q_{ij}=(X_iX_j+L^2)/(X_i-X_j)\in\mathbb Z.
\]

Put $\Delta_{ij}=X_i-X_j$ and
$D=\operatorname{lcm}_{i<j}|\Delta_{ij}|$.
For an integer translation $X_i'=X_i+t$, with the same $L$, the exact
update is

\[
 Q'_{ij}=Q_{ij}+\frac{t(X_i+X_j+t)}{\Delta_{ij}}. \tag{1}
\]

Thus the translated tuple is another integer cotangent clique exactly when

\[
 \Delta_{ij}\mid t(X_i+X_j+t)\qquad(i<j). \tag{2}
\]

Equation (2) is already in
[the whole-clique translation audit](integer_cotangent_translation_obstruction.md).
The refinement below uses its nonzero local solution.

## A clean collision prime

Let $p\equiv1\pmod4$ be odd, $p\nmid L$, and let $E$ be the
largest $p$-valuation of any $\Delta_{ij}$.  Suppose $E\ge1$ and there is a
single collision class $S$, $|S|\ge2$, with

\[
 X_i\equiv x\pmod{p^E}\quad(i\in S),
 \qquad p\nmid\Delta_{ij}\quad\text{if }\{i,j\}\not\subseteq S.
 \tag{3}
\]

The source divisibility for one inside pair gives
$x^2+L^2\equiv0\pmod{p^E}$: indeed $\Delta_{ij}$ divides
$X_i^2+L^2$, since
$X_i^2+L^2=\Delta_{ij}(Q_{ij}+X_i)$.
Choose $t$ by the Chinese remainder theorem so that

\[
 t\equiv-2x\pmod{p^E},\qquad
 t\equiv0\pmod{q^{E_q}}\quad
 (q\ne p, q^{E_q}\Vert D). \tag{4}
\]

For pairs inside $S$, $X_i+X_j+t\equiv0\pmod{p^E}$, so (2) holds
at $p$.  No other pair difference is divisible by $p$.  At every
other prime dividing $D$, $t\equiv0$ to the required depth, so (2)
holds there as well.  Hence **all** translated $Q'_{ij}$ are integers.

For each $i\in S$, $X_i+iL$ is divisible by one Gaussian prime
orientation above $p$.  After translation,
$X_i'\equiv-x\pmod p$, so $X_i'+iL$ is divisible by its conjugate
orientation and is a unit at the original one.  This is an exact
orientation flip at $p$.  To keep the cut support exactly $S$ modulo
$p$, also require every outside $X_j$ and $X_j-2x$ to be a
nonroot of $-L^2\pmod p$.  The construction does not assert that
valuation depths or the least radius are preserved.

The same proof works **simultaneously** for any finite collection of
clean collision primes.  At each chosen prime, independently prescribe
$t\equiv0$ or $t\equiv-2x$ at its maximal pair-difference power, and
prescribe zero at all other prime powers dividing $D$.  The Chinese
remainder theorem then gives one globally closed translation for every
choice of orientation flips.  Clean support after a flip requires the
outside nonroot condition above at each selected prime.  This is
independence of exact local signs, not a bound on the size of $t$.

## Two actual cut blocks

The following literal four-point clique has three finite cotangents:

\[
 L=1,\qquad (X_1,X_2,X_3)=(31,32,57),\qquad
 (Q_{12},Q_{13},Q_{23})=(-993,-68,-73). \tag{5}
\]

At $p=13$, the clean pair cut $\{1,3\}$ has residues
$(5,5)\mid6$; the two roots of $-1$ are $5,8$.  At $p=5$, the
clean pair cut $\{2,3\}$ has residues $(2,2)\mid1$; its roots
are $2,3$.  The pair-difference lcm is $D=650=2\cdot5^2\cdot13$.

Take $t=250$.  It satisfies $t\equiv-2\cdot5=3\pmod{13}$
and $t\equiv0\pmod{25}$; the factor two is also covered because
$t$ is even.  The translated actual clique is

\[
 (X'_1,X'_2,X'_3)=(281,282,307),\qquad
 (Q'_{12},Q'_{13},Q'_{23})=(-79243,-3318,-3463). \tag{6}
\]

The $13$-cut now has residues $(8,8)\mid9$, so its Gaussian
orientation flips while the outside row stays a nonroot.  The $5$-cut
retains residues $(2,2)\mid1$.  Both source and target all-edge least
squared radii retain $v_{13}(N)=1$ and $v_5(N)=3$, while the total
radii are

\[
 N=2{,}465{,}125,\qquad N'=455{,}260{,}346{,}125. \tag{7}
\]

Their exact normalized arc constants
$2N^{1/4}\arctan(1/\min X_i)$ are approximately $2.56$ and
$5.85$, respectively.  The growth in (7) is why this fixture is a
scope test, not an endpoint construction.

An endpoint-sized **source** fixture has two independently flippable
cuts.  Take $L=1$ and $X=(157,182,447)$.  All three pair cotangents
are integral: $(-1143,-242,-307)$.  The clean $29$-cut is
$\{1,3\}$ with residues $(12,12)\mid8$; flipping gives
$(17,17)\mid13$, and both outside residues are nonroots.  The clean
$53$-cut is $\{2,3\}$ with residues $(23,23)\mid51$;
flipping gives $(30,30)\mid5$, again with nonroot outside residues.
The pair-difference lcm is $D=76{,}850$.  Closed CRT translations
for the four flip choices are:

| Flip at 29 | Flip at 53 | $t$ | All-edge $N_t$ | $C_*(t)$ |
|---|---|---:|---:|---:|
| No | No | 0 | 212,298,125 | 1.538 |
| No | Yes | 65,250 | 102,866,995,511,848,031,667,125 | 17.317 |
| Yes | No | 29,150 | 842,320,937,321,019,957,125 | 11.626 |
| Yes | Yes | 17,550 | 41,549,115,737,611,768,325 | 9.068 |

All four configurations retain $v_{29}(N_t)=v_{53}(N_t)=1$.
This source itself has $C_*<2$, while every nonzero flip loses that
endpoint scale.  It makes the height cost of independent sign flips
explicit; it does not rule out a different, configuration-dependent
closed translation with better height.

Any proposed quadratic or quartic reciprocity restriction using only
the orientation signs of these two cut blocks must survive this flip.
The construction leaves open a restriction that also uses coefficient
heights, all other primes, or the endpoint arc condition.  See
[the checker](check_integer_cotangent_closed_translation_prime_orientation_flip.py)
for exact integer divisibility, both prime cuts, and the all-edge norms.
