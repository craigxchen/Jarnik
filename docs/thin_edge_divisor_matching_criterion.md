# Direct edge divisors require a matching before a thin Hadamard bridge is possible

Reanchored rows and the standard pair quotients are indexed by edges
between the original points. Four such words cannot yield a common-support
Hadamard quartet by Gaussian division if two of their edges share a
vertex, unless a changing group is paid for by the correcting factors.
The obstruction is primewise and retains correction overlap depth.

Thus this direct-divisor route must start with four vertex-disjoint
edges, requiring at least eight points. At eight points the incidence
patterns do allow the desired four columns, but the resulting divisors
are not automatically thin. Additive constructions and large ordinary
real factors are outside this incidence obstruction.

## 1. Edge words and the exact correction bound

Label the points by \(V=\{0,1,\ldots,m\}\), with distinguished anchor
\(0\). For every nonempty \(T\subseteq\{1,\ldots,m\}\), let
\(H_T\) be a conjugate-primitive Gaussian integer. Assume their norms
are pairwise coprime. Put \(\chi_T(0)=0\), with \(\chi_T\) the
indicator elsewhere. For an oriented edge \(e=(a,b)\), define

\[
E_e=\prod_{\chi_T(a)-\chi_T(b)=1}H_T
     \prod_{\chi_T(a)-\chi_T(b)=-1}\overline{H_T}.
                                                               \tag{1}
\]

These are exactly the core words of the anchor rows and their primitive
pair quotients, up to edge orientation. Take four sources

\[
W_j=K_jE_{e_j},\qquad 0\ne K_j\in\mathbb Z[i],
\quad j=0,1,2,3.                                      \tag{2}
\]

The \(K_j\) may overlap any core prime to any depth. Let nonzero
Gaussian integers \(U_j\) be direct divisors of these sources:
\(W_j=D_jU_j\), with \(D_j\in\mathbb Z[i]\). Suppose

\[
\begin{aligned}
U_0&=J_0DABC,& U_1&=J_1DA\overline B\overline C,\\
U_2&=J_2D\overline A B\overline C,&
U_3&=J_3D\overline A\overline B C,
\end{aligned}                                         \tag{3}
\]

where \(A,B,C\) are nonunit conjugate-primitive Gaussian integers,
\(D\ne0\) is conjugate-primitive and may be a unit, their norms are
pairwise coprime, and the \(J_j\) are arbitrary nonzero Gaussian
integers. No thinness assumption is needed for the following result.

**Shared-vertex obstruction.** If two source edges \(e_r,e_s\) share
a vertex, then at least one \(G\in\{A,B,C\}\) satisfies

\[
N(G)\mid N(K_r)N(K_s).                                \tag{4}
\]

To prove this, let \(v\) be a shared vertex, and define
\(\epsilon_e(v)=1\) if \(v\) is the first endpoint of the oriented
edge, and \(-1\) otherwise. Whenever a cut separates both edges,

\[
\chi_T(a)-\chi_T(b)
 =\epsilon_e(v)(2\chi_T(v)-1).
\]

The two core orientations therefore have the fixed ratio

\[
\eta=\epsilon_{e_r}(v)\epsilon_{e_s}(v)\in\{1,-1\}.   \tag{5}
\]

For any pair of different rows in (3), the three changing columns
have sign products consisting of one \(+1\) and two \(-1\).
Choose a group \(G\) whose product is \(-\eta\).

Fix a split rational prime \(p\) in \(N(G)\), and let
\(\pi^d\) be its full oriented power in \(G\). The core norms are
pairwise coprime, so at most one \(H_T\) contains this rational prime.
If that cut does not separate both source edges, at least one source
has no core factor at this prime. If it does separate both, (5) and
the chosen opposite target sign product imply that at least one
source lacks the oriented prime power required by its target in (3).
In either case the full depth \(d\) must occur in the corresponding
source correction \(K_r\) or \(K_s\). Therefore

\[
d\le v_p\bigl(N(K_r)N(K_s)\bigr).
\]

Multiplying these inequalities over the norm primes of \(G\) proves
(4). Primes absent from all cores are handled by the same correction
argument. This proof does not discard a whole core prime power when a
correction contains a smaller power of that prime. Extra factors in
the \(J_j\) cannot remove the required target divisibilities.

Repeated unoriented edges are covered: their two core orientation
rows are identical or opposite according to the chosen orientations,
so the same shared-vertex argument applies.

## 2. Consequence for extracted profiles and rational multipliers

For a balanced extracted profile, suppose

\[
\log N(H_T)=w+o(w),\qquad \log N(K_j)=o(w).
\]

If all three changing groups in (3) have logarithmic norm bounded
below by a positive multiple of \(w\), equation (4) rules out every
pair of source edges with a common vertex. Hence the four sources
must be a matching of four vertex-disjoint edges. In particular, no
such direct-divisor construction is possible from a five-, six-, or
seven-point full profile.

The groups in (3) need not be whole \(H_T\)'s: splitting a composite
core does not change its cut orientation column, and the primewise
proof already permits such splitting.

Equation (4) is stated for Gaussian-integer corrections and direct
Gaussian divisibility. If a source or output is additionally multiplied
by a Gaussian rational number of subpower numerator and denominator
norm, clear those denominators first. An equation

\[
U=\frac{a}{bD}E_e\quad(a,b,D\in\mathbb Z[i])
\]

implies \(U\mid aE_e\); the cleared numerator plays the role of
\(K\) in (4). Thus subpower rational corrections do not evade the
asymptotic consequence. Their numerator heights must be included in
the budget; arbitrary large rational or ordinary integer multipliers
are not free.

In particular, the assertion is about the edge words (1), with their
controlled corrections. A large ordinary real integer may contain
both Gaussian orientations and is not an edge core word. This note
does not exclude constructions beginning with such an integer, or
with a nontrivial additive combination of edge words.

## 3. At eight vertices the incidence obstruction disappears

For eight vertices choose oriented edges

\[
(1,0),\quad(3,2),\quad(5,4),\quad(7,6).
\]

There are eight cuts, with the anchor side fixed, that separate all
four edges. Their orientation columns give all eight sign lines in
\(\{1,-1\}^4/\{\pm1\}\). Four explicit cuts are

\[
\begin{array}{c|c|c}
\text{group}&T&\text{column across the four edges}\\ \hline
D&\{1,3,5,7\}&(+,+,+,+)\\
A&\{1,3,4,6\}&(+,+,-,-)\\
B&\{1,2,5,6\}&(+,-,+,-)\\
C&\{1,2,4,7\}&(+,-,-,+).
\end{array}                                           \tag{6}
\]

Taking the corresponding \(H_T\)'s as \(D,A,B,C\) gives exactly
the four incidence products in (3), each a Gaussian divisor of its
source core word. For more vertices the same cuts can be extended.
Thus four disjoint edges are sufficient at the level of incidence.

They are not sufficient for thinness. Dividing the near-real source
by its complementary Gaussian factor can rotate the quotient by an
uncontrolled phase. In fact
[the matching-norm gap](four_signed_witness_matching_norm_gap.md)
already proves that these four quotient cores cannot all have
subpower imaginary parts after subpower corrections when the changing
groups grow. A full-circle bridge would have to derive that extra
thinness from the original simultaneous arithmetic, rather than infer
it from common divisibility.

## 4. A separate half-height restriction on thin Gaussian division

There is a simple exact limitation when both the numerator and the
quotient are known to be thin. Suppose \(V=DU\ne0\) in
\(\mathbb Z[i]\), and put
\(Y_V=|\operatorname{Im}V|\), \(Y_U=|\operatorname{Im}U|\).
If \(D\) is nonreal, then its imaginary part is a nonzero integer.
From \(D=V\overline U/|U|^2\) one obtains

\[
1\le|\operatorname{Im}D|
\le\frac{Y_V}{|U|}+\frac{|V|Y_U}{|U|^2},
\qquad
|U|^2\le Y_V|U|+|V|Y_U.                               \tag{7}
\]

If \(|V|\le\exp(\alpha w+o(w))\) and both imaginary heights are
subpower, a nonreal divisor therefore forces
\(|U|\le\exp(\alpha w/2+o(w))\). An output above that scale can
only have a real Gaussian-integer divisor. This applies, in particular,
to a subpower real-integer linear combination of thin actual rows:
its imaginary height stays subpower. Complex coefficients require a
separate cancellation argument and are not automatically covered.

For the four-row balanced profile, row moduli are
\(\exp(4w+o(w))\). Such a nonreal division cannot leave a thin
output larger than \(\exp(2w+o(w))\). A four-group core of one
balanced block per group lies exactly at that critical scale.

The threshold is sharp even with conjugate-primitivity and coprime
factor norms. For every positive even integer \(a\), set

\[
D=a+2+i,\qquad U=a-i,\qquad
V=DU=(a+1)^2-2i.                                      \tag{8}
\]

All three Gaussian integers are conjugate-primitive. The two factor
norms are odd, and their gcd divides \(4(a+1)\). An odd common
prime would then divide \(a+1\) and \(a^2+1\equiv2\), impossible.
Thus \(\gcd(N(D),N(U))=1\). Both factors have imaginary height one
and modulus comparable to \(\sqrt{|V|}\), while \(Y_V=2\).
This is a sharp individual factorization example, not a full Boolean
endpoint profile or a compatible four-witness construction.

The graph obstruction and (7) delimit different operations. They do
not supply a positive residual-height exponent for arbitrary actual
profiles, and the general radius-uniform circle theorem remains open.

## 5. Verification

The [exact checker](check_thin_edge_divisor_matching_criterion.py) exhausts
28,035 four-edge sets in `K_5,...,K_8`: the Hadamard incidence is feasible
exactly for the 105 matchings in `K_8`. It also checks the eight explicit
sign lines, 1,000 positive-even-a sharp factorizations with their
primitivity and coprime norms, and 144 shared-vertex sign/depth cases.
Root, Astra and Sol audited the prose argument; the finite checks do
not replace its all-dimension proof or establish thinness of the outputs.
