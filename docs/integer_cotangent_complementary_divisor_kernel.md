# The complementary offset divisor kernel and its missing content factor

Complementary reanchoring pairs every positive offset divisor
$d_i=X_i-A$ with $c_i=(A^2+L^2)/d_i$.  The direct divisor count can
be made scale-invariant by removing the two common offset contents.
It gives a universal bound for any integer cotangent clique, but its
divisor parameter can exceed the all-edge squared radius through an
exact anchor-content factor.
This is a precise scope audit of the proposed direct analogue of the
conditional reciprocal-stretch divisor grid; it does not prove the
unproved exponent-four height target.

Let $A=X_0<X_1<\cdots<X_{m-1}$, with $m\ge2$, and assume

\[
 Q_{ij}=\frac{X_iX_j+L^2}{X_i-X_j}\in\mathbb Z
 \quad(0\le i<j<m),\qquad L>0. \tag{1}
\]

Put $S=A^2+L^2$, $d_i=X_i-A$, and $c_i=S/d_i$ for $i\ge1$.
The pair condition involving $A$ gives $d_i\mid S$ exactly.  Define

\[
 h_d=\gcd_{i\ge1}d_i,\qquad
 h_c=\gcd_{i\ge1}c_i,\qquad
 H=\frac{S}{h_dh_c}. \tag{2}
\]

Since $h_d\mid d_i$ and $h_c\mid c_i$ for every $i$,
$h_dh_c\mid d_ic_i=S$, so $H$ is a positive integer.  The two
primitive complementary factors satisfy

\[
 \frac{d_i}{h_d}\frac{c_i}{h_c}=H. \tag{3}
\]

The factors $d_i/h_d$ are distinct positive divisors of $H$.
Consequently the exact universal divisor bound is

\[
 \boxed{m-1\le\tau(H).} \tag{4}
\]

No multiplicity-two loss is needed: every original offset is distinct.

There is also an exact lcm form for the kernel:

\[
 H=\operatorname{lcm}_{i\ge1}(d_i/h_d)
  =\frac{\operatorname{lcm}_{i\ge1}(d_i)}{h_d}. \tag{4a}
\]

Indeed, for `D_i=d_i/h_d`, equation (3) gives `c_i/h_c=H/D_i`.
The latter integers have gcd one by the definition of `h_c`, while
`gcd_i(H/D_i)=H/lcm_i(D_i)`. Thus `lcm_i(D_i)=H`. This identity is
scale invariant and shows that the divisor parameter is exactly the
normalized lcm of the positive offsets, rather than merely a divisor of
`S`.

Under common integer scaling $(A,L,X_i)\mapsto u(A,L,X_i)$,
$S\mapsto u^2S$, while $d_i,c_i,h_d,h_c\mapsto u$ times themselves.
Thus $H$ and (4) are unchanged, as is the least squared radius $N$.

There is an exact comparison with $N$.  Put $g_A=\gcd(A,L)$,
$a=A/g_A$, $b=L/g_A$, and let $\epsilon_A=2$ if $a,b$ are both odd,
otherwise $\epsilon_A=1$.  The reduced anchor edge norm

\[
 n_A=\frac{a^2+b^2}{\epsilon_A}
\]

divides the all-edge least squared radius $N$.  Therefore

\[
 \boxed{H=\frac{g_A^2}{h_dh_c}\,\epsilon_A n_A,
 \qquad H\le 2N\frac{g_A^2}{h_dh_c}.} \tag{5}
\]

The ratio $g_A^2/(h_dh_c)$ is itself invariant under common scaling.
It is the exact missing content factor in an attempted replacement of
$H$ by a divisor parameter controlled by $N$.  Even when that ratio
is small, $\tau(N)$ is not radius-independent, so (4) alone does not
establish uniformity.

For a literal clique in the target positive range, take

\[
 L=70,\qquad X=(85,182,210),\qquad
 (Q_{01},Q_{02},Q_{12})=(-210,-182,-1540). \tag{6}
\]

Here $A=85\ge L$, $S=12{,}125$, $d=(97,125)$,
$c=(125,97)$, and $h_d=h_c=1$.  The reduced all-edge norms are

\[
 (485,97,5,5,97,485),\qquad N=485. \tag{7}
\]

The anchor content is $g_A=5$, $\epsilon_A=1$, and $n_A=485$,
so $H=25N$.  This is an actual integral cotangent clique showing
that $H\le 2N$ cannot be asserted without the content factor.
Its exact normalized arc constant is approximately $6.466$, so it
does **not** challenge an endpoint-arc improvement to (5).

For contrast, the short four-point source in
[the closed-translation obstruction](integer_cotangent_closed_translation_descent_obstruction.md)
has $A=157$, $L=1$, $d=(25,290)$, $c=(986,85)$,
$h_d=5$, $h_c=17$, and $H=290\ll N=212{,}298{,}125$.
The two cases show why the size of $h_dh_c$ relative to $g_A^2$
and $S$ must be derived from the endpoint and simultaneous pair
conditions, rather than assumed from complementarity.

The [checker](check_integer_cotangent_complementary_divisor_kernel.py)
verifies (2)--(7), the all-edge norms, common-scaling invariance, and
the two literal normalized arcs.
