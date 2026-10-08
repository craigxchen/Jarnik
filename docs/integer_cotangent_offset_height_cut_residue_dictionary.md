# The exact cut and residue height of complementary offsets

The normalized lcm of the offsets from the smallest positive cotangent has
an exact primewise description.  Its leading term is the width of the
Gaussian allocation interval visible between the circle anchor and that
smallest cotangent.  Primitive cotangent-chart imaginary residues account
for every remaining prime, including inert primes; a one-bit correction
accounts for the ramified prime two.  This yields an actual lower height
bound when the common integer cotangent parameter `L` is small.  It does
not establish the unproved exponent-four radius target, since that target
does not bound `L`.

Use the integer cotangent notation of
[the all-edge lcm target](integer_cotangent_lcm_height_target.md).  There
are `m>=2` distinct positive finite cotangents
`A=X_1<X_2<...<X_m`, an additional circle anchor `q_0=1`, and `L>0`.
Every pair cotangent

\[
 Q_{ij}=\frac{X_iX_j+L^2}{X_i-X_j}\in\mathbb Z
\]

is integral.  Let `N` be the all-edge lcm of the reduced pair norms,
`W=log N`, `d_i=X_i-A` for `i>1`, and

\[
 H=\frac{\operatorname{lcm}_{i>1}d_i}
          {\gcd_{i>1}d_i}.                                      \tag{1}
\]

This equals the complementary-divisor kernel `S/(h_d h_c)` with
`S=A^2+L^2` in
[the offset kernel note](integer_cotangent_complementary_divisor_kernel.md).

For an edge `e`, use its integer cotangent `Q_e`, with `Q_{0i}=X_i`, and
put

\[
 g_e=\gcd(|Q_e|,L),\quad r_e=L/g_e,\quad
 \eta_e=\mathbf 1_{\;Q_e/g_e\text{ and }r_e\text{ both odd}},\quad
 n_e=\frac{(Q_e/g_e)^2+r_e^2}{2^{\eta_e}}.             \tag{2}
\]

The positive integer `r_e` is the primitive imaginary numerator **in
this cotangent chart**.  More explicitly, with
`h_e=(Q_e/g_e+i r_e)/(1+i)^(eta_e)`, the factor `h_e` is a
conjugate-primitive Gaussian integer of norm `n_e`, and its phase is

```text
q_j/q_i=i^(eta_e) h_e/bar(h_e).
```

Thus `r_e=Im h_e` when `eta_e=0`, but
`r_e=Re h_e+Im h_e` when `eta_e=1`.  Reduction at `1+i`
introduces the displayed quarter-turn unit.  The residue correction
below keeps that unit exact rather than treating `r_e` as the
imaginary coordinate of every canonical physical pair factor.

## 1. Primewise formula

In the primitive physical Gaussian circle tuple of squared radius `N`,
only split odd primes occur in `N`.  For `p=pi bar(pi)` put
`e_j=v_pi(z_j)` at every physical row `j`.  Define, for `i>1`,

\[
 T_p(i)=\frac{|e_i-e_0|+|e_1-e_0|-|e_i-e_1|}{2}
       =\frac{v_p n_{0i}+v_p n_{01}-v_p n_{1i}}{2}.        \tag{3}
\]

The two forms agree by the all-edge primitive normalization.  `T_p(i)`
is a nonnegative integer: the distance from `e_0` to the projection of
`e_i` onto the segment `[e_0,e_1]`.  At an inert odd prime, an odd
split prime absent from `N`, and at `p=2`, set `T_p(i)=0`.

At `p=2` define

\[
 \delta_i=(\eta_{0i}+\eta_{01}-\eta_{1i})/2\in\{0,1\};
                                                               \tag{4}
\]

at odd primes put `delta_i=0`.  Then, for **every ordinary rational
prime**,

\[
 \boxed{\quad
 v_p(d_i)=T_p(i)+v_p(g_{01})
                  +v_p(r_{1i})-v_p(r_{0i})+\delta_i,
 \qquad
 v_p H=\operatorname{range}_{i>1}
       \bigl(T_p(i)+v_p(r_{1i})-v_p(r_{0i})+\delta_i\bigr).
 \quad}                                                       \tag{5}
\]

Here `range` means maximum minus minimum.  In particular, an inert
prime dividing `H` comes entirely from the exact residue difference;
no inert prime is hidden in a Gaussian allocation width.

To prove (5), the original cotangent identity is

\[
 (X_i^2+L^2)(A^2+L^2)
       =d_i^2(Q_{1i}^2+L^2).                              \tag{6}
\]

Substitute `Q_e^2+L^2=g_e^2 2^{eta_e} n_e`.  At odd primes,
half of the reduced-norm valuation difference is exactly `T_p(i)`;
at two, all reduced norms are odd and the remaining parity exponent
is (4).  Because the left and right sides of (6) are ordinary squares
times the displayed edge factors, the numerator in (4) is even, hence
`delta_i` is zero or one.  Finally
`v_p g_{0i}-v_p g_{1i}=v_p r_{1i}-v_p r_{0i}` because both `r` values
use the same `L`.  The anchor term `v_p g_{01}` cancels from the range.

## 2. Actual height bounds and their scope

Set

\[
 K_{01}=\sum_{\substack{p\mid N\\p\text{ split}}}
       \operatorname{range}_{i>1}T_p(i)\log p.           \tag{7}
\]

Every `r_e` divides `L`, so its valuation lies between zero and
`v_p L`.  The range of the residue correction in (5) differs from
the range of `T_p` by at most `2v_p L`; at two the parity correction
adds at most one.  Summing over all primes gives the unconditional
error bound

\[
 \boxed{|\log H-K_{01}|\le 2\log L+\log2.}            \tag{8}
\]

For a split prime, `range_i T_p(i)` is the length of the intersection
of the allocation intervals `[e_0,e_1]` and
`[min_{i>1}e_i,max_{i>1}e_i]`.  Consequently

\[
 0\le \log n_{01}-K_{01}\le\gamma_I,                \tag{9}
\]

where `I={2,...,m}` and `gamma_I` is the exact total threshold-layer
weight constant on all rows of `I`, as defined in
[projective orbit radius rigidity, Section 6](projective_orbit_radius_rigidity.md).
The lost part of `[e_0,e_1]` consists only of layers on which every
row in `I` is on one side of the threshold; (9) holds prime by prime.

Suppose the primitive physical tuple of radius `R=sqrt(N)>1` lies in an
arc of length at most `C sqrt(R)`.  Put

\[
 J=\lceil C/2\rceil,\qquad
 \ell=\left\lceil\frac{m-1}{4J}\right\rceil.            \tag{10}
\]

Partitioning `I` by short subarc and its four Gaussian units gives a
cohort of at least `ell` rows.  The obtuse sign-moment argument in the
cited note bounds `gamma_I<=W/ell`.  The literal endpoint chord from
`q_0=1` to `q_1=(A+iL)/(A-iL)` has length
`2sqrt(N) r_{01}/sqrt(2^{eta_{01}}n_{01})`.  Since `r_{01}>=1`, its
being at most `C N^(1/4)` implies

\[
 \log n_{01}\ge W/2+\log\frac{2}{C^2}.           \tag{11}
\]

Combining (8)--(11) proves the **actual endpoint lower bound**

\[
 \boxed{\log H\ge
   \left(\frac12-\frac1\ell\right)W
          -2\log L+\log\frac1{C^2}.}             \tag{12}
\]

For fixed `C` and a growing number of rows, this is
`log H>=W/2-O_C(W/m)-2log L-O_C(1)`.  Thus a small common `L`
in the selected chart forces `H` to have *large* height, of order
`sqrt(N)` up to the displayed errors.  Cut geometry alone cannot
make it `N^{O(1/m)}` while keeping `log L=o(W)`.  If `L` is
comparable with a power of `N`, (12) has no
positive leading consequence.  The integer target explicitly leaves
`L` unbounded, so this does not imply a general growth improvement.

There is a matching cut-only scope check.  In the formal binary
full-fair profile on the `m+1` rows, give every nonconstant unoriented
cut the same log weight `w`.  For any distinct rows `0,1`, the cut
contribution to (7) is `w` exactly when the cut separates `0,1`
and both signs remain among `I`.  There are `2^(m-1)-2` such cuts,
whereas `W=(2^m-1)w`; therefore

\[
 K_{01}=(2^{m-1}-2)w=(1/2+o(1))W.             \tag{13}
\]

This is a formal allocation calculation, not an asserted actual
endpoint family.  It explains why the exact allocation and chord
budgets do not themselves yield `log H=O(W/m)`.

For a literal check that both residue-only and ramified terms matter,
take `L=3` and `X=(3,4,9)`.  The three pair cotangents are
`Q_12=-21`, `Q_13=-6`, `Q_23=-9`, all integral.  Here
`d=(1,6)`, `H=6`, and `N=25`.  The inert prime `3` and ramified
prime `2` divide `H` while neither divides `N`: their allocation
terms `T_p` are zero.  At `3`, the difference of primitive chart
residues gives the exponent one; at `2`, `delta_3=1` gives it.

The [checker](check_integer_cotangent_offset_height_cut_residue_dictionary.py)
verifies (1)--(12) on 682 searched positive integral cliques and 16
literal/scaled fixtures, including an inert prime in `H` and a
ramified parity contribution.  It also checks 498 formal cuts for
(13).
