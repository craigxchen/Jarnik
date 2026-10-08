# Four fixed norm directions and the finite-cover barrier

For the disjoint elliptic integer family, the twenty-five oriented
contact congruences reduce exactly to four norm divisibilities at
**fixed rational directions**. This gives an explicit small-residual
system for the first frame. A separate Siegel argument proves that
an infinite family of integral frames cannot be parametrized by
any fixed finite algebraic cover of the elliptic source. The bound
of two short frames over each rational configuration does not
provide such a cover.

The integer family and its normalization are fixed as in
[the contact-profile lemma](disjoint_mixed_elliptic_integer_contact_profile.md).
Put `w_N=q(P)N^2`, so every contact integer has logarithm
`w_N+O(N)`. All primes below are outside a fixed exceptional set
containing two and every prime needed in the local coordinate models.

## 1. Four fixed directions replace twenty-five variable slopes

Group the contact integers by their ordinary collision direction:

```text
D_infinity = product_(seven contacts with anchor infinity) d_C,
D_0        = product_(seven contacts with anchor zero) d_C,
D_1        = product_(seven contacts with anchor one) d_C,
D_rho      = product_(four unanchored contacts) d_C,
rho=-3/2.
```

They are pairwise coprime and satisfy

```text
log D_infinity = log D_0 = log D_1 = 7w_N+O(N),
log D_rho = 4w_N+O(N),
log(D_infinity D_0 D_1 D_rho)=25w_N+O(N).              (1)
```

At a good prime power `p^k` belonging to a contact `C`, the primitive
coefficient rows in that cluster are unit multiples, modulo `p^k`,
of the corresponding fixed primitive vector

```text
(1,0),       (0,1),       (1,1),       (-3,2).          (2)
```

For anchored contacts this is immediate from the fixed anchor.
For the fourth fiber the coordinate difference from `rho` has
valuation `k`: the inverse slopes at that common root are nonzero
units outside the fixed exceptional set. The assertion therefore
retains the full contact depth, not only reduction modulo `p`.

Let an integral frame be `U=(z,w)`, with `z,w in Z[i]` and
`Im(bar(z)w)=1`. Across all choices of Gaussian orientations, its
shared contact congruences are equivalent to

```text
D_infinity | Norm(z),
D_0        | Norm(w),
D_1        | Norm(z+w),
D_rho      | Norm(-3z+2w).                            (3)
```

To see the converse as well as necessity, an `SL_2(Z)` matrix sends
each vector in (2) to a primitive integer row. At an odd split
prime its Gaussian value can be divisible by at most one of the
two conjugate primes. Divisibility of its norm by `p^k` therefore
selects one unique orientation to the full depth `k`. All rows in
the cluster share it by (2). At an inert prime such norm divisibility
is impossible for a primitive row. Thus (3) also detects the exact
inert-prime obstruction; when all four moduli are split-supported,
it is precisely the unoriented version of the shared congruence
system. It does not impose private factors or residual row sizes.

## 2. Exact residual equations at the proposed short-frame scale

Write

```text
A=Norm(z),        C=Norm(w),        B=Re(bar(z)w),
K=Norm(z+w),      L=Norm(-3z+2w).
```

The determinant-one identity and the four fixed directions give

```text
AC-B^2=1,        K=A+2B+C,
L=9A-12B+4C=15A+10C-6K.                              (4)
```

Hence every admissible frame determines positive integers
`R_infinity,R_0,R_1,R_rho` with

```text
A=D_infinity R_infinity,    C=D_0 R_0,
K=D_1 R_1,                 L=D_rho R_rho,

15D_infinity R_infinity + 10D_0 R_0
       -6D_1 R_1 -D_rho R_rho = 0,                    (5)

4D_infinity D_0 R_infinity R_0
       -(D_1 R_1-D_infinity R_infinity-D_0 R_0)^2=4.   (6)
```

Suppose `||U||_F <= exp(4w_N+o(w_N))`. Equations (1)--(4) yield

```text
log R_infinity, log R_0, log R_1 <= w_N+o(w_N),
log R_rho <= 4w_N+o(w_N).                             (7)
```

For example `K<=2||U||_F^2` and `L<=13||U||_F^2`.
The three smaller residuals are the norm-size budget left at the
three normalized anchors. Equation (6) retains determinant one;
equation (5) imposes the fourth, nonanchor modulus with its full
weight. This is a concrete simultaneous norm problem with coefficients
given by the elliptic contact sequence.

Conversely, positive integer solutions of (5)--(6), with
`K-A-C` even, give a positive definite integral binary Gram matrix
of determinant one. Such a binary form is integrally equivalent
to the identity and is the Gram matrix of an `SL_2(Z)` frame.
For completeness, reduce a positive integral determinant-one form
so that `|2B|<=A<=C`. If `A>=2`, then
`AC-B^2 >= 3A^2/4 >=3`, a contradiction. Therefore `A=1`,
and an integral shear makes the form the identity. Adjusting an
orthogonal sign gives determinant `+1`. The frame then satisfies
(3). Its exact height still requires `A+C` to meet the stated bound.

There is no proved contradiction in (5)--(7). The residuals are
allowed to grow exponentially in `N^2`; none can be declared a
fixed Gaussian factor or a fixed rational function after passing
to an infinite subsequence. The general CRT theorem still gives
an integral frame at unrestricted height for every split-supported
choice of these four moduli.

The grouped product in (1) exceeds `H^6=exp(24w_N+o(w_N))`, so the
[coupled Gram theorem](coupled_orientation_frame_packing.md)
allows at most two Gaussian-unit classes in the prescribed height
window at each fixed `N`. This counts exceptional candidates; it
does not eliminate the first one.

## 3. A fixed finite-cover parametrization would force finiteness

Here is a precise restriction on any attempted algebraic extraction.

**Lemma.** Let `E/Q` be an elliptic curve. Let `Q_n in E(Q)` be
distinct, and suppose `U_n in SL_2(Z)` are integral frames satisfying

```text
M_n | Norm(first Gaussian column of U_n),    M_n -> infinity.
```

Then every fixed algebraic curve in `E times SL_2` contains at most
finitely many pairs `(Q_n,U_n)`. In particular, if there are infinitely
many pairs, their Zariski closure has dimension at least two. They
cannot all arise from rational maps on a fixed finite algebraic
cover of `E`, or on a fixed finite collection of such covers.

**Proof.** A vertical curve contains at most one of the pairs,
since the base points are distinct. Consider an irreducible curve
dominating `E`, and take its smooth projective normalization `C`.
The induced nonconstant morphism `C -> E` is finite. Riemann--Hurwitz
gives `g(C)>=1`.

If one matrix entry is nonconstant on `C`, call that rational
function `f`. At every pair under consideration, its value is an
integer. Siegel's theorem for a nonconstant rational function on a
curve of positive genus says that only finitely many rational points
of `C` have `f` integral. Removing the finitely many singular points
does not change the conclusion. The precise rational-function
form used here is [Levin, Theorem 1.1](https://jtnb.centre-mersenne.org/item/10.5802/jtnb.818.pdf).

If every entry is constant, the frame is fixed. Its first column
has a fixed positive norm, which cannot be divisible by `M_n` for
unbounded `M_n`. This also gives only finitely many pairs. A fixed
finite union of curves is handled component by component. This
proves the lemma.

Applying it with `M_n=D_infinity(N)` proves the assertion for any
infinite first-frame family on the disjoint elliptic curve. The
lemma uses no assumption that prime orientations are rational
functions, and no ineffective bound is presented as quantitative.

This identifies the exact gap in a tempting finite-cover argument.
A pointwise bound of two on a selected set of integral points is
not a bound on the fibers of its **Zariski closure**. An arithmetic
graph can be Zariski dense in a surface. Here any hypothetical
infinite family would necessarily have that higher-dimensional
behavior, by the lemma. Proving containment in finitely many fixed
algebraic multisections would be a new exclusion theorem, not a
formal consequence of uniqueness or of the elliptic base.

The norm system (5)--(7) and this finite-cover obstruction do not
prove that inert mass, private factors, or Gaussian corrections
must grow. They isolate additional exact arithmetic and explain
why the short-frame counting theorem alone cannot be promoted to
an infinite-family exclusion. The global uniform lattice-arc bound
remains unresolved.

A subsequent [fixed rational-frame height theorem](fixed_rational_frame_contact_height.md)
gives a quantitative exclusion for rational matrix formulas on this
elliptic source, even with subpower common denominators and contact
errors. It does not exclude pointwise choices of varying frames in
the system above.
