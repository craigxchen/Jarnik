# Complement symmetry and the three fixed-pair cuts

The common Gaussian matching-content lower bound for a standalone
six-point circle is exactly invariant under complementing an old core
cut. This is a covariance statement for the actual degree-seven lifts,
not an invariance of their generic minimum orders. It also preserves
the total local correction error. The cuts `12,13,23` and their three
complements each force contribution at least `e/2-E_p` in units of
`log p`, for arbitrary real nonnegative valuation errors.

This note proves these local statements. Other cut cases and their
global radius consequences require separate arguments. The additional
continuous local bounds are proved in
[joint_involution_mixed_cut_content.md](joint_involution_mixed_cut_content.md)
and [joint_involution_two_pair_cuts.md](joint_involution_two_pair_cuts.md).

## 1. Standalone matching content and its gap formula

For six signed relative phase valuations `d_1,...,d_6`, let
`d_(1)<=...<=d_(6)` be their increasing rearrangement. On the least
standalone integral circle the local allocation width is

```text
n=d_(6)-d_(1).
```

The marked old anchor is not included in this range. By the two
Gaussian-orientation matching bound from
[joint_circle_matching_content.md](joint_circle_matching_content.md),
the prime contributes at least `Gamma(d) log p` to the logarithmic
modulus of the common Gaussian factor, where

```text
Gamma(d)=(3n+d_(1)+d_(2)+d_(3)-d_(4)-d_(5)-d_(6))/2.
```

Writing `g_j=d_(j+1)-d_(j)` gives the exact identity

```text
Gamma=g_1+g_2/2+g_4/2+g_5.                              (1)
```

Thus `Gamma` is nonnegative, and

```text
Gamma>=(d_(3)-d_(1))/2,
Gamma>=(d_(6)-d_(4))/2.                                 (2)
```

It is invariant under a common translation of all six phases and under
their simultaneous negation. These assertions hold for arbitrary real
phase valuations; no discreteness or bounded range is required.

## 2. An exact complement operation on the inputs

Fix an odd split prime `p=pi bar(pi)` and a core exponent `e>0`.
Write `s_i=1` on a cut `S` and zero outside. Allow the following local
data for the six Gaussian inputs:

```text
v_pi(P_i)=e s_i+alpha_i,
v_barpi(P_i)=beta_i,
v_p(Delta_ij)=e s_i s_j+delta_ij,                        (3)
```

where all `alpha_i,beta_i,delta_ij` are nonnegative. In particular we
do not require conjugate primitivity in this local formulation. Set

```text
E_p=sum_i(alpha_i+beta_i)+sum_(all15 edges)delta_ij.
```

For the actual centrally truncated inputs,
`alpha_i+beta_i=v_p(N(K_i))` and `delta_ij=v_p(b_ij)`.
Consequently `sum_p E_p log p<=12s+15beta` under the previous global
height bounds. Adding `v_p(2)` is harmless if using the error convention
in the other notes; all core primes under consideration here are odd.

Let `H=pi^e`, and transform each row by

```text
R_i=H bar(P_i)/p^(e s_i).                               (4)
```

These rows are locally integral: inside `S`, the conjugate `bar(P_i)`
is divisible by `bar(H)`, so (4) equals `bar(P_i)/bar(H)`; outside `S`
it is simply `H bar(P_i)`. Their new local data are

```text
v_pi(R_i)=e(1-s_i)+beta_i,
v_barpi(R_i)=alpha_i,
Delta_ij(R)=-p^(e(1-s_i-s_j)) Delta_ij(P),
v_p(Delta_ij(R))=e(1-s_i)(1-s_j)+delta_ij.               (5)
```

Thus the cut becomes its complement, the orientation errors
`alpha_i,beta_i` swap, and every bracket error remains unchanged.
In particular `E_p` is preserved exactly. The transformed rows can
have common conjugate factors; this is why imposing conjugate
primitivity again would be an unnecessary and potentially incorrect
extra step.

## 3. Covariance gives complement invariance of the actual output

Let `Q_j`, for `j=4,5,6`, be the raw degree-seven covariants in
[joint_circle_conductor_lift.md](joint_circle_conductor_lift.md).
Each has degree two in row `j` and degree one in every other row.
A common real binary transformation `A` gives

```text
Q_j(A P)=det(A)^3 A Q_j(P).
```

This also follows directly because every summand of `Q_j` contains
three real brackets and one Gaussian row. In (4), the common real
linear transformation is `z -> H bar(z)`, with determinant `-p^e`,
followed by row scalars `p^(-e s_i)`. Therefore

```text
Q_j(R)=-p^(e(3-|S|-s_j)) H bar(Q_j(P)).                 (6)
```

The real scalar in (6) has equal valuations at the two Gaussian
orientations even when its exponent is negative. Thus each of the six
output signed phases, including the unchanged first three rows, obeys

```text
d_i(R)=e-d_i(P).                                        (7)
```

Equations (1) and (7) prove exact complement invariance of `Gamma`.
Accordingly, any local inequality proved under (3) for a cut `S`
transfers to its complement with the same correction bound. The proof
concerns the standalone six-point allocation range. The range including
an additional fixed zero-phase anchor is not invariant under (7).

## 4. The three pure fixed-pair cases

We use both forms of the exact lifts:

```text
Q_4=Hform P_4-2 Delta_34 Delta_56 Delta_41 P_2
   =V P_4-2 Delta_34 Delta_56 Delta_42 P_1,
Q_5=F P_5-2 Delta_35 Delta_46 Delta_51 P_2
   =V P_5-2 Delta_35 Delta_46 Delta_52 P_1,
Q_6=Hform P_6+2 Delta_36 Delta_45 Delta_61 P_2
   =F P_6+2 Delta_36 Delta_45 Delta_62 P_1.               (8)
```

Here `Hform=B+C` is the Segre linear form, distinguished from the
Gaussian number `H=pi^e` in (4). The other forms are
`F=C-B` and `V=2A+B+C`.

For the cut `12`, the core orders of `(A,B,C)` are `(e,e,0)`.
The first and third forms in the original `P_2` expressions are
`Hform`, and the second is `F`. Their valuations are bounded by
`sum delta_ij`: each is dominated by the core-order-zero matching `C`.

For `23`, the core orders are `(0,0,e)`. The same three expressions
use forms dominated by the core-order-zero matching `B`.

For `13`, use the alternate `P_1` expressions in (8). Their forms are
`V,V,F`. The exact matching identities

```text
A+C=Delta_13 Delta_24 Delta_56,
A+B=Delta_12 Delta_35 Delta_46,
F=(A+C)-(A+B),       V=(A+C)+(A+B)                       (9)
```

have core orders `e,0` in the first two lines. Thus `F,V` are
dominated by the core-order-zero matching `A+B`.

Suppose `e>E_p`. In each case every moving row is outside the cut.
The first summand of each chosen expression (8) has first-orientation
valuation at most

```text
sum delta_ij+alpha_j<=E_p.
```

The other summand contains an inside reference row (`P_2` or `P_1`),
so its first-orientation valuation is at least `e`. These are strictly
unequal orders. It follows that `v_pi(Q_j)<=E_p`, and hence
`d_j<=E_p`, since the other orientation is nonnegative.

The one retained row outside the cut also has phase at most `E_p`.
The two retained rows inside have phase at least `e-E_p`. There are
therefore at least four phases at most `E_p`, and some phase at least
`e-E_p`. Applying (2) gives

```text
Gamma>=e/2-E_p.                                         (10)
```

If `e<=E_p`, (10) follows immediately from `Gamma>=0`. Thus it is
valid without any restriction on the correction size. By (7), the same
bound holds for the complementary cuts `3456,2456,1456`, respectively.
Together these six blocks contribute at least `3(1-eta)w` minus their
summed correction heights. Their role in a stronger all-cut lower bound
must still be combined with separately proved cases.

## Verification scope

The degree-seven polynomial identities are certified by
[check_joint_circle_lift.py](check_joint_circle_lift.py).
The complement transformation, the error preservation, and the arguments
in Sections 1--4 are symbolic identities and strict valuation comparisons;
they do not infer actual valuations from a finite experimental grid.
The symmetry agent independently audited the full covariance exponent,
the error swap, the invariance of `Gamma`, and all three pure fixed-pair
cases. The mixed-cut and two-pair error proofs must use the same
two-orientation hypothesis (3) when transferring to complementary cuts.
