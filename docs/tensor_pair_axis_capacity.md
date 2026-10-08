# Paired Gaussian factors have bounded simultaneous axis capacity

A literal paired factorization gives a uniform bound on its base factor,
with arbitrary block weights and arbitrary Gaussian row units. Applied to
the critical eight-row `H12 tensor H2` restriction, this proves the desired
intrinsic weight inequality under the disjoint-block hypotheses below.
It does not extract this tensor profile from a general endpoint tuple.

The projection argument extends the axis-capacity calculation in
[pair_reflection_domain_capacity.md](pair_reflection_domain_capacity.md).
The common-unit small-constant tensor exclusion also already follows from
the repeated pair squareclass obstruction in
[angular_pair_squareclass_independent_audit.md](angular_pair_squareclass_independent_audit.md).
For arbitrary units and `C<1`, the tensor exclusion also follows from
augmented four-wise label separation in
[squareclass_baseline_and_power_groups.md](squareclass_baseline_and_power_groups.md):
the four paired edges share their allocation squareclass and have only two
possible unit-parity bits, so two edges give a four-point augmented relation.
The quantitative application here bounds the base factor for every fixed
endpoint constant and explains the simultaneous projection mechanism. It is
not a new qualitative small-constant exclusion or a general growth improvement.

## Literal paired-factor theorem

Suppose, for `1<=i<=m`,

```
z_(i,+) = d epsilon_(i,+) B_i A_i,
z_(i,-) = d epsilon_(i,-) B_i conjugate(A_i),
|B_i|=b,  |A_i|=a,  R=|d|ba,
```

where all displayed factors are nonzero Gaussian integers and the
`epsilon` are arbitrary Gaussian units. Assume that no two of the `A_i`
are associates or conjugate associates. Suppose all source points lie in
an arc of angular width `Delta`, with `Delta sqrt(R)<=C`.
If `m>=3`, then

```
                         |d| b <= C^2/8.                 (1)
```

In particular such a family cannot occur for `C<sqrt(8)`.
The hypothesis about the `A_i` is essential; no distinctness assertion is
being inferred merely from distinct source points.

For the proof first take `Delta<pi`. Choose source argument lifts in the
containing arc. The ratio of each pair gives

```
z_(i,+)/z_(i,-) = u_i A_i/conjugate(A_i),  u_i in mu_4,
dist(arg(A_i), (pi/4) Z) <= Delta/2.                       (2)
```

After multiplying each `A_i` by a Gaussian unit, assign it to either the
axis direction zero or the diagonal direction `pi/4`, at angular distance
at most `Delta/2`. Any witnessing target in (2) is sufficient; uniqueness
of the nearest target is unnecessary. Projection on the assigned axis
lies in the interval

```
[a cos(Delta/2), a],  of length s=a(1-cos(Delta/2)).        (3)
```

In the coordinate class projection levels are integers. In the diagonal
class they are `(x+y)/sqrt(2)` for Gaussian integers `x+iy` of fixed norm
`a^2`. Since `x+y == x^2+y^2 == a^2 (mod 2)`, their spacing is at least
`sqrt(2)`. At a fixed projection level the two possible points on the
circle are reflected across that axis. They are conjugate associates:
`z` and `bar(z)` for the coordinate axis, and `z` and `i bar(z)` for the
diagonal. Consequently at most one selected index occupies each level.
This proves the more precise capacity bound

```
m <= 2 + floor(s) + floor(s/sqrt(2)).                      (4)
```

If `m>=3`, at least one class has two levels, so `s>=1`. Using
`1-cos t<=t^2/2` gives

```
1 <= s <= a Delta^2/8 <= C^2/(8 |d| b),
```

which proves (1). If `Delta>=pi`, the endpoint assumption instead gives
`|d|b<=R<=C^2/pi^2<C^2/8`, completing the proof without an angle restriction.
The argument is unchanged when factors share primes, provided the literal
paired factorization and conjugation-injectivity hypotheses hold.

## Application to the critical tensor profile

Select four rows of a normalized order-twelve Hadamard matrix, and retain
both rows of the order-two tensor factor. Delete the constant tensor
column and one balanced tensor column. Suppose the resulting eight-row
code has eight minority-size-two columns and fourteen balanced columns,
as in the selected rows `(0,1,2,4)` of the Paley fixture.

Use a nonunit Gaussian block for every retained column, each coprime to
its conjugate, with pairwise coprime rational norms across columns. Allow
an arbitrary common Gaussian factor `d` and independent Gaussian row
units. Group the columns according to whether their `H2` factor is
constant or alternating. Their products give exactly the `B_i` and `A_i`
of the theorem. All pair cuts belong to the base group; every alternating
column is balanced.

Two distinct order-twelve Hadamard rows differ in six positions. The
retained alternating sign rows have length twelve and mutual distance
six if a base column was deleted; if an alternating column was deleted,
they have length eleven and mutual distance five or six. Thus no two
are equal or complementary. Disjoint conjugate-primitive nonunit block
supports imply that their Gaussian products are neither associates nor
conjugate associates: compare the valuations at a prime in a differing
coordinate, or in an agreeing coordinate, respectively.

Hence (1) applies with `m=4`, for arbitrary positive block weights
`w_j=log Norm(H_j)`. Write `W_pair` for the sum of the eight pair-cut
weights and `W_bal` for the remaining fourteen weights. The base contains
all pair cuts, so

```
W_pair <= log(b^2) <= 2 log(C^2/8),
2 W_pair - W_bal <= 4 log^+(C^2/8).                       (5)
```

For `C<sqrt(8)` there is no such factorization at all. Neither equal block
sizes, rational linear templates, nor moment expansions enter (5).

For nested prime layers, disjointness cannot simply be dropped from the
injection proof: signed exponents may cancel. Equation (5) still follows
whenever the literal grouping produces conjugation-injective `A_i`.
No claim is made here that arbitrary nested layers guarantee that property.

## An exact barrier for scalar integer-row averaging

For the same Paley rows, delete tensor columns zero and one. There are
four column classes: eight base pair cuts, three base balanced cuts,
eight alternating columns coming from singleton four-row cuts, and three
alternating columns coming from balanced four-row cuts. Give each column
in these classes respective formal weights

```
                         (5,4,6,4).
```

Then `W=112`, `W_pair=40`, and `W_bal=72`, so `2W_pair>W_bal`.
Nevertheless every nonzero zero-sum integer row vector `c` satisfies

```
v=c^t S/2 in Z^22,    V(c)=sum_j |v_j|w_j >=56=W/2.       (6)
```

Here is an exhaustive proof with a bounded finite part. The sixteen
columns in the first and third classes form a doubled order-eight
Hadamard matrix, up to column orientations, so their row Gram is `16I`.
Hadamard inversion gives

```
sum_(these columns) |c^t S_j/2| >= 8 max_i |c_i|.
```

If the maximum is at least two, their weights are at least five and
`V(c)>=80>56`. If the maximum is one, there are exactly 1106 nonzero
zero-sum vectors in `{-1,0,1}^8`. The exact checker gives minimum heights
`56,56,84,72` at supports `2,4,6,8`, respectively, proving (6).

For actual independent blocks, an integer-row character has the exact
phase relation

```
product_i z_i^c_i = u (Beta/conjugate(Beta)),
Beta=product_j H_j^(v_j positive part)
               conjugate(H_j)^(v_j negative part),
log Norm(Beta)=V(c).
```

The common factor cancels because `sum c_i=0`. Its argument approaches
the rational-slope grid `(pi/4)Z`; conjugate primitivity and nonunit
support prohibit exact equality. The nonzero integral axis/diagonal
coordinate gives the usual `V(c)>=W/2-O_(C,c)(1)` bound. Thus (6) is a
leading-order dual barrier to deducing (5) by nonnegative averaging of
these scalar integer-row half-height bounds. It is not an actual endpoint
example, a barrier to all rational characters with branch information,
or a contradiction to the simultaneous projection theorem.

Run [the checker](check_tensor_pair_axis_capacity.py). It verifies
the tensor profile, all deletion locations, injection on actual Gaussian
products and the complete scalar barrier with integer arithmetic. Its
additional geometric projection-capacity fixtures use floating-point
angles and supplement the proof above.
