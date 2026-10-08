# Averaging balanced Gaussian characters with arbitrary block weights

The thirty-five balanced sign characters give a direct arbitrary-weight
inequality for eight-point profiles having only pair and balanced cuts.
This is a simultaneous arithmetic restriction, beyond the fixed rational
linear-template moment exclusion. For actual endpoint tuples, the existing
affine allocation rigidity theorem supplies all required nonvanishing at
`C<=2` with a common unit, or `C<=sqrt(2)` with arbitrary units. Thus the
result applies to arbitrary nested even-cardinality threshold layers in
these ranges. Even total row-unit parity is still required for the
arbitrary-unit version; it cannot be silently removed by taking fourth
roots.

## Exact hypotheses and conclusion

Suppose eight actual Gaussian points have the literal factorization

```
z_i = d epsilon_i product_j H_j^((1+s_ij)/2)
                           conjugate(H_j)^((1-s_ij)/2),
R=|d| product_j |H_j|,   Delta<=C/sqrt(R).
```

Assume each block `H_j` is coprime to its conjugate, as for canonical
split-prime threshold blocks. In particular all block norms are odd;
ramified and inert common factors belong in `d`. Every column of `S`
has minority size two or four. Write

```
P=sum_(pair cuts) log Norm(H_j),
B=sum_(balanced cuts) log Norm(H_j),
D=log Norm(d),  W=log R^2=D+P+B.
```

The source units may vary, but assume `product_i epsilon_i` is `+1` or
`-1`. For each balanced sign vector `lambda in {+1,-1}^8`, let

```
v_j(lambda)=(lambda^t S_j)/4.
```

These coefficients are integers. Assume that the conjugate-primitive
reduction `Beta_lambda` of the signed block product with exponents `v_j`
is a nonunit, for all thirty-five characters modulo global negation.
Then

```
                  B >= 5P + 35D - 140 log(sqrt(2) C).       (1)
```

In particular `2P<=B+140 log^+(sqrt(2)C)`.
If all source units agree, replace `sqrt(2)C` by `C` in (1).
The precise displayed inequality is valid even when its logarithmic
constant is negative.

A sufficient clean nonvanishing hypothesis is that all blocks are
conjugate-coprime nonunits with pairwise coprime rational norms, and no
balanced `lambda` annihilates `S`. More generally the theorem permits
shared primes and nested layers whenever the stated primitive products
are nonunit. It does not infer nonvanishing from the formal sign matrix
when its block primes overlap.

## Arithmetic of the fourth-root character

A balanced `lambda` has four plus and four minus entries. Since every
column also has even plus-count, `lambda^t S_j` is divisible by four.
The exact multiplicative identity is

```
product_i z_i^lambda_i
  = u_lambda (Beta_lambda/conjugate(Beta_lambda))^2,
u_lambda=product_i epsilon_i^lambda_i.                    (2)
```

The common factor cancels. Rational content discarded when forming the
conjugate-primitive `Beta_lambda` cancels from its conjugate ratio.
Writing `V_lambda=log Norm(Beta_lambda)`, one has

```
V_lambda <= sum_j |v_j(lambda)| log Norm(H_j),              (3)
```

with equality under the clean disjoint-support hypotheses.

Choose source argument lifts in their containing arc. Their balanced
signed sum has absolute value at most `4 Delta`. The unit parity of
`u_lambda` equals that of `product_i epsilon_i`, because all entries of
`lambda` are odd. Consequently `u_lambda=+/-1`, and (2) gives

```
dist(arg(Beta_lambda), (pi/4) Z) <= Delta.                 (4)
```

A conjugate-primitive nonunit Gaussian integer lies on neither a
coordinate axis nor a diagonal. After a Gaussian unit rotation, its
nonzero perpendicular coordinate is therefore an integer for an axis,
or a nonzero integer divided by `sqrt(2)` for a diagonal. For a closest
target direction in (4), this gives

```
1/sqrt(2) <= |Beta_lambda| sin(dist)
          <= |Beta_lambda| Delta.
```

Thus (3) implies the individual finite inequality

```
sum_j |v_j(lambda)| log Norm(H_j)
  >= W/2 - 2 log(sqrt(2) C).                              (5)
```

If the row units agree, `u_lambda=1`; the target grid is `(pi/2)Z`, and
the perpendicular coordinate has absolute value at least one. This gives
(5) with `C` in place of `sqrt(2)C`.

When the total unit parity is odd, all `u_lambda` are `+/-i`. Their target
directions instead lie in the odd multiples of `pi/8`, with irrational
quadratic slopes. The rational-axis lower bound used in (5) is then
unavailable. No inequality (1) is claimed for that remaining branch.

## The exact averaging certificate

Choose the thirty-five balanced `lambda` with first entry `+1`. Global
negation does not affect absolute coefficient values. For a fixed pair
cut, `|lambda^t S_j|/4` is one exactly when the two minority positions have
the same `lambda` sign, and zero otherwise. There are fifteen such
characters. For a fixed balanced cut the coefficient distribution is

```
absolute value       0    1    2
number              18   16    1.
```

Its total is eighteen. Therefore the exact average over all characters is

```
sum_lambda sum_j |v_j(lambda)| log Norm(H_j) = 15P+18B.
```

Summing (5) gives

```
15P+18B >= (35/2)(D+P+B)-70 log(sqrt(2)C),
```

which rearranges to (1). This is an explicit positive dual certificate;
there is no optimization tolerance, equal-weight assumption, or limiting
argument.

## Actual even-layer profiles need no independence hypothesis

Normalize an actual eight-point tuple by its complete common Gaussian
factor `d`, and represent each varying split-prime allocation by its
nested threshold layers. Suppose every nonconstant layer has even
cardinality, hence minority size two or four. Grouping repeated layers
is optional. The factorization above is literal with `H_j` equal to a
power of the corresponding chosen Gaussian prime; blocks may overlap
in rational prime support.

For distinct actual endpoint points, [affine allocation rigidity](linear_allocation_affine_rigidity.md)
proves that their varying-prime allocation vectors are affinely independent
when `C<=2` and all units agree, or `C<=sqrt(2)` with arbitrary units.
If any `Beta_lambda` were a unit, (2) would be a multiplicative relation
modulo units; prime valuations would give

```
sum_i lambda_i a_i(p)=0 for every split prime p,
sum_i lambda_i=0,
```

contradicting that theorem. The unit-radius exception contains fewer than
eight distinct points in these arc ranges. This proves nonvanishing for
all thirty-five characters, without using the rank of the formal layer
matrix or disjoint block supports.

Consequently every actual common-unit eight-point endpoint tuple with
`C<=2` and only even nonconstant threshold layers obeys

```
B >= 5P + 35D - 140 log C.                               (6)
```

With arbitrary units, `C<=sqrt(2)`, and even total row-unit parity, it obeys
(1). Both conclusions allow arbitrary prime powers, nested layers,
unequal weights, and complete common content. They prove the intrinsic
inequality for this actual subprofile. They do not claim that a general
eight-point endpoint tuple has only even layers.

## Positive order-twenty-four profiles

Let `T` be an eight-row restriction of an order-twenty-four Hadamard
matrix, and suppose its columns consist of one constant column, eight
pair cuts, and fifteen balanced cuts. Delete its constant column and one
balanced column `b` to obtain `S`. Then

```
SS^t=24 I-1 1^t-b b^t.
```

The vectors `1,b` are orthogonal and have squared norm eight. Thus the
smallest eigenvalue of `SS^t` is sixteen. Every nonzero `lambda` has
`lambda^t S != 0`, proving the matrix part of clean nonvanishing.
Independent conjugate-primitive nonunit blocks therefore satisfy all the
hypotheses above whenever their row units have even total parity.

This covers every one of the 11385 positive Paley order-twenty-four codes
in [the exact fixture](paley24_even_column_moment_fixture.json), and the
critical tensor profile, with arbitrary positive block weights. Their
full clean endpoint realization consequently satisfies (1), a stronger
bound than the requested intrinsic inequality `2P<=B+O_C(1)`.

The general uniform lattice-circle endpoint problem remains open. In
particular the argument does not extract these codes, remove the odd
unit-parity branch, or establish nonvanishing outside the stated small-constant
affine-rigidity regimes for arbitrary overlapping prime layers.

[The checker](check_even_column_balanced_character_height.py) verifies the
exact averaging table for every even nonconstant eight-row sign column,
the Hadamard rank identity for all 11385 stored codes, and (2) on literal
Gaussian products with independent units.
It also checks 840 signed-character identities on repeated-prime nested
layers, including primitive reduction with nontrivial rational content.
These latter fixtures verify the exact cancellation and height inequality;
they do not assert that the generated points lie on endpoint arcs.
