# A nonlinear three-bit permutation gives a balanced aligned foliation

This note concerns only the repeated-Walsh model. It gives an explicit
nonlinear assignment with polynomial-size aligned flats and an
arbitrary-label restriction balance. An explicit polynomial growth
exponent under arbitrary physical prime weights remains open for this
family.

Write three-bit vectors as integers `0,...,7`, with the usual binary dot
product, and define the permutation

```text
u:       0 1 2 3 4 5 6 7
f(u):    0 1 3 6 7 4 5 2.                                    (1)
```

Its coordinate polynomials, where `u=u0+2u1+4u2`, are

```text
f0=u0+u1+u2+u1u2,
f1=u1+u2+u0u2,
f2=u2+u0u1.                                                (2)
```

The quadratic part is the gradient of `u0u1u2`. Direct expansion gives

```text
(f(u)+f(u+w)) dot w
 = w0+w1+w2+w0w1+w0w2+w1w2+w0w1w2
 = 1_(w!=0),                                                 (3)
```

independently of `u`. In particular `f` is genuinely nonlinear and has
no aligned two-row pair in any nonzero direction: the two restriction
bits differ by one.

It is quantitatively far from affine. For each nonzero `w`, the eight
values `f(u)+f(u+w)` consist of four values, each occurring twice.
Consequently five points of agreement with an affine map are impossible:
among their ten unordered pairs, two share a difference `w` (there are
only seven nonzero differences); those two disjoint pairs form a
parallelogram, but an affine map has the same derivative on both pairs
whereas `f` gives each derivative value on only one unordered pair.
Thus any affine three-bit map agrees with `f` at most four times.
Fixing all but one input block of (4) shows that any global affine
map agrees with `a0` on at most `M/2` rows. The zero-label repair can
increase this count by at most one, so `a` differs from every affine
assignment on at least `M/2-1` rows. In particular the polynomial-size
foliation below is outside the sublinear-exception near-affine theorem.

Now take `2m` independent three-bit blocks, `m>=1`, so
`t=6m`, `M=2^(6m)`. Before repairing the zero label, assign

```text
a0(u1,...,u_(2m))=(f(u1),...,f(u_(2m))).                      (4)
```

Choose any perfect matching of the blocks. On each matched pair
`(i,j)`, choose **independent** nonzero directions `wi,wj` in `F_2^3`
and let `v_(i,j)` have `wi,wj` in those blocks and zero elsewhere.
The `m` vectors generate a row subspace `V` of size

```text
h=2^m=M^(1/6).                                               (5)
```

For every row `x`, (3) on the two blocks gives

```text
[a0(x+v_(i,j))+a0(x)] dot v_(i,j)=1+1=0.                     (6)
```

Adding any other generator does not change this scalar, because its
block support is disjoint. Thus `a0(x)|V` is constant on **every**
row coset `x+V`, including cosets with nonzero restriction. Since `a0`
is a permutation of all `M` labels and restriction to `V` is onto,
each `ell in V*` occurs on exactly `M/h` rows, or

```text
n=M/h^2=2^(4m)                                               (7)
```

whole row cosets. Every nonzero `ell` gives `n` disjoint aligned cosets;
together they contain `(1-1/h)M` rows. This holds for every choice of
matching and independent block directions.

The repeated-Walsh assignment cannot use label zero. Replace only
`a0(0)=0` by any fixed nonzero label `c`, obtaining `a`. Every nonzero
label is used at most twice. The changed row lies in the former
zero-restriction coset, so all the `n(h-1)` nonzero-restriction aligned
cosets and their exact counts survive. Independent physical prime
columns can be chosen for the duplicate label when `b>=5`.

The freedom to choose the two directions independently matters. Given
any nonzero label `beta=(beta_1,...,beta_(2m))`, choose a matched pair
with a nonzero component, set its direction so that `beta|V` is nonzero,
and choose the remaining directions arbitrarily. For example, if
`beta_i!=0`, choose `wi` with `beta_i dot wi=1`; if the other component
is also nonzero, choose `wj` with `beta_j dot wj=0`. Hence every
nonzero label belongs to a nonzero restriction fiber for **some**
perfectly aligned foliation. Coverage across varying `V` does not by
itself balance the physical weights in one certificate mixture.

Under the actual near-equal physical prime-log hypothesis
`u<=log p_j<u+1/r`, `u>=log 5`, `r=b(M-1)`, the existing aligned-flat
phase certificate applies to any of the surviving flats. For fixed
arc constant `C` and sufficiently large `m`, its exact height
`bM/2-h` forces `b>h`. Distinct physical primes then give

```text
W=log R^2 >= c_C h M log M = c_C M^(7/6) log M,             (8)
M=O_C((log R/loglog R)^(6/7))                                (9)
```

above a fixed logarithm cutoff. The same explicit exponent does not
follow for arbitrary prime weights: the aligned nonzero fibers together omit the
zero-restriction label block for any fixed `V`, leaving the baseline
mass short of `W0/2`. A complementary pair or outer character family
would have to fill that block while preserving the flip reduction.

The bounded checker
[`check_walsh_nonlinear_triple_gadget_foliation.py`](check_walsh_nonlinear_triple_gadget_foliation.py)
verifies (1)--(7), independent directions, the zero-row repair, and
copy capacity for `m=1,2`. It does not assert that its fixtures are
endpoint arcs.

For arbitrary prime weights this family also satisfies the general
[Ramsey growth bound](walsh_arbitrary_weight_ramsey_growth.md).
Its stronger explicit exponent `6/7` above still requires the
near-equal-prime hypothesis.
