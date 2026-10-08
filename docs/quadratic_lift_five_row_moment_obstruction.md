# A real moment obstruction to quadratic lifting a full profile

The certified full four-row profile in
[four_row_real_polynomial_profile_certificate.md](four_row_real_polynomial_profile_certificate.md)
cannot be extended to a full five-row profile by pulling its four rows back
through any real quadratic polynomial and adding a fifth row. This excludes
the entire quadratic lifting mechanism for that profile, including every
reanchoring, real affine change of its base variable, and global conjugation. It does not
exclude an arbitrary full five-row profile: its four-row restrictions need
not be quadratic pullbacks of this particular four-row profile, or admit a
common quadratic composition at all.

The obstruction is over the real numbers and uses three exact Newton
moments. It does not use rationality or an expected-dimension argument.

## 1. The general quadratic lifting test

Let `m>=3`. Suppose a full `m`-row degree-one profile has distinct nonreal
roots `z_T`, indexed by the nonempty subsets of `[m]`, with all roots and
their conjugates distinct. Thus each row product

```
P_i(u) = product_(T contains i) (u-z_T)
```

has all its nonconstant coefficients real and has a nonzero constant
imaginary part. The degree of each row is `2^(m-1)`.

Let `q(t)` be any real quadratic polynomial. Suppose that a full
`(m+1)`-row degree-one profile has its first `m` rows proportional, by real
nonzero constants, to `P_i(q(t))`. The new rows have degree `2^m`.
For each nonempty `T subset [m]`, the two preimages of `z_T` under `q`
must then be precisely the roots of the blocks labelled `T` and
`T union {m+1}`. They are distinct by the full-profile hypothesis. The
last row contains exactly one of them. It also contains the single root
`r` belonging to its private block `{m+1}`. Put `w=q(r)`.

Newton's identities for the last row imply that the sum of any real
polynomial of degree at most `2^m-1` over its roots is real: each of its
root power sums of orders `1,...,2^m-1` is real. Applying this fact to
the real polynomial `q(t)^k` gives the necessary conditions

```
Im(sum_(T nonempty) z_T^k + w^k) = 0,
1 <= k <= 2^(m-1)-1.                                  (1)
```

These conditions hold for every possible choice of one preimage from
each quadratic fiber. No choice of branches, coefficients of `q`, or
private root can bypass them.

In particular define

```
A = Im sum_T z_T,
B = Im sum_T z_T^2,
C = Im sum_T z_T^3,
E = 4 A C - 3 B^2 + 4 A^4.                            (2)
```

Writing `w=x+iy`, the first three equations in (1) are

```
A+y=0,
B+2xy=0,
C+3x^2 y-y^3=0.                                      (3)
```

If `A!=0`, these force `y=-A`, `x=B/(2A)`, and then `E=0`.
If `A=0`, they force `y=B=C=0`, which again implies `E=0`.
Consequently:

> If the root-union invariant `E` is nonzero, no full degree-one
> `(m+1)`-row profile can extend any real quadratic pullback of this
> full degree-one `m`-row profile.

The converse is not claimed. Even `E=0` only passes this necessary
three-moment test; the remaining equations in (1), the odd moments of
the last row, and all nondegeneracy conditions remain to be satisfied.

## 2. Reanchoring, affine invariance, and conjugation

Include the original anchor as label `0`, so that the profile has `m+1`
labels and its roots represent cuts oriented away from label `0`.
Reanchoring at a row label `a` replaces exactly the roots `z_T` with
`a in T` by their conjugates, then relabels each cut so that its chosen
side omits the new anchor. This relabeling permutes the nonempty cut
blocks. It gives a full polynomial profile: its phase rows are the old
pair quotients `B_j/B_a` and the old-anchor row `1/B_a`, whose reduced
numerators have constant imaginary part by the exact pair identities in
[boolean_polynomial_map_geometry.md](boolean_polynomial_map_geometry.md).

For `S_k=sum_T z_T^k`, conjugating the incident roots changes the global
power sum by exactly

```
S_k_new-S_k = -2i Im(sum_(T contains a) z_T^k) = 0,
1 <= k <= 2^(m-1)-1.                                  (4)
```

The last equality is the existing row moment equation. In particular
`A,B,C`, and therefore `E`, are unchanged before any affine
renormalization. Anchor changes and row permutations generate the full
permutation action on the `m+1` labels. Thus `E` is invariant under all
these reanchorings. For the certified four-row profile, (4) holds through
`k=7` and covers the entire `S_5` action; it is an exact consequence of
the certified equations, not a numerical comparison of anchor choices.

Replacing every base root by `a z_T+b`, with `a,b` real and `a!=0`, gives

```
A_new = a A,
B_new = a^2 B + 2ab A,
C_new = a^3 C + 3a^2 b B + 3ab^2 A,
E_new = a^4 E.                                       (5)
```

Global conjugation negates `A,B,C` and leaves `E` unchanged. Row
permutation also leaves the root union unchanged. Thus nonvanishing of
`E` is not an artifact of the anchor or normalization of the certified
profile. A normalization after reanchoring only multiplies `E` by a
nonzero real fourth power, so it cannot remove the obstruction.

## 3. Exact application to the certified four-row profile

Use the rational center `z_T^0` and radius `rho=10^(-30)` from
[four_row_real_polynomial_certificate.json](four_row_real_polynomial_certificate.json).
Its certified algebraic root tuple `z_T^*` lies in the coordinate cube
of radius `rho`; the full root is fixed at `i`. Compute `A_0,B_0,C_0,E_0`
by (2) at that rational center. For orientation only,

```
A_0 =  1.921619906674024...
B_0 =  1.223723176180431...
C_0 = -0.983810868143956...
E_0 = 42.487322991849254...
```

The exact certificate already places every root in `|Re z|+|Im z|<2`
throughout the cube. Only fourteen roots move. The real-coordinate
gradient bounds for the three moment summands are respectively

```
1,   2(|x|+|y|)<4,   3(|x|+|y|)^2<12.
```

The mean value theorem therefore gives

```
|A-A_0| <= 14 rho,
|B-B_0| <= 56 rho,
|C-C_0| <= 168 rho.                                   (6)
```

Exact rational comparisons give `|A|,|B|,|C|<3` on the cube.
For `E(A,B,C)=4AC-3B^2+4A^4`, the absolute partial derivatives
on that box are at most `444,18,12`. Hence

```
|E-E_0| <= (444*14+18*56+12*168) rho = 9240 rho,
E(z_T^*) >= E_0-9240 rho > 42.                        (7)
```

Every comparison in (6)--(7) is checked using exact rational arithmetic
by [check_quadratic_lift_five_row_moment_obstruction.py](check_quadratic_lift_five_row_moment_obstruction.py).
The checker first runs the existing contraction and nondegeneracy
certificate, then verifies the new invariant bound. Thus the conclusion
concerns the actual algebraic solution, not its rounded displayed roots.

## 4. The entire vertical three-row family also fails the test

The same obstruction applies to the rational family in
[the vertical triangle construction](quartic_root_completion_and_quadratic_split_criterion.md).
Its seven roots are

```
i;
X+i(s-1),        X(1+1/s)+i s,       s in {p,q,r},
p+q+r=1,        pqr=X^2,            X!=0.
```

Assume the nondegeneracy conditions of that construction. Directly summing
the imaginary coordinates gives

```
A=1+(p+q+r-3)+(p+q+r)=0.
```

The full root contributes zero to the imaginary second moment. The pair
and private roots give, respectively, `2X(sum s-3)` and
`2X(sum s+3)`. Hence

```
B=4X,       E=-3B^2=-48X^2!=0.
```

Thus no member of this entire family admits a full four-row degree-one
extension through any common real quadratic pullback. By Section 2 this
also holds after reanchoring, real affine changes and conjugation. This
differs from the earlier rational splitting obstruction for one seed:
the present conclusion permits arbitrary real preimages but requires a
new row with constant imaginary part. It does not assert that the
quadratic factors themselves cannot split over `Q(i)`.

## 5. Remaining scope

A general full five-row degree-one profile restricts to a full four-row
profile with degree-two blocks after forgetting one label. Equal block
degrees do not imply that those blocks are fibers of one real quadratic
polynomial. The missing statement needed to promote this result to a
general five-row exclusion would be a classification forcing such a
common composition, together with a restriction on its base profile.
Neither statement is established here. The result also supplies no
converse from integer endpoint configurations to polynomial profiles.
