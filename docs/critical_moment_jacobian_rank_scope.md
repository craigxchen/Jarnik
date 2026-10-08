# Generic full Jacobian rank for the admissible Euler words

Every word in the existing admissible Euler family has a full-column-rank
moment Jacobian for generic magnitude nodes. The proof explicitly partitions
its columns into independent sets and identifies a nonzero determinant
coefficient. This does not exclude solutions at exceptional nodes. A small
Boolean example shows why the specialization step still requires proof.

## The Jacobian and its signs

Let `S` have eight sign rows and `n=4s+2` columns. Set
`B_ij=(S_ij-S_0j)/2`, for `1<=i<=7`. The full common-row odd-moment
equations are

```
B a^(2r+1)=0,       0<=r<s,       0<a_1<...<a_n.
```

After dividing its degree-`r` row block by `2r+1`, their Jacobian is

```
J_s(x)=[B; B diag(x); ...; B diag(x^(s-1))],
x_j=a_j^2.
```

Every full moment solution satisfies `J_s(a^2)a=0`. Full column rank at
every ordered positive node tuple would therefore exclude a word.

Every nonconstant column has `B_j=epsilon_j c_j`, with `epsilon_j` a
sign and `c_j` a nonzero Boolean vector of length seven. Complementary
eight-bit masks have the same `c_j`. Normalizing columns to `c_j`
preserves rank, but changes the required kernel vector to
`(epsilon_j a_j)_j`. The signs are not discarded from the moment problem.

## An explicit partition and the actual-family budget

Represent nonzero seven-bit vectors by integers `1,...,127`. The binary
linear transformation

```
T(v)=((v<<1)&127) xor (3 if v&64 else 0)
```

has an orbit of length 127 starting at one. Its first seven iterates are
`1,2,4,8,16,32,64`. Thus every seven consecutive orbit vectors form a
basis over `F_2`, being the images of this basis under an invertible power
of `T`. Their zero-one integer columns are also real-independent: their
determinant is odd. Split the orbit into eighteen groups of seven and one
final singleton, obtaining nineteen independent groups.

If each Boolean direction occurs at most `m0` times, partition all
columns into `19m0` independent sets as follows. For occurrence number
`ell=0,...,m0-1` and group `g=0,...,18`, take the `ell`th occurrence of
each direction in group `g` that has such an occurrence. Each set is a
subset of an independent group.

The [admissible Euler words](critical_moment_compressed_graph_even_q2.md)
have `s=2+20,805,120 K`, `K>=1`, and no constant columns. The
[proved multiplicity bound](critical_moment_christoffel_relaxation_limit.md)
is `403,200 K+5` per signed eight-bit mask. A Boolean direction merges at
most two complementary masks. Consequently

```
m0 <= 806,400 K+10,
s-19m0 >= 5,483,520 K-188 > 0.                         (1)
```

Assign set `(ell,g)` to evaluation degree `r=19ell+g`. Equation (1)
places every used degree in `0,...,s-1`. This explicitly partitions the
actual word columns into independent sets, one per permitted degree.

## A nonzero maximal minor

For each degree's column set `I_r`, choose `|I_r|` rows `R_r` of the
original signed matrix `B` with `det B_(R_r,I_r)!=0`. Select these rows
from degree block `r` of `J_s(x)`, obtaining a square minor of size `n`.
In its row-block determinant expansion, the coefficient of

```
product_r product_(j in I_r) x_j^r
```

is, up to sign, `product_r det B_(R_r,I_r)`, which is nonzero. Indeed,
the exponent of each labeled column fixes its degree block, so no other
block partition contributes the same monomial. Within-block expansions
are exactly the indicated determinants. No matroid-partition theorem is
needed for this explicit argument.

Thus every actual Euler word has generic full Jacobian column rank.
Substituting `x_j=a_j^2` preserves nonzeroness, since it doubles exponent
vectors injectively. Full rank holds on a dense open subset of the
strictly ordered positive magnitude chamber. Any full moment solution
must lie in the simultaneous zero set of all maximal minors, including
the nonzero polynomial just constructed. No argument here excludes that
exceptional algebraic set from the ordered positive chamber.

## A Boolean rank drop at integer-square nodes

The specialization issue already occurs with distinct Boolean directions.
For `s=2`, take

```
C=[1 0 1 0 1 0]
  [0 1 0 0 1 1]
  [0 0 1 1 0 1],
a=(1,3,4,5,6,9),       x=(1,9,16,25,36,81).
```

The six columns are distinct and nonzero. With zero-based indices,
`{0,1,3}` and `{2,4,5}` are bases, with determinants one and two.
They give generic full column rank six by the same proof. Yet

```
lambda=(-32,15,56,-65,-24,9),
C lambda=C diag(x) lambda=0,
rank J_2(x)=5.                                       (2)
```

The exact determinant polynomial is

```
det J_2(x)=-[(x_2-x_0)(x_4-x_1)(x_5-x_3)
            +(x_4-x_0)(x_5-x_1)(x_2-x_3)].             (3)
```

The two bracketed terms at the square nodes are `22680` and `-22680`.
At `x=(1,2,3,4,5,6)`, the determinant instead equals four.

This is a Boolean Jacobian fixture, not an orthogonal eight-row template.
Moreover `lambda` is not proportional to `(epsilon_j a_j)_j` for any
choices of signs, so it is not a full moment solution. It shows only
that Boolean directions, no repetitions, a generic independent-set
partition, and distinct ordered positive square nodes do not together
give universal full rank. The special orthogonal word or the specific
signed radial kernel could still yield a stronger exclusion.

## Verification

The [exact checker](check_critical_moment_jacobian_rank_scope.py) verifies
the orbit, all 127 cyclic seven-vector bases, the nineteen-group partition,
and the all-`K` degree budget. It checks both Boolean ranks, the kernel,
and (3) on 84 ordered integer-node tuples. It does not enumerate the
giant Euler word or claim that its exceptional locus contains a solution.
