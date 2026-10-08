# A parity gap for four thin signed witnesses on eight pattern groups

The [four-witness matching-norm gap](four_signed_witness_matching_norm_gap.md)
excludes a closed four-sign pattern whose three changing groups all have
large norms. There is a quantitative extension for a quartet of pairwise
orthogonal sign words even when their product is not constant. It bounds the
imbalance between the even and odd *fourth-order* sign patterns. The bound
does not control the number of sign words independently of the number of
blocks.

## 1. Exact eight-group divisor lemma

Let `E_0,...,E_3,O_0,...,O_3` be nonzero conjugate-primitive Gaussian
integers whose positive integer norms are pairwise coprime. Units are
allowed. Put `e_j=N(E_j)`, `o_j=N(O_j)`, and
`P=prod_j e_j o_j`. Form four Gaussian products `B_0,...,B_3` by using
each of the eight groups once, with the following orientations. A `+`
means the group itself and a `-` means its conjugate.

| group | `B_0` | `B_1` | `B_2` | `B_3` |
|---|:---:|:---:|:---:|:---:|
| `E_0` | + | + | + | + |
| `E_1` | + | + | - | - |
| `E_2` | + | - | + | - |
| `E_3` | + | - | - | + |
| `O_0` | + | + | + | - |
| `O_1` | + | + | - | + |
| `O_2` | + | - | + | + |
| `O_3` | + | - | - | - |

The labels `E` and `O` refer to the parity of the three signs after
`B_0`: their products are respectively `+` and `-`. Every `B_i` has norm
`P` and is conjugate-primitive. For nonzero Gaussian corrections `K_i`,
write

```text
U_i=K_i B_i=X_i+iY_i,
M=max_i |K_i|,                 Y=max_i |Y_i|,
Delta_ij=X_i Y_j-X_j Y_i.
```

The common oriented Gaussian factor of each pair has the following norm.
The two factors in each row are *opposite edges* of the quartet.

| edges | common-factor norms `G_ij` |
|---|---|
| `01`, `23` | `G_01=e_0 e_1 o_0 o_1`; `G_23=e_0 e_1 o_2 o_3` |
| `02`, `13` | `G_02=e_0 e_2 o_0 o_2`; `G_13=e_0 e_2 o_1 o_3` |
| `03`, `12` | `G_03=e_0 e_3 o_1 o_2`; `G_12=e_0 e_3 o_0 o_3` |

Thus `G_ij | Delta_ij` in the ordinary integers, including when the
corrections overlap core factors. Define

```text
Q_*=min_(i<j) sqrt(P/G_ij),
O=max(o_0,o_1,o_2,o_3),       G=min_(i<j) G_ij.
```

If `M<Q_*`, every `Y_i` and `Delta_ij` is nonzero. Hence the six numbers
`t_ij=Delta_ij/G_ij` are nonzero integers. With
`T=max_(i<j)|t_ij|`, they satisfy the **exact finite bounds**

```text
max(e_1,e_2,e_3) <= 2 Y^2 T^2 O^2,
T <= 2 M sqrt(P) Y/G,
max(e_1,e_2,e_3) <= 8 M^2 P Y^4 O^2/G^2.          (1)
```

The first two inequalities retain more information than the last one.
The same first inequality applies whenever all the `Y_i` and
`Delta_ij` are nonzero, without the sufficient cutoff `M<Q_*`.

**Nonvanishing and divisibility.** If `U_i` were real, then
`K_i B_i=bar(K_i) bar(B_i)` and conjugate-primitivity would give
`B_i | bar(K_i)`. Therefore `M>=sqrt(P)>=Q_*`. For a pair `i,j`, write
their cores as `B_i=F E` and `B_j=F bar(E)`, where `F` is their common
oriented factor and `N(E)=P/G_ij`. Collinearity of `U_i,U_j` would give
`E^2 K_i bar(K_j)=bar(E)^2 bar(K_i)K_j`. Since `E` and `bar(E)` are
coprime, `E^2 | bar(K_i)K_j`, whence
`P/G_ij=N(E)<=|K_i||K_j|<=M^2`. This contradicts `M<Q_*`.
Finally `F|U_i,U_j`, so `N(F)=G_ij` divides their determinant.

**Rank-two proof of (1).** For any three indices,
`Y_k Delta_ij-Y_j Delta_ik+Y_i Delta_jk=0`. Divide these four
identities by their common factor `e_0 o_0`, `e_0 o_1`, `e_0 o_2`,
`e_0 o_3`, respectively. The resulting integer rows annihilate the
primitive vector `(e_1,e_2,e_3)^t`:

```text
R_012=(Y_2 t_01 o_1, -Y_1 t_02 o_2,  Y_0 t_12 o_3),
R_013=(Y_3 t_01 o_0,  Y_0 t_13 o_3, -Y_1 t_03 o_2),
R_023=(Y_0 t_23 o_3,  Y_3 t_02 o_0, -Y_2 t_03 o_1),
R_123=(Y_1 t_23 o_2, -Y_2 t_13 o_1,  Y_3 t_12 o_0).       (2)
```

The first row is the `012` identity, and similarly for the others.
For a rank check, divide each row by the product of its three nonzero
`Y` coordinates and abbreviate

```text
a=t_01/(Y_0Y_1), b=t_02/(Y_0Y_2), c=t_03/(Y_0Y_3),
d=t_12/(Y_1Y_2), e=t_13/(Y_1Y_3), f=t_23/(Y_2Y_3).
```

The rows then become

```text
(a o_1,-b o_2,d o_3), (a o_0,e o_3,-c o_2),
(f o_3,b o_0,-c o_1), (f o_2,-e o_1,d o_0).               (3)
```

All six letters and all four `o_j` are nonzero. If (3) had rank one,
proportionality of its first two rows would force
`c=-d o_0 o_3/(o_1 o_2)`. Proportionality of its first and third rows
would force `c=+d o_0 o_3/(o_1 o_2)`, a contradiction. Its rank is
therefore at least two; the nonzero null vector makes the rank exactly
two. A cross product of two independent integer rows in (2) is a
nonzero integer multiple of `(e_1,e_2,e_3)`, because these three norms
are pairwise coprime. Each row entry has modulus at most `YTO`, so each
cross-product entry is at most `2Y^2T^2O^2`. This proves the first
inequality. The determinant bound
`|Delta_ij|<=2M sqrt(P)Y` proves the remaining inequalities.

## 2. Consequence for balanced orthogonal sign words

Fix an integer `r` and let `H_1,...,H_r` be conjugate-primitive,
rational-prime-disjoint Gaussian blocks with

```text
log N(H_j)=w+o(w)                 (j=1,...,r),
```

as `w` tends to infinity. Suppose four signed products of all `r`
blocks have pairwise Hamming distance `r/2`, and Gaussian corrections
make all four imaginary coordinates and all corrections subpower:
`log max(1,M)=o(w)` and `log max(1,Y)=o(w)`.

Reorient each block so the first word has sign `+`. Pairwise
orthogonality makes the remaining three signs unbiased and pairwise
independent over the `r` block positions. Walsh inversion on three
bits shows that every even-parity triple pattern occurs the same
number `alpha` of times, and every odd-parity pattern the same number
`beta` of times. In particular `r=4(alpha+beta)`. Grouping the
corresponding blocks produces exactly the table in Section 1, with

```text
log e_j=alpha w+o(w),       log o_j=beta w+o(w).
```

Each pair agrees on `r/2` blocks. Thus `log G_ij=rw/2+o(w)`,
`log P=rw+o(w)`, and `log Q_*=rw/4+o(w)`. The cutoff holds for large
`w`, while (1) gives `log T=o(w)`. Taking logarithms of the first
inequality in (1) yields `alpha<=2 beta`. Conjugating one entire
witness preserves its correction height, absolute imaginary height,
pairwise orthogonality, and all nonvanishing conditions. Its effect on
the product of the four sign words is to exchange even and odd pattern
classes. Applying the same argument again gives `beta<=2 alpha`.
Therefore

```text
|alpha-beta|/(alpha+beta) <= 1/3.                    (4)
```

The left side is the absolute normalized fourth-order sign correlation
`|(1/r) sum_j sigma_0(j)sigma_1(j)sigma_2(j)sigma_3(j)|`.
In particular, an XOR-closed quartet (one of `alpha,beta` equal to zero)
cannot have all four corrected witnesses thin under these assumptions.
For `r=8` or `r=16`, integrality and (4) force `alpha=beta` for every
quartet. Then all pairwise products of distinct sign words, together
with the constant word, are orthogonal, giving
`1+binom(k,2)<=r` for any `k` such thin words on these blocks.

The dependence on `r` cannot be eliminated by this sign constraint.
For `r=2^m`, index the block positions by `x in F_2^m` and take the
`m` sign words

```text
sigma_i(x)=(-1)^(x_i),                 i=1,...,m.
```

They are pairwise orthogonal, and every nonempty product of distinct
words is a nonconstant Walsh character. Hence they have no closed XOR
circuit of **any** size; every quartet has `alpha=beta=r/8` and
fourth-order correlation zero. This is an abstract sign-code example,
not a construction of thin Gaussian witnesses. It shows that pairwise
orthogonality and exclusions of closed four-, six-, or higher-order
sign circuits alone cannot give a radius-uniform count. The endpoint
profile also has differing sparse supports, so applying this common
signed-support lemma to its actual rows requires a separate arithmetic
construction.

The [exact checker](check_four_signed_witness_parity_gap.py) tests the
Gaussian sign table, determinant divisibilities, rank-two integer kernel,
finite bounds, and the Walsh sign counts. These finite fixtures supplement
the proof; they do not construct thin Gaussian witnesses or prove a
uniform cluster bound.
