# The remaining common rotation in an integral norm realization

The norm and determinant equations determine a configuration only up to
a common rotation. This note classifies all its integral rotations exactly.
For the full Boolean profile, the possible common factors have the size of
one norm block. The small-imaginary-coordinate requirement selects at most
one such rotation, but the argument does not exclude its existence and
does not bound the number of rows.

## 1. The norm of the common Gaussian divisor is a rational gcd

Let `P_1,...,P_m` be nonzero Gaussian integers, and write

```text
N_i=Norm(P_i),
Delta_ij=Im(bar(P_i)P_j),
S_ij=Re(bar(P_i)P_j),
H=gcd_Z[i](P_1,...,P_m),       Q_i=P_i/H.
```

Then the following identity is exact:

```text
Norm(H)=gcd_Z({N_i}_i union {Delta_ij}_{i<j}).            (1)
```

To prove it, set `d=Norm(H)`. The `Q_i` generate the unit ideal in
`Z[i]`. Choose Gaussian integers `a_i` with `sum a_i Q_i=1`.
Multiplying this identity by its conjugate shows that the products
`bar(Q_i)Q_j` also generate the unit ideal. Consequently all the
Gram entries `bar(P_i)P_j` generate exactly the ideal `(d)`.

Certainly `d` divides every `N_i` and `Delta_ij`. Conversely let `d_0`
be their ordinary positive gcd. The identity

```text
S_ij^2=N_i N_j-Delta_ij^2
```

implies `d_0 | S_ij`: this follows prime by prime from divisibility of
the right side by `d_0^2`. Thus `d_0` divides every Gaussian Gram
entry and hence divides their Gaussian linear combination `d`.
As both are ordinary integers, this proves (1).

## 2. All integral realizations of the fixed oriented Gram matrix

Assume some `Delta_ij` is nonzero. Fix the `P_i`, and hence the
Gaussian-coprime tuple `Q_i`, as above. Every other Gaussian-integer
tuple with the same `N_i`, `S_ij`, and `Delta_ij` is exactly

```text
P_i'=H' Q_i,          H' in Z[i],       Norm(H')=d.       (2)
```

Indeed, two spanning real configurations with the same Gram matrix
and orientation differ by an orientation-preserving orthogonal map.
In complex notation this is multiplication by `zeta` with
`|zeta|=1`. Since one nonzero `P_i` and its image are Gaussian
integers, `zeta=P_i'/P_i` belongs to `Q(i)`. If every `zeta H Q_i`
is integral, Gaussian Bezout gives

```text
zeta H=sum_i a_i P_i' in Z[i].
```

Thus `H'=zeta H` has norm `d`. The converse in (2) is immediate.
This classification fixes the full Gram matrix, not only the diagonal
norms and unsigned pair areas. Three noncollinear labels determine its
rational off-diagonal entries when the norm and oriented determinant
data satisfy the compatibility equations.

In particular there are exactly `r_2(d)` labelled integral realizations
of this fixed oriented Gram matrix, where `r_2` counts signed ordered
representations as a sum of two squares. Dividing out by a common unit
gives `r_2(d)/4` possibilities.

## 3. Only one norm block is available for this rotation

Suppose the realized data have the full Boolean factorizations

```text
N_i=k_i product_(T containing i) n_T,
Delta_ij=t_ij product_(T containing i,j) n_T,
```

with positive integers `k_i,n_T`, nonzero integers `t_ij`, and
pairwise disjoint rational-prime supports of the `n_T`. No Gaussian
orientation hypothesis is needed for the following conclusion. Put
`g=n_[m]`. Equation (1) gives

```text
g | d,                d/g | product_i k_i.              (3)
```

For the second divisibility, consider each prime separately. At a
prime of the full-set block, the exponent in `d` is at most the
common-block exponent plus the minimum correcting exponent. At a
prime of a proper block `T`, choose any row outside `T`; its norm
contains that prime only through `k_i`, so it bounds the exponent
in `d`. A prime absent from every block is handled by any row.
Each bound is at most the exponent in `product k_i`, proving (3).

Under the extracted equal-block profile,

```text
log n_T=w+o(w),       log k_i=o(w),
```

one therefore has `log d=w+o(w)`. The number of possible rotations
is subpower, since `r_2(d)<=4 tau(d)=exp(o(w))`. This is a count for
one fixed Gram configuration; it does not bound the number of possible
Gram configurations or the number of points in the original arc.

## 4. Tiny imaginary coordinates select a unique rotation

Consider two tuples in (2) with every real coordinate positive and
every imaginary coordinate of absolute value at most `B`. For any
fixed row `i`, if

```text
4 d B^2 < N_i,                                          (4)
```

the tuples coincide. On the circle of radius `sqrt(N_i)`, the cap
with positive real coordinate and `|Im|<=B` has diameter at most
`2B`. Indeed its arguments lie between
`-arcsin(B/sqrt(N_i))` and `arcsin(B/sqrt(N_i))`, and the longest
chord is the one joining the two endpoints. Thus

```text
|H-H'|=|P_i-P_i'|/|Q_i| <= 2B sqrt(d/N_i) < 1.
```

Distinct Gaussian integers have distance at least one, so `H'=H`.
Here (4) ensures `B<sqrt(N_i)`, as needed for the displayed angular
interval. This chord argument improves the weaker bound
`sqrt(N_i)>2Bd` obtained by using an integer determinant between
`H,H'`. Neither bound gives a nonexistence theorem.

For `m>=2` in the full profile, `sqrt(N_i)=exp(2^(m-2)w+o(w))`,
`d=exp(w+o(w))`, and the desired `B=exp(o(w))` satisfy (4).
The common rotation cannot be varied freely while keeping the
endpoint imaginary-coordinate bound. Nevertheless uniqueness of a
possible realization is not nonexistence. In fact, the
[balanced two-row family](small_imaginary_row_lattice_transference.md)
and the [balanced three-row family](balanced_three_row_gaussian_polynomial_family.md)
have such rotations at arbitrarily large block heights. An exclusion
of sufficiently large balanced configurations would therefore have to
use some fixed `m>=4` or additional hypotheses not captured here.

The [exact checker](check_gaussian_common_rotation_classification.py)
checks the gcd identity, rotation classification on small fixtures,
and the displayed finite endpoint-profile configuration. The general
claims above follow from the proofs, not from finite enumeration.
