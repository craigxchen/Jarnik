# Exterior coefficient floors: valid transfer and a strict-gain barrier

For a rank-`h` coherent kernel, one can retain all rank constraints at
a coalition by taking an `h`-by-`h` determinant of raw source-monomial
coefficients. The resulting generic valuation gives an actual arithmetic
divisibility estimate, with the residual-discriminant index loss. However,
its full Boolean sum never exceeds the existing fusion covolume upper
exponent. In the smallest higher-rank example `(n,d)=(2,3)`, the minimum
over nonzero raw-coordinate determinants is exactly the sum of the four
individual coefficient floors.

## 1. Exterior coordinates and their generic orders

Fix `n,d>=2`, `m=nd`. Let `M_Z` be the integral multilinear invariant
lattice in a fixed integral basis, and let `K` be the
[higher-rank coherent kernel](higher_rank_terminal_fusion_grid.md) of
generic binary directions. Write `h=rank K`. The
[primitive polynomial determinant tuple](terminal_kernel_determinant_line_height.md)
is denoted by `W`, with the common binary row degree `delta` and affine
total degree `E=m*delta/2`. Normalize its integer polynomial coefficients
to have gcd one. Every component is an `SL_2` invariant.

The constant integral injection from `M_Z` into the full source-monomial
coefficient lattice induces an injection on `h`th exterior powers. For an
`h`-element set `S` of raw monomial coordinates, let `F_S` be the resulting
scalar polynomial component of `W`. Every nonzero `F_S` remains an
`SL_2` invariant of row degree `delta`. A constant injective linear map
preserves the minimum order of a tuple over the characteristic-zero
function field, so the common coalition order is still the fusion order
`nu_T=nu_(|T|)`. Set

```text
r_(S,T)=ord_T F_S-nu_T >= 0.                              (1)
```

Zero polynomial components are excluded from (1). The orders include
empty, singleton, full, and near-full coalitions with the affine
conventions in the
[invariant full Boolean order cap](invariant_full_boolean_coalition_order_cap.md).
In particular, `r_(S,empty)=r_(S,[m])=0`.

## 2. A coefficientwise transfer to the actual saturated lattice

First use the universal affine collision parameters

```text
z_i=s*a_i  (i in T),       z_j=b_j  (j outside T),
R_0=Z[s,(a_i),(b_j)].
```

Let `G_T` be the matrix of all local determinant generators after division
by the support-dependent powers `s^max(0,|I intersection T|-n)`. Its
entries lie in `R_0`. Put

```text
W_T=s^(-nu_T) W((1,z_i)_i).                              (2)
```

This is an integral polynomial tuple and has gcd one in `R_0`.
Indeed, any common affine factor would rehomogenize to a common
multihomogeneous factor of `W`. Dehomogenization preserves integer
coefficient content too, since the exponent of one binary coordinate
determines the other from the fixed row degree. After inverting `s`, the parameter change is the
isomorphism `a_i=z_i/s`, so it also creates no gcd. Before inverting `s`
any new irreducible common factor must therefore be `s`. The exact common
order `nu_T` removes that factor. Substitution sends distinct affine
monomials to distinct `(a,b)` monomials, so it creates no common integer
content either.

For every `h`-element set `J` of generator columns, their wedge is
collinear with `W_T` over the fraction field. Hence

```text
exterior^h G_(T,J)=c_J W_T,          c_J in R_0.           (3)
```

To justify polynomiality, a denominator of the rational scalar `c_J`
would divide every coordinate of the primitive tuple `W_T` in the UFD
`R_0`. There is no such denominator. Applying the raw-coordinate
projection `S` to (3) proves the polynomial divisibility

```text
s^r_(S,T) divides det P_S G_(T,J)   for every J.          (4)
```

Thus (4) holds coefficientwise, before choosing residual parameters. It
does not assume generic actual residues.

Now specialize to a characteristic-zero DVR `R`, with `s!=0`, distinct
generic-fiber directions and nonzero residual discriminant `D_T(a,b)`.
Let `N_R` be the generated lattice and `K_R` its saturation. Choose a
basis matrix `B_K` of `K_R` and write `G_T=B_K C`. If
`ell=length_R(K_R/N_R)`, the ideal of all `h`-minors of `C` is
`(pi^ell)`. Therefore, for the fixed raw-coordinate set `S`,

```text
ideal_J(det P_S G_(T,J))=(pi^ell det P_S B_K).           (5)
```

Write `w_S=det P_S B_K`. Equations (4)--(5) and the
[residual-discriminant local index bound](residual_discriminant_local_index_bound.md)
give

```text
v_R(w_S) >= r_(S,T) v_R(s)-ell
         >= r_(S,T) v_R(s)-v_R(A)-H v_R(D_T).           (6)
```

The zero coordinate has infinite valuation, so causes no exception.
Equation (5) uses the ideal of all generator minors; no selected large
minor is substituted for a lattice covolume. Also, (6) does not assume
that the evaluated global tuple has exactly its generic common content.

There is a uniform global version under the hypotheses of the
[primitive local-generator index theorem](global_primitive_local_generator_index.md).
Use its notation

```text
Delta_ij=t_ij product_(T containing i,j) n_T,
B_res=product_(i<j)|t_ij|,       Gamma=product_i gcd(|X_i|,|Y_i|).
```

The `n_T` have disjoint rational-prime support. At a core prime of depth
`e`, that theorem constructs the chart with
`v(s)=e'=max(0,e-2u)`, where `u` is the maximum row-content valuation,
and

```text
v(D_T)<=v(B_res)+m(m-1)v(Gamma).
```

Since `0<=r_(S,T)<=E` and `u<=v(Gamma)`, (6) implies

```text
product_T n_T^r_(S,T)
    divides A B_res^H Gamma^(H*m*(m-1)+2E) w_S.        (7)
```

Here `w_S` is a coordinate of the primitive integer basis wedge of the
actual saturated coherent kernel. A simultaneous binary basis change
does not change the kernel; individual row rescalings only multiply its
polynomial determinant tuple by a common scalar. Consequently the chart
changes used in deriving (7) preserve these actual coefficient
coordinates. All exponents on the right depend only on `n,d`.

For nonzero `w_S`, the full Boolean profile with `log n_T=w+o(w)`,
`log B_res=o(w)`, and `log Gamma=o(w)` therefore gives the rigorous floor

```text
log |w_S| >= (sum_T r_(S,T)) w-o(w).                    (8)
```

## 3. The generic floor cannot exceed the fusion upper exponent

The [sharp invariant order cap](invariant_full_boolean_coalition_order_cap.md)
applies separately to every nonzero scalar polynomial `F_S` and proves

```text
sum_T ord_T F_S <= delta*m*2^(m-3).
```

Subtracting the common fusion orders yields, for every such `S`,

```text
sum_T r_(S,T)
 <= delta*m*2^(m-3)-sum_T nu_T
 = B_fusion(n,d).                                      (9)
```

Thus even the largest generic exterior-coordinate floor is at most the
existing covolume upper exponent; the minimum over all nonzero
coordinates certainly cannot be strictly larger. The fixed conversion
between an invariant basis and raw monomial coordinates changes
Archimedean norms only by constants, so it does not alter this exponent
comparison.

Some exact parameters are:

| `(n,d)` | kernel rank | `delta` | raw degree budget | common fusion content | `B_fusion` |
|:---:|---:|---:|---:|---:|---:|
| `(2,3)` | 4 | 1 | 48 | 30 | 18 |
| `(3,3)` | 21 | 8 | 4,608 | 2,520 | 2,088 |
| `(3,4)` | 341 | 55 | 337,920 | 166,320 | 171,600 |

In particular no generic coalition-order computation for the twelve-row
kernel can make `min_S sum_T r_(S,T)>171600`. The obstruction applies to
arbitrary constant linear combinations of determinant coordinates too,
because a surviving linear combination is still a binary invariant of
the same row degree.

## 4. Exact raw-coordinate enumeration at `(n,d)=(2,3)`

There are five independent multilinear binary invariants on six rows.
Use the noncrossing matching basis

```text
f_1=[12][34][56],       f_2=[12][36][45],
f_3=[14][23][56],       f_4=[16][23][45],
f_5=[16][25][34].
```

The coherent map is their scalar evaluation, so its kernel has rank four.
Its primitive source determinant tuple is the signed dual of
`(f_1,...,f_5)`. The raw coefficient embedding has 20 rows, indexed by
the three-element sets receiving the second binary coordinate.

The [exact checker](check_exterior_coefficient_coalition_floor_barrier.py)
enumerates all `binom(20,4)=4845` raw-coordinate projections. Of these,
`1725` polynomial coordinates vanish identically. The remaining counts
are

| `sum_T r_(S,T)` | number of raw-coordinate sets `S` |
|---:|---:|
| 8 | 240 |
| 18 | 2,880 |

Consequently the minimum is exactly `8`, equal to four times the scalar
coefficient floor `M(2,3)=2`. This is an exact polynomial calculation,
not a random finite-field specialization.

One witness for the minimum is

```text
S={{1,2,3},{1,2,4},{1,3,5},{1,4,5}},
F_S=f_3-f_2=[14][23][56]-[12][36][45].
```

For this coordinate, `r_(S,T)=1` exactly when `T` is one of those four
three-element sets or its complement, and is zero otherwise. The
common order is `nu_T=max(0,|T|-3)`.

The result rules out a strict gain from these generic determinant
valuations under the uniform full Boolean weighting. It does not cap
additional valuations at special arithmetic residual tuples, prove the
existence of an actual circle configuration, or settle the radius-uniform
arc bound. A sufficient new determinant input would have to control
extra actual divisibility beyond (8), or improve the actual upper
exponent below the applicable generic floor.
