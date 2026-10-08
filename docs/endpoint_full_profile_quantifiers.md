# Endpoint extraction: full-profile quantifiers and determinant residuals

This is an independent audit and quantitative supplement to
[central truncation](endpoint_central_truncation.md),
[allocation rounding](endpoint_allocation_rounding.md), and
[profile extraction](uniform_profile_extraction.md). The extraction is
valid for every fixed number of rows. A sharper pair-residual estimate
below also preserves its full growing-dimension range. Products of all
residuals require a further polynomial factor in the dimension.

An alternative [exact nested-layer extraction](exact_nested_profile_residual_extraction.md)
now retains all prime powers with no correction factors and aggregate
residual height `O(k^2 W/M)`. Its blocks can share rational primes along
nested chains. The present note retains the prime-disjoint formulation
and its `O(W/sqrt(M))` loss; an arithmetic theorem requiring disjoint
supports cannot silently use the alternative's smaller height budget.

The input is the repository's centered-circle problem: distinct Gaussian
integers of common modulus `R` on an arc of length `C sqrt(R)`.
No reduction for arbitrary translated circle centers is asserted here.

## 1. Normalization and the original cardinality

To disprove a radius-uniform bound for a fixed `C>0` would be to have
such clusters with cardinalities tending to infinity. Subdivide the arcs
into at most `ceil(C/C_0)` pieces for any fixed
`0<C_0<min(C,sqrt(2))`; one piece still has unbounded cardinality.
Call its cardinality `M`.

Remove this entire cluster's Gaussian gcd. Division preserves Gaussian
integrality and distinctness. If its modulus is `d>=1`, the radius becomes
`R/d`, while the endpoint constant becomes at most `C_0/sqrt(d)<=C_0`.
Every inert and ramified prime in the common norm is removed. Rename the
remaining radius `R`; its square `N=R^2` is odd and split-supported. For unbounded `M`
we have `N>1`. Write

```text
W=log N,       lambda=log(sqrt(2)/C_0)>0.
```

The primitive equal-norm pair difference has squared modulus at least
two. Thus the allocation distances obey

```text
d_ij >= W/2+2lambda.
```

The threshold sign vectors have pairwise inner products at most
`-4lambda/W`. Nonnegativity of the squared norm of their sum gives

```text
W >=4(M-1)lambda.                                       (1)
```

This bound applies before selecting any smaller tuple. It is the
reason block heights cannot remain bounded when the original cluster
cardinality tends to infinity.

The audited interior-allocation bound is
`sum_i E_i<=W sqrt(M)/2`. At least `M/2` rows therefore have
`E_i<=W/sqrt(M)`. Their largest common Gaussian-unit class has
size `L>=M/8`. No further common factor is removed at this stage.
In particular, all subsequent weights still use this same `W`.

## 2. Actual primitive rows, with an explicit common scale

Choose `k=m+1` rows, including anchor `0`, by simultaneous pattern
and pair-moment extraction. For `k<=L`, its exact bounds are

```text
eta=2(2^k-1)/sqrt(L-k+1),
zeta=2sqrt(binom(k,2)/(L-1)),
W/2<=d_ij<=(1+zeta)W/2.
```

Central truncation removes total threshold mass `2s`, where
`s<=kW/sqrt(M)`. Put

```text
b=2^m,                    w=W/b.
```

This note consistently uses the original weight `w`, not the weight
after empty-block removal. The retained blocks `H_T`, including
`T=empty`, satisfy

```text
|log Norm(H_T)-w| <= epsilon w,
epsilon=eta+2b s/W.                                     (2)
```

For (2), the inherited merged cut has weight within `eta w` of
`w`, and at most `2s` of its original threshold layers are removed.
This is the sharper old-scale bound already available from the literal
retained-layer interpretation of central truncation.

For each nonanchor row take its actual primitive anchor numerator

```text
P_i=product_p pi_p^((a_ip-a_0p)_+)
                  bar(pi_p)^((a_0p-a_ip)_+)=X_i+iY_i.
```

The common unit cancels exactly, so
`z_i/z_0=P_i/bar(P_i)`. The central construction gives

```text
P_i=K_i product_(T containing i)H_T,
gcd_G(P_i,bar(P_i))=1,
log|K_i|<=s,
0<|Y_i|<=(C_0/2)exp(zeta W/4).                           (3)
```

Each block is conjugate-primitive; different blocks have disjoint
rational-prime support. The correcting factors need not be coprime to
the incident core blocks. Their compatibility with the opposite
orientations follows from the exact primitivity in (3).

The empty block appears in none of the prescribed core row products
and can now be omitted. Equations (2)--(3) still hold for all `2^m-1` nonempty
blocks. Its omission does not reset the angular precision.

Distinct original points also give

```text
Delta_ij=Im(bar(P_i)P_j)!=0.
```

Indeed a zero determinant makes `P_i/P_j` real, hence makes their
conjugate quotients equal and forces `z_i=z_j`. A nonzero imaginary
part in (3) follows in the same way from `z_i!=z_0`.

## 3. A sharper bound for every pair residual

Set

```text
G_ij=product_(T containing i,j)Norm(H_T),
t_ij=Delta_ij/G_ij in Z\{0}.
```

The following estimate uses the actual extracted rows:

```text
log|t_ij| <= log(C_0/2)+zeta W/4+2s
          <= zeta W/4+2s.                               (4)
```

To prove it, let `D_ij` be a Gaussian gcd of `P_i,P_j`. Primewise,
the two anchor differences have a common oriented power exactly on
the threshold layers where both rows differ from the anchor. Therefore

```text
log Norm(D_ij)=(d_0i+d_0j-d_ij)/2.                       (5)
```

The core common factor keeps exactly those layers that survive central
truncation. Hence

```text
G_ij divides Norm(D_ij),
0<=log(Norm(D_ij)/G_ij)<=2s.                             (6)
```

All these are exact integer divisibilities. Removing the common
Gaussian gcd from the two rows gives the nonzero integer
`r_ij=Delta_ij/Norm(D_ij)`. If `theta_ij` is the original point
angle difference, then

```text
|r_ij|=exp(d_ij/2)|sin(theta_ij/2)|
      <=(C_0/2)exp((d_ij-W/2)/2)
      <=(C_0/2)exp(zeta W/4).
```

Here (5) cancels the two anchor distances exactly. Multiplying by
`Norm(D_ij)/G_ij` and using (6) proves (4).

The generic block-error estimate for `t_ij` remains valid, but it
would introduce an avoidable extra factor `2^m` when `m` grows.
Equation (4) avoids that loss. It does not require constant imaginary
coordinates or a new endpoint bound for the core circle.

## 4. Fixed and growing dimension, including aggregate heights

For `m>=2`, `k=m+1`, and `M>=16k^2`, the inequalities `L>=M/8`
give

```text
eta<=16b/sqrt(M),            zeta<=4sqrt(2)k/sqrt(M).
```

Equations (2)--(4) imply the convenient explicit bounds

```text
epsilon<=8kb/sqrt(M),
log|K_i|/w <=kb/sqrt(M),
log max(1,|Y_i|)/w <=sqrt(2)kb/sqrt(M),
log|t_ij|/w <=(2+sqrt(2))kb/sqrt(M)<4kb/sqrt(M).           (7)
```

Thus every fixed `m` gives the claimed genuine profile with
`log Norm(H_T)=w+o(w)` and all individual correction and residue
heights `o(w)`. More generally the same conclusion holds whenever

```text
(m+1)2^m/sqrt(M) ->0.                                    (8)
```

This includes `m=floor(c log_2 M)-1` for every fixed `c<1/2`.
A closer approach to the boundary permitted by (8) is
`m<=0.5 log_2 M-log_2 log_2 M-h(M)` with `h(M)->infinity`.
These are sufficient ranges from the estimates, not impossibility
claims beyond them.

For the products used in subsequent arithmetic theorems, put

```text
K_all=product_i Norm(K_i),
Y_all=product_i |Y_i|,
B_all=product_(i<j)|t_ij|.
```

Then, with no omitted dependence on the number of rows,

```text
log K_all <=2m s,
log Y_all <=m zeta W/4,
log B_all <=binom(m,2)(zeta W/4+2s).                     (9)
```

Consequently their aggregate heights divided by `w` are respectively
`O(k^2 2^m/sqrt(M))`, `O(k^2 2^m/sqrt(M))`, and
`O(k^3 2^m/sqrt(M))`. All three are subpower if

```text
(m+1)^3 2^m/sqrt(M) ->0,                                (10)
```

which again holds throughout every fixed `c<1/2`. Near the boundary,
one sufficient range is
`m<=0.5 log_2 M-3log_2 log_2 M-h(M)` with `h(M)->infinity`.
An arithmetic theorem whose constants themselves depend on `m` still
needs those constants controlled before it can be applied uniformly
in growing dimension. For a fixed `m`, no such issue arises.

Finally, (1) gives

```text
w >=4(M-1)lambda/2^m ->infinity                          (11)
```

in the fixed-dimensional case and throughout (8). The intrinsic
pre-empty-removal scale is `W_core/2^m=(1-2s/W)w`; it is
`(1+o(1))w` in these ranges, so either convention gives the same
asymptotic profile.

## 5. A fixed positive exponent obstruction suffices

The following is a sufficient new arithmetic theorem; it remains
unproved:

> There exist a fixed `m>=2`, a fixed `delta>0`, and a threshold
> `w_0`, such that for no `w>=w_0` do there exist conjugate-primitive
> Gaussian rows `P_i=X_i+iY_i`, oriented blocks `H_T` for every
> nonempty `T subset [m]`, and Gaussian correcting factors `K_i`,
> with disjoint rational-prime block supports and all of
>
> ```text
> P_i=K_i product_(T containing i)H_T,
> |log Norm(H_T)-w|<=delta w,
> 0<|Y_i|,             Delta_ij!=0,
> Delta_ij=t_ij product_(T containing i,j)Norm(H_T),
> t_ij in Z\{0},
> max_i log max(1,|K_i|,|Y_i|)<=delta w,
> max_(i<j) log|t_ij|<=delta w.
> ```

The formulation is existential in `m`; it cannot hold at every row
count. An exact [two-row family](small_imaginary_row_lattice_transference.md)
has the full three-block profile with `Y_1=Y_2=1`. The unbounded
[three-row polynomial family](balanced_three_row_gaussian_polynomial_family.md)
has all seven blocks and fixed `Y=(16,90,14)`. Both have primitive
rows, disjoint block supports, and bounded residuals. Hence a universal
exclusion theorem of precisely the displayed form must use `m>=4`.
The three-row counterexample permits nonconstant imaginary coordinates,
as the endpoint extraction requires; it does not refute a theorem
restricted to `Y_i=1`.

The theorem may use the exact simultaneous norm, determinant, and
orientation equations; they all hold in the extracted data. The
correcting factors must be allowed to share primes with the blocks.
Nonconstant `Y_i`, arbitrary split-prime support, and arbitrary prime
powers must also be allowed unless a further reduction is proved.

It is enough to prove this for one fixed positive `delta`, however
small. Shrink it to at most `1/2`. With `k=m+1`, any endpoint cluster
of cardinality at least

```text
M_0=max(16k^2, (8k2^m/delta)^2,
        1+2^m max(w_0,0)/(4lambda))                      (12)
```

contradicts that theorem by (7) and (11). Integer ceilings and the
fixed arc-subdivision factor turn (12) into a bound independent of
the radius. If the proposed theorem instead bounds aggregate heights
by `delta w`, (9) gives the same implication with a larger fixed
constant in (12).

No convergence rate of `M` relative to the radius is needed. The
extraction works even if the original unbounded counts grow more
slowly than every suggested unbounded function of the radius. The
strict-gap estimate (1) and the fixed tuple selection are what make
this quantifier order valid.

The [primewise checker](check_endpoint_full_profile_quantifiers.py)
verifies the exact gcd-depth identity and its central-truncation loss
for exhaustive bounded allocation tuples, including ties. The
Archimedean and asymptotic claims above follow from the displayed
inequalities and the audited extraction, not from these finite checks.
