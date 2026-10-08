# A sharp cap on summed coalition orders of binary invariants

Let `F` be a nonzero `SL_2`-invariant polynomial in `m>=2`
homogeneous binary rows, of degree `delta` in each row. For a subset
`T subset [m]`, let `nu_T=ord_T F` be its generic order when the
directions indexed by `T` coalesce. Set `nu_empty=nu_{i}=0`, as the
empty and singleton conditions impose no collision. The full
coalition is included, using an affine collision family after a
common projective coordinate change.

**Theorem.** With `E=m delta/2`,

```text
sum_(T subset [m]) nu_T <= 2^(m-2) E
                         = delta m 2^(m-3).                 (1)
```

The bound is sharp: every nonzero bracket monomial of row degree
`delta` attains equality. It applies to an arbitrary linear
combination of bracket monomials, including cancellations.

## Affine support and collision orders

Put `f(z)=F((1,z_i)_(i=1)^m)`. Invariance under binary translations
gives `f(z_1+c,...,z_m+c)=f(z)`. The diagonal torus gives
homogeneity of total degree `E`; in particular `m delta` is even
when `F!=0`. Every individual degree is at most `delta`.

Write `A=Supp(f) subset {0,...,delta}^m`. Every `alpha in A`
has `|alpha|=E`. Moving a generic collision center to zero by the
translation action and substituting `z_i=t u_i` for `i in T`
shows exactly that

```text
nu_T=min_(alpha in A) alpha(T).                              (2)
```

The coefficient of the first nonzero power of `t` is a nonzero
polynomial in the generic `u_i` and outside directions, so there is
no generic cancellation in (2). This also covers `T=[m]`, where
`nu_[m]=E`. A singleton has order zero: if every term of `f` were
divisible by `z_i`, `F` would vanish with row `i=(1,0)` and then,
by `SL_2` transitivity, identically.

Define the centered support function

```text
h(S)=max_(alpha in A)(alpha(S)-delta |S|/2).                 (3)
```

Since `|alpha|=E`, equation (2) is equivalently

```text
nu_T=E-max_(alpha in A) alpha(T^c)
    =delta |T|/2-h(T^c).                                    (4)
```

## The coefficient-support grid bound

The [polynomial grid lemma and invariant coefficient construction](sl2_invariant_average_cut_width.md)
prove, for a nonzero
translation-invariant homogeneous polynomial of total degree `E`
and individual degree bounds `delta`, that

```text
average_(S subset [m]) h(S) >= E/4.                         (5)
```

For completeness, polarize each variable into `delta` Boolean
variables. Translation invariance and homogeneity make the
polarized polynomial vanish at Boolean vertices of weight other
than `E`. Its separately symmetric grid-value polynomial has
total degree at most `E` and nonzero grid locus exactly `A`.
The grid lemma bounds the average distance from a random corner
`B(S)` to `A` by `E/2`; that distance is `E-2h(S)`, giving (5).
The [multiplicity extension](polarized_coefficient_support_multiplicity.md)
gives the same bound as its special case with multiplicity `E` at
the common diagonal.

Averaging (4) over all subsets gives

```text
average_T nu_T=E/2-average_S h(S)<=E/4.
```

Multiplication by `2^m` proves (1).

The binary coordinate swap sends the support `A` to
`delta*1-A`. It follows that

```text
nu_T-nu_(T^c)=delta (|T|-m/2).                              (6)
```

Thus `nu_[m]=E` and `nu_([m] minus {i})=E-delta`; these are the
full- and near-full-coalition conventions used in the
[fusion-content note](terminal_kernel_fusion_content.md). They are
already counted in (1), with no boundary-divisor exception.

For equality, take `F` to be a product of `E` binary brackets
`[ij]=P_iY_j-P_jY_i`, with each label occurring `delta` times.
For every edge `ij`, its bracket has order one precisely at the
`2^(m-2)` coalitions containing both endpoints. Summing over the
`E` edges gives equality in (1).

## Consequence for kernel Pluecker coordinates

The primitive polynomial Pluecker tuple `W` in the
[determinant-line note](terminal_kernel_determinant_line_height.md)
has row degree `delta`, and every coordinate is an `SL_2`
invariant. Let `nu_T(W)` be the minimum generic `T`-order among its
coordinates. For every nonzero coordinate `W_S`, (1) gives

```text
sum_T (ord_T W_S-nu_T(W))
 <= delta m 2^(m-3)-sum_T nu_T(W).                            (7)
```

In the full Boolean profile the right side is exactly the generic
height numerator `B(n,d)` (called `C_k` in the terminal
fusion-content calculation). Hence
generic coalition orders of one exterior coordinate cannot yield
a strictly larger coefficient floor. Equation (7) does not bound
extra valuations caused by arithmetic specialization of the
residual directions, and it does not prove a radius-uniform result.
