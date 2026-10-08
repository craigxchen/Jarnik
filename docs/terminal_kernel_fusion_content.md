# Fusion residues give the terminal kernel's coalition content

The terminal coherent-kernel Pluecker vector has a computable exact
order at every generic coalition collision. This is a statement about
the actual polynomial numerator, not just the class of its determinant
line bundle. It supplies additional universal core divisibility beyond
the majority bound. The resulting height estimate is still above the
original single-pair threshold in the checked cases `k=1,...,15`; no uniform
endpoint theorem follows here.

## 1. The bundle, normalization, and residue data

Put `m=6k=3d`, `d=2k`, and use the level `ell=d-1` for `sl_3`, with
all `m` labels equal to the defining fundamental weight `omega_1`.
Let `A` be the constant coinvariant space, and let `V` be the
conformal-block **quotient** bundle. The
[exact highest-root identity](terminal_coherent_conformal_block_identity.md)
identifies the coherent source kernel with `V^* subset A^*` at every
distinct configuration. Write

```text
h=rank V=r_k-rho_k,
delta=h(6k-1,2k),       q=m/2.
```

Use the primitive integral polynomial Pluecker tuple `W` constructed
in [the determinant-line note](terminal_kernel_determinant_line_height.md).
It has degree `delta` in every homogeneous binary row. On the affine
chart `(1,z_i)`, every coordinate is homogeneous of total degree
`m*delta/2`. The tuple is nowhere zero on the distinct-point locus.
This last fact uses the constant-rank determinant bundle and the
principal-open/Picard-group argument in that note; polynomial gcd
one alone would not establish it.

The source facts used below are [Fakhruddin, Lemma 2.5 and Section
3.1.4, equation (3.4)](https://arxiv.org/pdf/0904.2918): the constant
coinvariant quotient is surjective also for stable genus-zero curves,
and the quotient KZ residue is scalar on each fusion summand.
The [affine KZ description and its compatibility with conformal-block
quotients](https://arxiv.org/pdf/2302.00798) are in Section 4 of
Belkale–Fakhruddin. We use the **dual** of Fakhruddin's quotient
connection throughout, so determinant residues change sign.

For a coalition of `r` labels, define the quotient determinant residue

```text
tau_r = [sum_mu c(mu) f_r(mu) f_(m-r)(mu^*)
         - r*h*c(omega_1)] / [2*(ell+3)],              (1)
```

where `f_r(mu)` is the fusion multiplicity of `mu` in `r` copies of
`omega_1`. For `sl_3`, in Dynkin coordinates `mu=(a,b)`,

```text
mu^*=(b,a),
c(a,b)=2(a^2+ab+b^2)/3+2(a+b),
c(omega_1)=8/3,       ell+3=d+2.
```

Factorization makes `f_r(mu) f_(m-r)(mu^*)` the rank of the
corresponding summand. Thus the determinant residue of `V^*` is
`-tau_r`, in the affine coalition gauge of (1).

For the direct use of Fakhruddin's type-(4) chart below we restrict to
`2<=r<=q`. There are then at least three labels outside the coalition,
since `m>=6`. Normalize three outside labels to `0,1,infinity`.
The normalization and its coordinate changes are regular units along
the collision, so they do not change this residue. Larger coalitions
will be handled by binary invariant symmetry, without extending that
chart beyond its stated range.

## 2. Exact scalar connection on the polynomial determinant section

On the distinct affine configuration space, `W` spans `det(V^*)`.
The induced determinant connection therefore has the form

```text
nabla W = alpha * W                                  (2)
```

for a rational one-form `alpha` in the coordinates `z_i`.
We claim

```text
alpha = -tau_2 * sum_(i<j) dlog(z_i-z_j).              (3)
```

Here is the argument fixing both the normalization and the sign.
The ambient KZ connection on the constant space has only logarithmic
poles at pair diagonals. The polynomial tuple `W` is primitive, so
at the generic point of each pair diagonal at least one coordinate
is a unit. Computing (2) with that coordinate shows that `alpha`
has at most a logarithmic pole there. It has no other poles: on
the distinct locus `W` is a nonvanishing frame, and any possible
remaining pole would have a codimension-one component.

The quotient `A -> V` extends surjectively over the stable pair
collision. Its dual and determinant are subbundles of the constant
ambient spaces. The nonzero limit of the primitive tuple `W` is
therefore a local frame of this determinant subbundle. Its residue
in (2) is `-tau_2`. All labels have the same weight, so this residue
is the same at every pair diagonal.

Subtract the right side of (3) from `alpha`. The difference has no
poles in codimension one, hence is a polynomial one-form on affine
space. Both `W` and the KZ connection are compatible with simultaneous
dilation `z_i |-> c*z_i`: `W` is multiplied by the constant
`c^(m*delta/2)`, and the KZ logarithmic forms are unchanged.
Consequently the difference is dilation-invariant. Every nonzero
polynomial one-form has strictly positive dilation weight, so the
difference is zero. This proves (3), including the absence of a
hidden regular one-form.

## 3. Exact coalition order

Fix a coalition `T` of size `r<=q` and a generic collision family

```text
z_i=s+t*u_i  (i in T),       z_j fixed (j outside T),
```

where the `u_i` are distinct and the outside points avoid `s`.
Let `nu_r` be the minimum `t`-order of the coordinates of `W` in
this family. Along the stable family, choose a nonvanishing local
frame `F` of `det(V^*)` inside the constant ambient exterior space.
Because `A -> V` remains a quotient of vector bundles, at least
one coordinate of `F` is a unit. Absorbing a scalar unit into `F`,
we have

```text
W=t^(nu_r) F.
```

The determinant connection residue in this frame is `-tau_r`.
Thus the residue of (2) is `nu_r-tau_r`. On the other hand, exactly
`binom(r,2)` pair differences in (3) vanish to first order along
the family. Equating the two residues proves

```text
nu_r=tau_r-binom(r,2)*tau_2,       2<=r<=q.             (4)
```

In particular these rational fusion expressions are nonnegative
integers, because they are polynomial vanishing orders. We also
have `nu_0=nu_1=nu_2=0`.

For completeness, the orders on complementary coalitions satisfy

```text
nu_r-nu_(m-r)=delta*(r-q).                             (5)
```

This follows directly from invariant monomial support. Total degree
in the second binary coordinates is `m*delta/2`. Swapping the two
binary coordinates simultaneously preserves the support, up to
one constant sign, and replaces the degree on `T` by `delta*r`
minus that degree. Thus the minimum degree on `T` is
`delta*r-m*delta/2` plus the minimum on its complement. Generic
collision order is exactly this minimum degree, proving (5).
Together (4) and (5) determine every coalition order, including
`r=m-2,m-1,m`, without an unstable use of the type-(4) chart.

Equation (4) is an equality of numerator orders, obtained by comparing
connections on actual embedded determinant frames. No inference from
an equality of divisor classes is used.

## 4. Universal arithmetic content

The generic collision order has a coefficientwise consequence.
For every subset `T`, every monomial in every coordinate of `W`
has total first-coordinate degree on `T` at least `nu_(|T|)`:

```text
sum_(i in T) a_i >= nu_(|T|)
for each monomial product_i P_i^(a_i) Y_i^(delta-a_i).  (6)
```

Indeed, substitute `P_i=t*A_i` for `i in T`, keeping all other
row coordinates generic. A coefficient of a smaller power of `t`
would be a nonzero polynomial in the generic `A_i,Y_i` and outside
rows, contradicting the exact generic order. Such a coefficient
cannot vanish only generically without vanishing identically.
Binary invariance identifies the common collapsed direction with
the affine collision used in Section 3.

Use the primitive odd split core factorization from
[the higher-kernel content note](higher_kernel_evaluation_heights.md):

```text
P_i=X_i+iY_i=K_i product_(T containing i) H_T,
n_T=Norm(H_T).
```

The oriented core blocks and their conjugates have disjoint support.
The determinant-one shear gives `W(P,Y)=W(X,Y)`, an integer tuple.
Applying (6) to each block makes
`product_T H_T^(nu_(|T|))` divide each coordinate termwise.
Conjugating the ordinary integer values then gives the universal
integer divisor

```text
G_fusion=product_T n_T^(nu_(|T|))
       divides gcd of the coordinates of W(X,Y).     (7)
```

This lower bound holds for every actual residue tuple, including
special tuples where still more cancellation occurs. It is not an
exact upper bound for the specialized integer gcd. The monomial
argument also adds over nested allocation layers with the same
oriented prime, when such layers are retained before truncation.

The majority bound is recovered from (5) and `nu_r>=0`:
`nu_r>=delta*max(0,r-q)`. The excess is symmetric in `r,m-r`.
Its first nonzero values are:

| `k` | excess over `delta*max(0,r-q)` |
|---:|:---|
| 1 | none |
| 2 | `r=6: 15` |
| 3 | `r=8,9,10: 210,1008,210` |
| 4 | `r=10,11,12,13,14: 3003,17589,54450,17589,3003` |

## 5. Refined height bound and finite comparison

In the full core profile, with `log n_T=w+o(w)` and correction
height `o(w)`, the raw bound for the determinant numerator remains
`delta*m*2^(m-3)w+o(w)`. Dividing by (7) yields

```text
log covol(ker_Z E_P) <= C_k*w+o(w),
C_k=delta*m*2^(m-3)-sum_(r=0)^m binom(m,r)*nu_r.        (8)
```

The first minimum consequently has upper exponent `C_k/h`.
Exact arithmetic gives:

| `k` | kernel rank `h` | `C_k` | `C_k/h` | single-pair threshold |
|---:|---:|---:|---:|---:|
| 1 | 1 | 18 | 18 | 8 |
| 2 | 341 | 171,600 | 503.225806 | 32 |
| 3 | 83,028 | 815,673,600 | 9,824.078624 | 128 |
| 4 | 23,193,775 | 3,196,140,304,896 | 137,801.643109 | 512 |
| 5 | 7,638,717,450 | 11,469,050,955,648,000 | 1,501,436.730802 | 2,048 |

The checker verifies that this estimate remains above threshold
through `k=15`. No all-`k` comparison for this refined estimate is
claimed. Additional arithmetic content at special actual configurations
or a different use of the determinant line remains possible.

The later [higher-rank calculation](higher_rank_terminal_fusion_grid.md)
combines all coordinate pairs to strengthen the arithmetic threshold
to `U(3,2k)`. In particular `U(3,2)=18`, so the six-row comparison
becomes equality; the table above records the original single-pair
comparison only.

For the finite fusion computation, the fundamental minuscule Pieri
rule is the alcove-truncated sum of fundamental weight steps; see
[Andersen–Stroppel, equation (2.2) and Section 2.2](https://people.math.uni-bonn.de/stroppel/Fusion.pdf).
In `sl_3` Dynkin coordinates it gives the transitions

```text
(a,b) -> (a+1,b), (a-1,b+1), (a,b-1),
```

retaining precisely the states with `a,b>=0` and `a+b<=ell`.
Starting at `(0,0)`, this computes the multiplicities in (1).
The exact implementation is
[check_terminal_fusion_collision_candidate.py](check_terminal_fusion_collision_candidate.py).
The historical filename is retained; the connection proof in
Sections 2–3 supplies the source bridge independently of that checker.

The independent [six-row connection-sign checker](check_terminal_kz_sign.py)
works directly in ternary polynomial coefficients. With
`Omega_ij=(row transposition ij)-1/3`, `kappa=4`, and the explicit
kernel polynomial `Q_z=sum_j M_j(z)G_j`, it verifies
`(partial_i-sum_(j!=i)Omega_ij/(4(z_i-z_j)))Q_z
 = (sum_(j!=i)1/(3(z_i-z_j)))Q_z`
at all six directions in three distinct rational configurations.
The opposite connection sign fails. This calibrates the dual convention
against the actual polynomial model, rather than treating different
authors' horizontal-equation conventions as interchangeable.
