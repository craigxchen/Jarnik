# Map geometry of a full Boolean polynomial profile

This note records exact geometric identities for a hypothetical rational
polynomial full-Boolean profile. They apply for every `m>=2` and every
factor degree `e>=1`. None excludes a profile in four or more rows.
The arithmetic significance of finding such profiles is explained in
[the specialization note](boolean_polynomial_integer_specialization.md).
The corresponding system over real algebraic coefficients now has a
[certified full four-row solution](four_row_real_polynomial_profile_certificate.md).
Its factors are not known to be rational over `Q(i)`. Thus a four-row
obstruction must retain arithmetic information beyond the real identities
and positivity discussed here.

## Setup and pair differences

For every nonempty `T subset [m]`, let `H_T in Q(i)[t]` be monic of
degree `e`. Assume the `2(2^m-1)` polynomials `H_T, bar(H_T)` are
pairwise coprime. Suppose

```text
P_i=K_i product_(T contains i) H_T = X_i+iY_i,
X_i in Q[t],  Y_i in Q\{0},  K_i in Q\{0}.
```

The leading coefficient `K_i` is real because the imaginary part of
`P_i` is constant. Put

```text
d=2^(m-1)e,       F_i=X_i/Y_i,       a_i=K_i/Y_i,
A_i=product_(T contains i) H_T,
U_i=F_i+i=a_i A_i,     V_i=F_i-i=a_i bar(A_i).
```

Thus `F_i` is a real polynomial of degree `d`, with leading coefficient
`a_i`. For `i!=j`, set

```text
G_ij=product_(T contains i,j) H_T bar(H_T).
```

It has degree `2*2^(m-2)e=d` and divides `F_i-F_j`: at each root of
`H_T` with `T` containing both labels, both `F_i,F_j` equal `-i`,
and at each conjugate root both equal `i`, with multiplicities. The
difference cannot vanish identically, because `U_i,U_j` have distinct
incident factor sets. Consequently

```text
F_i-F_j=lambda_ij G_ij,        lambda_ij=a_i-a_j != 0.       (1)
```

In particular the leading coefficients `a_i` are pairwise distinct.

## Linearly nondegenerate map and its phase functions

The morphism `f:P^1 -> P^m` defined on the affine line by
`f(t)=[1:F_1(t):...:F_m(t)]` is linearly nondegenerate. Indeed, if
`b+sum_i c_i F_i=0`, evaluation at a root of `H_[m]` gives
`b-i sum_i c_i=0`; evaluation at a root of `bar(H_[m])` gives
`b+i sum_i c_i=0`. Hence `b=sum_i c_i=0`. For any `j`, at a root
of `H_([m]\{j})`, all `F_i=-i` except `F_j`, while `F_j+i` is
nonzero by coprimality. The relation reduces to `c_j(F_j+i)=0`,
so every `c_j=0`. This argument does not require the factors to be
squarefree.

Let `L_0` be the hyperplane `z_0=0` and let `L_i^+` and `L_i^-`
be `z_i+i z_0=0` and `z_i-i z_0=0`. The pullback of `L_0` is
`d[Infinity]`. The pullback of `L_i^+` is the zero divisor of
`U_i`, namely `sum_(T contains i) div(H_T)`; for `L_i^-` it is the
conjugate divisor. These are exactly the Boolean zero/pole cuts.

Define `B_i=U_i/V_i`. Since `U_i,V_i` are coprime, `B_i` has degree
`d`, and

```text
B_i-1 = 2i/V_i.
```

Thus `B_i` has contact order `d` with the value `1` at infinity.
For a pair, cancel the shared norm product `G_ij` explicitly:

```text
B_i/B_j = N_ij/D_ij,
N_ij = U_i V_j/G_ij
     = a_i a_j product_(i in T,j notin T) H_T
               product_(j in T,i notin T) bar(H_T),
D_ij = V_i U_j/G_ij
     = a_i a_j product_(i in T,j notin T) bar(H_T)
               product_(j in T,i notin T) H_T.
```

There are `2^(m-1)` separating blocks, so both reduced polynomials
have degree `d`; they are coprime by hypothesis. Equation (1) gives

```text
N_ij-D_ij
 = ((F_i+i)(F_j-i)-(F_i-i)(F_j+i))/G_ij
 = 2i(F_j-F_i)/G_ij = -2i lambda_ij.                 (2)
```

Therefore `B_i/B_j` also has degree `d`, tends to `1` at infinity,
and has contact order `d` there. In particular the pair contact is
exactly the degree, after accounting for every shared denominator
factor; using the uncancelled degree `2d` would overcount it.

By Riemann--Hurwitz, either degree-`d` map `B_i` or `B_i/B_j` has
total ramification `2d-2`, of which infinity uses `d-1`. The
remaining finite budget is `d-1`. If every `H_T` is squarefree,
all finite zeros and poles of these phase maps are simple, so this
budget lies over values other than `0` and `Infinity`. Squarefreeness
is an additional hypothesis here; a constant imaginary row does
not itself imply it. The budget and Boolean cuts alone yield no
contradiction.

## Pair power sums add no linear moment equations

For `1<=k<=d-1`, write

```text
S_T(k)=Im(sum_(alpha root of H_T, with multiplicity) alpha^k).
```

The nonleading coefficients of `U_i=F_i+i` in degrees `1,...,d-1`
are real. Newton's identities say equivalently that

```text
sum_(T contains i) S_T(k)=0             (1<=k<=d-1).       (3)
```

The reduced pair numerator `N_ij` has roots from `H_T` when
`i in T,j notin T`, and from `bar(H_T)` when `j in T,i notin T`.
Its `k`th imaginary root-power sum is therefore

```text
sum_T (1_(i in T)-1_(j in T)) S_T(k),
```

the difference of the `i` and `j` instances of (3). Thus the
constant-imaginary property of `N_ij` in (2) provides no additional
*linear* condition on the block power sums. This does not address
nonlinear compatibility among the roots.

## Rank of the target hyperplane arrangement

The normals of `L_0,L_i^+,L_i^-` are `e_0,e_i+i e_0,e_i-i e_0`.
For any subcollection, let `J` be its set of used indices `i`, and
write `k=|J|`. If no index appears with both signs, its rank is
`k+1` when `L_0` is present and `k` otherwise. If some index
appears with both signs, those two normals span `e_0,e_i`;
every other used index adds one dimension and `L_0` adds none.
The rank is then `k+1`.

The largest subcollection with nonempty projective intersection
(rank at most `m`) has `2m-1` members: `L_0` and both signs on
`m-1` indices achieve it. Every `2m` hyperplanes have rank `m+1`.
Thus the exact subgeneral-position parameter of this arrangement
is `N=2m-1`, for `q=2m+1` hyperplanes. The standard Cartan--Nochka
theorem quoted as Theorem 1.2 in [Dethloff--Tan--Thai,
*An Extension of the Cartan--Nochka Second Main Theorem for
Hypersurfaces*](https://arxiv.org/pdf/0911.2562) assumes
`q>=2N-m+1=3m-1`. Our arrangement misses that hypothesis for
`m>2`; for `m=2` it meets the threshold but its left coefficient
`q-2N+m-1` is zero. Hence that theorem, applied directly to this
full arrangement, yields no exclusion. The rank computation itself
is elementary and does not depend on the theorem.

For `m>=3`, deleting hyperplanes does not produce an ordinary
general-position subarrangement of size `m+2`. If such a selection
omits `L_0`, its `m+2` signed members include both signs for at least
two indices. Choose those four and any `m-3` other selected members:
these `m+1` normals use at most `m-1` indices and have rank at most
`m`. If the selection includes `L_0`, its other `m+1` members
include both signs for some index. Choose `L_0`, that pair, and any
`m-2` other selected members; again `m+1` normals have rank at most
`m`. Thus no size-`m+2` subcollection is in ordinary general
position. This observation concerns only this fixed hyperplane
arrangement, not every possible target-space argument.
