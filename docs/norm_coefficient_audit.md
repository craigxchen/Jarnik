# Exact norm coefficients, coherent conductor components, and auxiliary forms

This note proves two precise statements about the central system in
[endpoint_central_truncation.md](endpoint_central_truncation.md). First,
the norm and pair-determinant equations generate the full coherent local
conductor condition, including every power of its defining ideal. Second,
conductor-sized **row norms** can be eliminated from auxiliary forms
without discarding their divisibility; the degree-versus-divisibility
comparison also survives generic common norm denominators. This extends
the moving-coefficient comparison to a useful class of norm resultants
and primitive Gram expressions. It does not eliminate the independent
cut norms, and does not prove the uniform endpoint bound.

## 1. Arithmetic setting

Use one common-unit tuple from the central truncation. With anchor zero,
write its primitive Gaussian numerators as

```text
h_i=x_i+i t_i=K_i A_i,
A_i=product_(S contains i) H_S,
q_i=N(h_i)=x_i^2+t_i^2,
n_S=N(H_S),
q_i=N(K_i) product_(S contains i) n_S.
```

Every `t_i` is a nonzero rational integer. The `H_S`, including the
empty block, are independent Gaussian blocks and are coprime to all
their conjugates. Put `m=k-1`, `B=2^m`, `W=log R^2`, and
`W_core=sum_S log n_S`. The proved bounds are

```text
log n_S=(1+O(eta_core))W_core/B,
W_core=W-O(kW/sqrt(M)),
log|K_i|<=kW/sqrt(M),
log max(1,|t_i|)=O(kW/sqrt(M)),
log|h_i|=W/4+O(kW/sqrt(M)).
```

Here `k=floor(c log_2 M)` for a fixed `c<1/2`, and
`eta_core=O(k 2^k/sqrt(M))=o(1)`. In particular
`log q_i=W/2+o(W)`: the row norms have conductor-sized height and
cannot simply be declared small polynomial coefficients.

For any two rows set

```text
g_ij=gcd_G(h_i,h_j),
G_ij=N(g_ij),
D_ij=x_i t_j-x_j t_i.
```

The exact primitive reanchoring identity gives
`D_ij=G_ij t_ij` up to a consistent sign convention, where `t_ij`
is the primitive pair residue. In particular the core norm product
`product_(S contains i,j) n_S` divides `D_ij` as a rational integer.

## 2. The norm and pair minors generate the coherent local condition

Let `J` be a nonempty set of row indices and choose `s in J`. Work in
the polynomial ring

```text
A=Z[i, 1/(2 product_(i in J) t_i)] [x_i : i in J].
```

The `t_i` can either be fixed nonzero integers or invertible
indeterminates. Define the two coherent isotropic ideals

```text
I_+=(x_i+i t_i : i in J),
I_-=(x_i-i t_i : i in J),
J_coh=(x_s^2+t_s^2,
       t_s x_i-t_i x_s : i in J\{s}).
```

Then the exact ideal identities are

```text
J_coh=I_+ intersect I_-=I_+ I_-,
J_coh^r=I_+^r intersect I_-^r   for every r>=1.       (1)
```

To prove the first equality, change coordinates to
`x_s` and `u_i=x_i-(t_i/t_s)x_s`. This is invertible over `A`.
The ideals become

```text
I_+=(u_i, x_s+i t_s),
I_-=(u_i, x_s-i t_s).
```

Their sum is the whole ring, since their difference contains the unit
`2i t_s`. Modulo the `u_i`, their intersection is therefore generated
by `(x_s+i t_s)(x_s-i t_s)=x_s^2+t_s^2`. This proves the first
identity. Powers of two comaximal ideals are comaximal and their
intersection is their product, proving the second.

In particular every other row norm is already in `J_coh`, because

```text
t_s^2 (x_i^2+t_i^2)-t_i^2 (x_s^2+t_s^2)
  =(t_s x_i-t_i x_s)(t_s x_i+t_i x_s).               (2)
```

Norm vanishing without the pair minors would allow all `2^|J|`
independent choices of isotropic signs. The minors restrict this to
the two coherent choices. Conversely (1) says that, after this
restriction, the norm equations do not conceal an additional local
component or an additional generic order of divisibility.

### No high-prime-power exception is needed

Suppose the split prime `p=pi bar(pi)` occurs in a core block `H_J`.
For each `i in J`, `pi` divides `h_i`, while `bar(pi)` does not,
because `h_i` is conjugate-coprime. It follows that **`p` does not
divide `t_i`**. Otherwise `h_i-bar(h_i)=2i t_i` would imply that
`pi` also divides `bar(h_i)`, a contradiction. The prime is odd,
so all the inversions in (1) are valid at this conductor prime.

More concretely, for every exponent `e>=1`, the rational congruences

```text
x_s^2+t_s^2=0 mod p^e,
t_s x_i-t_i x_s=0 mod p^e   for i in J\{s}          (3)
```

are equivalent to a common isotropic sign modulo `p^e`: choose either
root `r` of `r^2=-1 mod p^e`, and all coordinates satisfy
`x_i=r t_i mod p^e`. The anchor norm has precisely the two choices
because `t_s` and `2` are units; the minors propagate that choice.
Thus this statement keeps the entire conductor exponent. It does
not discard a large prime power merely because a coefficient might
have shared its prime support.

An exact finite enumeration checked (3) on 136,050 residue tuples for
split primes `5,13,17`, prime-power depths one and two, and sets of
two or three rows, omitting only the largest enumeration boxes. No
exception occurred. The ideal proof supplies the unrestricted result.

## 3. Exact elimination of conductor-sized row norm coefficients

Let `F(X,Q)` be a polynomial with Gaussian rational coefficients.
The letters `Q_i` are row-norm variables, not independent cut norms.
Let

```text
F_red(X)=F(X, X_1^2+t_1^2,...,X_m^2+t_m^2).        (4)
```

The evaluation identity is exact:

```text
F(x,q)=F_red(x).
```

Thus there are two possibilities.

* If `F_red=0` as a polynomial, the proposed auxiliary value is a
  universal norm identity. It cannot be a nonzero integer in a
  product-formula contradiction.
* Otherwise let `d=deg F_red`. Its coefficients have small height
  whenever the original coefficients and the `t_i` do and the
  original weighted degree is bounded, assigning weight one to
  `X_i` and weight two to `Q_i`. Exact leading cancellations are
  incorporated by using `d`, rather than the unreduced weighted
  degree.

This substitution also preserves prescribed norm divisibility.
At a prime of `H_S`, each `Q_i` with `i in S` is replaced by

```text
(X_i+i t_i)(X_i-i t_i).
```

The first factor has the prescribed oriented block divisibility,
and the second is a unit at this prime. The roles reverse at the
conjugate prime. Consequently a row norm contributes exactly one
copy of the rational block norm, as its actual factorization requires.
No conductor-sized coefficient has been replaced by a coefficient
whose prime factors are forgotten: it has been replaced by its exact
polynomial factorization in the large coordinates.

For `F_red!=0`, let `nu_S^+` and `nu_S^-` be its generic orders
along `(X_i+i t_i:i in S)` and `(X_i-i t_i:i in S)`. The moving
two-point-grid argument gives

```text
sum_S (nu_S^+ + nu_S^-) <= d B/2.                 (5)
```

Therefore its prescribed core divisor has logarithmic modulus at most

```text
d(1+eta_core)W_core/4,                            (6)
```

before small coefficient-prime corrections. Its ordinary size bound
is `exp(dW/4+o(W))`. Products of row norms attain the same leading
exponent on both sides: for example `F_red=X_i^2+t_i^2` has `d=2`,
and the actual leading logarithmic norm divisor and size are both
`W/2`. This gives a sharp obstruction for this enlarged class of
auxiliaries, rather than an argument that norm factors are absent.

### Growing dimension does not invalidate the bounded-degree statement

The comparison remains uniform for the central choice of growing `k`
when the degree cap `d<=D` is fixed and coefficient logarithmic heights
are `O_D(k^a W/sqrt(M))` for any fixed power `a`. There are at most
`binom(m+D,D)` monomials. Translating a coefficient expansion to any
one of the moving coordinate subspaces costs at most a fixed-degree
multiple of the residue heights plus the logarithm of this monomial
count. Choose one nonzero minimal-order translated coefficient for
each of the `2B` subspaces. Their combined prime-divisibility cost is

```text
O_D(2^k k^a W/sqrt(M))=o(W),                      (7)
```

with a harmless increase of `a` if rational denominators are cleared.
This follows for every fixed `c<1/2`. Profile errors in (6) are also
`o(W)`. Thus allowing a bounded-degree norm auxiliary to vary with
the extracted growing tuple does not introduce a hidden leading
coefficient-prime gain.

## 4. Common norm denominators cancel the same budget

There is a further exact extension. Let

```text
D_norm=product_S n_S^(a_S),        a_S>=0,
E=F_red(x)/D_norm.
```

Suppose the denominator is justified by the generic block orders:
`a_S<=nu_S^+` and `a_S<=nu_S^-` for every `S`. Then `E` is
integral after the same small coefficient clearing as above. The
remaining prescribed divisor has logarithmic modulus at most

```text
dW/4-log D_norm+o(W),                            (8)
```

because subtracting `a_S` at both conjugate places subtracts exactly
`a_S log n_S`. The ordinary quotient size bound subtracts the same
quantity:

```text
log|E|<=dW/4-log D_norm+o(W).                     (9)
```

Thus passing to primitive numerators through a common norm denominator
does not itself improve the exponent comparison. Examples include
`h_i bar(h_j)/N(gcd_G(h_i,h_j))`, its imaginary part, and products
of such expressions that share one generic common denominator. The
small correcting parts of the gcd have only negligible height and
can be tracked in the coefficient error.

The generic denominator hypothesis is substantive. If a numerator
is divisible by the proposed denominator only because of a special
cancellation at the actual outside coordinates, that is extra
arithmetic information, not a consequence of (5). Likewise sums with
different norm denominators may introduce independent cut-norm
coefficients upon clearing; they are not automatically covered by (8).

## 5. What remains beyond this comparison

The individually independent cut norms `n_S` have not been eliminated.
There are `B` of them and only `m` equations
`q_i=N(K_i) product_(S contains i) n_S`. Allowing a polynomial in
the `n_S` therefore goes beyond merely allowing row norms as
coefficients. Its relations can connect the actual unit values left
after division by one core block. Neither the moving-grid inequality
nor the local ideal identity (1) bounds all such special global
cancellations.

For example the exact pair equations produce, for every triple,

```text
t_i G_jk t_jk - t_j G_ik t_ik + t_k G_ij t_ij=0,
```

with small residue coefficients and `G_ij` equal to a prescribed
intersection product of the `n_S` times a small correction. These
changing intersection products cannot generally be expressed as
polynomials in the row norms alone. Studying their simultaneous unit
values remains a concrete possible route.

What has been proved here is narrower and exact: norm identities and
all their generic coherent local powers are already accounted for by
the row-coordinate substitution and pair minors; bounded-degree row
norm auxiliaries and generic common norm denominators do not provide
a positive leading exponent margin. No incompatibility theorem for
the entire norm-block system, and no uniform endpoint bound, is claimed.
