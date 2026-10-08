# Exact cotangent spectra and the missing height restriction

The integer cotangent matrix has a fixed characteristic polynomial,
independent of the circle configuration. This supplies exact rational
Gaussian spectral projectors with controlled denominators. It does not
prove the lcm height target in
[integer_cotangent_lcm_height_target.md](integer_cotangent_lcm_height_target.md).
The matrix is highly nonnormal on a short arc, and its extreme
projectors are not positive operators.

## 1. Consistent signs and the spectral identity

Use `n` circle labels, including the anchor at infinity. Thus there
are `n-1` finite cotangents `X_i`, and put

```text
Q_0i=X_i,          Q_i0=-X_i,
Q_ij=(X_i X_j+L^2)/(X_i-X_j),          Q_ji=-Q_ij.
```

These signs are simultaneous: changing only the anchor signs would
destroy the following identity. Let

```text
s_i=sum_(j!=i) Q_ij,
B_ii=s_i,         B_ij=Q_ij+iL  (i!=j).                (1)
```

Then

```text
det(tI-B)=product_(r=0)^(n-1) [t-iL(2r-(n-1))].        (2)
```

The proof does not require integrality. Put

```text
r_0=1,      r_i=(X_i-iL)/(X_i+iL),
w_i=1/product_(j!=i)(r_i-r_j).
```

All nodes are distinct, and

```text
Q_ij=iL(r_i+r_j)/(r_i-r_j).
```

On polynomials of degree less than `n`, the operator `z d/dz` has
eigenvalues `0,...,n-1`. In the basis of values at the nodes, its matrix
`D` has entries

```text
D_ij=r_i w_j/[w_i(r_i-r_j)]       (i!=j),
D_ii=sum_(j!=i) r_i/(r_i-r_j).
```

The off-diagonal formula follows by differentiating the Lagrange
interpolation polynomials. With `W=diag(w_i)`, identity (1) is precisely

```text
B=2iL W D W^(-1)-iL(n-1)I.
```

This proves (2). An eigenvector for `lambda_r=iL(2r-(n-1))` is

```text
v_r=(w_i r_i^r)_i.                                    (3)
```

## 2. Exact integral projectors and their extreme forms

If the cotangents are integers, `B` is a Gaussian integer matrix.
Since its eigenvalues are distinct, its spectral projector is

```text
P_r=product_(s!=r)(B-lambda_s I)/(lambda_r-lambda_s).
```

Its denominator divides, in `Z[i]`,

```text
D_n=(2iL)^(n-1)(n-1)! .                               (4)
```

Indeed the actual denominator is
`(2iL)^(n-1) (-1)^(n-1-r) r!(n-1-r)!`, and the quotient
of `(n-1)!` by the two factorials is an integer. Thus every `D_n P_r`
is Gaussian integral.

The extreme projectors are particularly explicit. Put `U=product_i r_i`.
Inverting the eigenvector matrix in (3) by polynomial interpolation gives

```text
(P_(n-1))_ij = w_i r_i^(n-1),
(P_0)_ij = (-1)^(n-1) U w_i/r_j.                      (5)
```

For example, the first identity uses that the leading coefficient of
`product_(l!=j)(z-r_l)` is one. Its column sum is one because
`sum_i w_i r_i^(n-1)=1`.

Unit-circle conjugacy is retained exactly. If `J=diag(r_i)`, then

```text
B=J conjugate(B) J^(-1),
P_r=J conjugate(P_(n-1-r)) J^(-1).                    (6)
```

In particular opposite extreme projectors have conjugate diagonal
entries. This is a diagonal unitary conjugacy, not self-adjointness or
positivity of the projectors.

## 3. Why the fixed spectrum does not control entry height

The Hermitian part is

```text
B+B*=2 diag(s_i).                                    (7)
```

It is generally large and indefinite. In fact, with the Frobenius norm,

```text
||B||_F^2-sum_r |lambda_r|^2 = 2 sum_i s_i^2.          (8)
```

To prove this, write `B=S+K`, with `S=diag(s_i)` Hermitian and `K`
skew-Hermitian. Since `K` has zero diagonal, `tr(SK)=0`. Therefore
`tr(B^2)=tr(S^2)+tr(K^2)` and
`||B||_F^2=tr(S^2)-tr(K^2)`. The eigenvalues in (2) are purely
imaginary, so `tr(B^2)=-sum_r |lambda_r|^2`, giving (8).

Already for two labels,

```text
B=((X, X+iL), (-X+iL, -X)),
B^2=-L^2 I,
P_+=(B+iL I)/(2iL).
```

The spectrum remains `{-iL,iL}`, and the projector denominator is fixed
when `L` is fixed, while both `||B||` and `||P_+||` grow without bound
as the integer `X` grows. The projector's diagonal entries are
`1/2-iX/(2L)` and `1/2+iX/(2L)`, so it is not Hermitian. Thus a
contractivity or positivity bound on the projectors cannot be inserted
into the height argument.

Unbounded integer spectral matrices exist with any fixed number of
labels, not just two. For `k=n-1`, choose a positive integer `L` divisible
by every integer `1,...,k-1`, and set

```text
X_j=L(t+j),          j=0,...,k-1,       t=1,2,... .     (9)
```

Every pair cotangent is the integer

```text
L [(t+i)(t+j)+1]/(i-j).
```

Hence (2), (4), and (6) all hold with fixed `n,L` and arbitrarily large
entries. This is not an endpoint counterexample. The rational ratios
in (9) equal those of `t+j+i`. Gaussian gcds between different such
numerators divide the bounded nonzero integers `i-j`; their individual
conjugate gcds have norm at most two. The least primitive radius is
therefore `R asymp_k t^k`, with constants independent of `t`, by the
Gaussian denominator lcm formula. The normalized arc including the
infinity anchor has size asymptotic, up to fixed constants, to
`t^(k/2-1)`. For `k>=3` it grows.

The fixed-`L` Pell triple in the height-target note is a different,
stronger warning: with four labels it also satisfies every spectral
identity here while its endpoint constant tends to zero.

## 4. The remaining arithmetic

The exact least radius is still `sqrt(N)`, where `N` is the all-edge
ordinary lcm from the height-target note. Neither the eigenvalues nor
the common projector denominator in (4) bound `N` from below in terms
of `min X_i/L` with exponent four.

The eigenvectors and their primitive integral normalizations retain
that arithmetic: they involve the actual rational nodes and their
Vandermonde denominators. A possible continuation would need a new
height restriction on these eigenvectors, compatible with the extreme
projector formulas (5). It cannot use positivity or bounded operator
norms that fail even in the two-label example. No such restriction,
finite-label lcm inequality, or improvement in the general growth bound
is proved here.

## 5. Primitive left/right eigenvectors have bounded pairings

Choose primitive Gaussian integral right and left eigenvectors `u_r`
and `v_r` for each eigenvalue, and put `h_r=v_r^t u_r`. Here left
eigenvectors use transpose, not conjugate transpose. Distinct eigenvalues
give `v_r^t u_s=0` for `r!=s`, while `h_r!=0` because the spectrum is
simple. The exact rank-one projector is

```text
P_r=u_r v_r^t/h_r.
```

Let

```text
Delta_r=product_(s!=r)(lambda_r-lambda_s)
       =(2iL)^(n-1)(-1)^(n-1-r) r!(n-1-r)!.
```

Then

```text
h_r | Delta_r       in Z[i].                          (10)
```

Indeed `Delta_r P_r` is Gaussian integral. The entries of
`u_r v_r^t` have Gaussian gcd one: multiply Gaussian Bezout identities
for the primitive coordinates of `u_r` and `v_r`. Applying that same
integer linear combination to `Delta_r u_r v_r^t/h_r` proves that
`Delta_r/h_r` is Gaussian integral.

With `U` and `V` the matrices of these right and left eigenvectors,

```text
V^t U=diag(h_r),
det(U) det(V)=product_r h_r | product_r Delta_r.       (11)
```

Thus these determinants have bounded Gaussian norm in terms of `n,L`.
This is a determinant bound, not a bound on eigenvector height.

## 6. Every conductor valuation profile cancels outside `2L`

The weighted Vandermonde eigenbasis makes the scale limitation in
(11) explicit. Fix a Gaussian prime `pi` not dividing `2L`. Write
`v` for its valuation. Since every `Q_ij` is integral,

```text
v(r_i-r_j)=min(v(r_i),v(r_j)).                          (12)
```

When the two node valuations differ, (12) is automatic. When they
are equal, a larger valuation of the difference would leave the sum
with the smaller valuation, since the residue characteristic is not
two. The identity `Q_ij=iL(r_i+r_j)/(r_i-r_j)` would then contradict
integrality. Multiplying every node by one common nonzero Gaussian
rational number does not change `B`, and permits a common shift of
all these valuations.

Sort the node valuations as `e_1<=...<=e_n`, allowing repeated or
negative values. Formula (12) gives

```text
v(w_i)=-sum_(j<i)e_j-(n-i)e_i.
```

In raw eigenvector number `r`, the coordinate valuation is therefore

```text
f_i=(r-n+i)e_i-sum_(j<i)e_j.
```

Its successive difference is

```text
f_(i+1)-f_i=(r-n+i+1)(e_(i+1)-e_i).
```

The coefficient changes from nonpositive to nonnegative at
`i=n-r-1`. The minimum, with the natural endpoint interpretation,
is exactly

```text
c_r=-sum_(j=1)^(n-r-1)e_j.
```

The raw weighted-Vandermonde determinant has valuation

```text
v(det(v_0,...,v_(n-1)))=-sum_(i=1)^n (n-i)e_i.
```

Indeed, the product of the weights contributes minus twice the
valuation of the Vandermonde determinant. But the sum of all column
minima is precisely

```text
sum_(r=0)^(n-1)c_r=-sum_(i=1)^n (n-i)e_i.
```

Making each rational eigenvector Gaussian integral and primitive
subtracts exactly `c_r` from its coordinate valuations. Consequently

```text
v(det U)=0                                            (13)
```

at every Gaussian prime not dividing `2L`. Thus the Gaussian prime
support of `det U` is contained in the support of `2L`, not merely
the larger spectral-gap product. No assumption on a two-level profile,
a fair cut, prime-power multiplicity, or label ordering is needed.

For a two-level cut with `s` valuations equal to `a>0` and the others
zero, both the raw determinant valuation and the sum of primitive
column contents equal `-s(s-1)a/2`. This illustrates the exact
cancellation on every fair cut, with its full conductor exponent
retained. The general calculation shows that mixed layers do not
restore a positive conductor contribution.

This explains why bounded spectral gaps and bounded eigenbasis
determinants alone do not give the missing radius power. A useful
continuation must control the archimedean sizes of the primitive
eigenvectors, or their simultaneous integer realization, rather than
charging the raw Vandermonde determinant without its column contents.

There is a useful denominator corollary. At every Gaussian prime not
dividing `2L`, the integral matrix `U` is unimodular, so `U^(-1)` is
integral and each of its rows is primitive over the local ring. Its
row numbered `r` is a left eigenvector. The globally primitive left
eigenvector differs from that row by a local unit. Consequently

```text
v(det V)=0,          v(h_r)=0,          P_r is integral locally
                                      whenever pi does not divide 2L.
```

Thus all spectral-projector denominators are supported on `2L`;
factorial primes outside that support cancel. Combining this with
(10) yields the common integer denominator

```text
(2L)^(n-1) * product_(p | 2L) p^(v_p((n-1)!)).
```

This improves the denominator accounting in (4), without bounding
the sizes of the projectors or the primitive eigenvectors.

Run `python3 docs/check_integer_cotangent_spectral.py` for exact rational
checks on twelve actual cliques, including the Pell example, a
five-point clique, and an eight-label clique. It checks the complete
eigenbasis, the extreme projectors and their integral denominators,
unit-circle conjugacy, and (8). It normalizes the actual rational right
and left eigenvectors to primitive Gaussian integer vectors, checks
(10), and computes the exact primitive right-eigenbasis determinant
by the Vandermonde formula to verify its prime support in `2L`.
It checks the same support restriction for the left determinant and
every primitive eigenvector pairing.
It also checks all column minima and their sum on 1,281 ordered
valuation profiles, including negative and repeated levels. Each
circle sample passes the independent all-edge least-radius and
reanchoring checker. These finite checks supplement the general
proofs above.
