# Exact Smith form of the whole squared-distance matrix

The whole integer distance matrix has an exact Smith form. Its last
invariant factor retains the circle conductor, multiplied by an explicit
square. In the primitive common-unit model the first two factors remove
exactly the common primitive chord residue. This is an arithmetic identity
for actual circle points, not a proof of the eight-row source-gcd inequality
`A^6 F^3 <= exp(B) N` or a new bound for endpoint arcs.

The rank-three factorization and individual minor identities already appear
in [the distance-matrix audit](distance_laplacian_denominator_audit.md).
The result here retains the gcds of **all** minors, hence the complete
Smith form. The chord and triangle normalizations are those of
[the residue-content note](residue_content_descent.md) and
[the primitive-triangle theorem](affine_shape_radius_divisibility.md).

## 1. General integral circle tuples

Let `z_i=(x_i,y_i)`, `1<=i<=m`, be distinct integer points of norm
`N>0`, with `m>=3`. Define

```text
E_ij = |z_i-z_j|^2/2 = N-x_i x_j-y_i y_j,
d = gcd_(i<j) E_ij,
D_ijk = det(z_j-z_i,z_k-z_i),
g = gcd_(i<j<k) |D_ijk|.
```

All displayed gcds are positive. Every triangle determinant is nonzero,
since a line intersects the circle in at most two points. Then the three
nonzero Smith invariant factors of `E` are exactly

```text
d, d, N g^2/d^2.                                      (1)
```

In particular `d^3 | N g^2`. This assertion allows independent row units,
even norm, and arbitrary common Gaussian factors.

**Proof.** For the matrix `Y` with rows `(1,x_i,y_i)`,

```text
E = Y diag(N,-1,-1) Y^t.
```

Thus `E` has rank three, and its row/column three-minors are

```text
det E[I,J] = N D_I D_J.                               (2)
```

The first determinantal divisor is `d`. Every two-minor is divisible
by `d^2`, while the principal two-minors are `-E_ij^2`; their gcd is
exactly `d^2`. Hence the second determinantal divisor is `d^2`.
Equation (2) makes the third determinantal divisor exactly `N g^2`.
Successive quotients of these divisors prove (1).

One also always has `d | g`. Anchor at `z_1`. The Gram entries of the
differences are `E_1i+E_1j-E_ij`, hence are divisible by `d`. The Gram
determinant of any two differences is their determinant squared, so
`d^2` divides every anchored triangle determinant squared. Thus `d`
divides every anchored determinant, and then every triangle determinant
by the usual three-term anchored area identity.

If the source tuple is Gaussian-primitive, then

```text
gcd(d,N)=1,       d^3 | g^2.                           (3)
```

Indeed a primitive equal-norm tuple has odd norm with only split prime
factors. At `p=pi bar(pi) | N`, choose rows with allocations zero and
`v_p(N)`. In the identity
`2E_ij=2N-z_i bar(z_j)-bar(z_i)z_j`, one cross product is a `pi`-unit
and the other is divisible by `pi`. Since `p` is odd, `E_ij` is a
`p`-unit. This proves the coprimality; (1) then gives the second assertion.

## 2. A primitive normalization allowing all row units

For a Gaussian-primitive tuple there are integers `h,a>=1` and
`kappa in {1,2}` such that

```text
d = kappa h^2,       g = 2 kappa h^3 a,
Smith(E/d) = diag(1,1,(4/kappa) N a^2,0,...,0).          (4)
```

Here `kappa=1` exactly when the tuple has both coordinate-parity types
`(odd,even)` and `(even,odd)`. In that case `h` is odd. Otherwise
`kappa=2`.

To prove the square assertion at an odd prime not dividing `N`, put
`w=z_i-z_j`. Equal norms give the exact Gaussian identity

```text
bar(w) = -N w/(z_i z_j).
```

At a split prime both rows are units at both orientations, so the two
orientation valuations of `w` agree. At an inert prime a norm valuation
is automatically even. Thus every nonzero `E_ij` has even valuation at
each odd prime outside `N`. The primes of `N` were excluded by (3).

At two, opposite coordinate-parity types give an odd `E_ij`. If all
points have the same parity type, fix Gaussian generators of the split
primes of `N`. Their row products have the same parity type because a
generator and its conjugate agree modulo two. The independent row units
therefore differ only by signs. After a pair gcd is removed its two rows
are `epsilon H` and `+/- epsilon bar(H)`, with `Norm(H)` odd. Hence
`E_ij` is twice an odd integer times the square of either `Re(H)` or
`Im(H)`. Every such valuation at two is odd, and so is their minimum.
This proves `d=kappa h^2` with the stated parity rule.

All points satisfy `x_i+y_i=1 mod2`, so every triangle determinant is
even. Combining this with `d^3 | g^2` proves `2h^3 | g` when `kappa=1`
and `4h^3 | g` when `kappa=2`. Substitution in (1) gives (4).

## 3. The common-unit primitive residue dictionary

Assume in addition the primitive common-unit model. Put

```text
n_ij = Norm(gcd_G(z_i,z_j)),
|z_i-z_j|^2 = 4 n_ij t_ij^2,       t_ij in Z_(>0),
h = gcd_(i<j) t_ij.
```

The exact global content is

```text
d = 2 h^2.                                           (5)
```

At a conductor prime the extremal-allocation pair in the previous proof
has both `n_ij` and `t_ij` units at that prime. Away from the conductor every
`n_ij` is a unit. Thus, for every odd prime outside `N`, the minimum
valuation of `E_ij=2 n_ij t_ij^2` is twice the minimum residue valuation.
At two all `n_ij` are odd, so the same equality has the additional factor
two. This proves (5) at every prime, including all prime-power depths.

The exact common-unit triangle identity is

```text
|D_ijk| = 4 Norm(gcd_G(z_i,z_j,z_k)) t_ij t_ik t_jk.
```

Consequently `g=4 h^3 a` for a positive integer `a`. Substitution in
(1) gives the particularly simple primitive distance matrix

```text
B = E/(2 h^2),
Smith(B) = diag(1,1,2 N a^2,0,...,0).                  (6)
```

The integer `a` retains every extra triangle content, including accidental
content at conductor primes. It is not asserted to be coprime to `N`.
Formula (6) exposes the conductor in one Smith factor; replacing `a` by
one without justification would discard actual evaluated arithmetic.

## 4. What the endpoint hypothesis supplies

Suppose the tuple lies on an arc of length at most `C sqrt(R)`, where
`R=sqrt(N)`. Every entry satisfies

```text
0 <= B_ij <= C^2 sqrt(N)/(4 h^2).                      (7)
```

For an ordered triangle, the two consecutive arc gaps sum to at most
the arc length. Their product is at most one quarter of its square.
The circumradius identity therefore gives

```text
|D_ijk| <= C^3 N^(1/4)/8,
h^3 a <= C^3 N^(1/4)/32.                              (8)
```

Thus the Smith identity is compatible with the existing triangle bound.
Bounding three-minors by their entries likewise supplies only that same
type of radius-dependent estimate. A squareclass constraint on the last
Smith factor is not a height bound.

For the unit-free normalization in (4), the corresponding bounds are
`max(E/d)<=C^2 sqrt(N)/(2 kappa h^2)` and
`h^3 a<=C^3 N^(1/4)/(16 kappa)`.

In particular, (1)--(8) do not control the singleton and pair source
deletion contents `A,F`. The exact Smith factors compress those detailed
cut locations into `N` and the common triangle content. A successful
continuation through this matrix would need an additional arithmetic
inequality on its simultaneous subset contents or its actual short
integer row relations; the rank, complete Smith factors, and entry
bound alone supply no strict radius exponent gain.

## Verification

Run `python3 -B docs/check_primitive_distance_smith_form.py`. The checker
computes all two- and three-minor gcds on actual equal-norm Gaussian
tuples, including nested allocations, independent units, and nonprimitive
common factors. In the common-unit cases it computes literal Gaussian
source gcds and primitive residues and checks (5)--(6) and every triangle
identity. Separate rational half-angle fixtures have nontrivial common
residue `h=6` and mixed coordinate parity with `kappa=1,h=3`. These
finite checks verify the identities and do not claim
endpoint localization for the fixtures.
