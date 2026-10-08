# Seven-point inverse metric: square stars and exact triple contents

This is an exact arithmetic supplement to
[the low-moment lattice note](eight_point_low_moment_lattice.md). It does
not prove a radius-height inequality. It replaces the rational determinant
of the reconstructed metric by two explicitly defined integer gcds and
locates all possible denominator primes of that determinant.

Let `B` be an integer basis matrix of the saturated rank-two low-moment
kernel of seven distinct rational unit-circle points. Write its rows as
`b_i`, its primitive minors as `p_ij=det(b_i,b_j)`, and
`D_i=prod_(j!=i)p_ij`, with antisymmetric indices. The inverse construction
gives a positive rational quadratic form `q` with square determinant and
coprime positive integer values `q(b_i)=c_i`, satisfying
`|D_i|=d0 c_i^2`. All minors are nonzero.

## 1. The common star content is a square

In fact every `D_i` is negative and `d0=s^2` for an integer `s>0`.
The inverse construction says that the seven stars have one sign. Their
product is `(-1)^21 prod_(i<j)p_ij^2`, so this sign is negative. Moreover,

```text
d0^7 prod_i c_i^2 = prod_(i<j)p_ij^2.
```

Taking a valuation shows that `v_p(d0)` is even at every prime. Thus

```text
h_i=s c_i=sqrt(-D_i) is an integer,
s=gcd_i h_i,
prod_i h_i=prod_(i<j)|p_ij|.
```

In particular every local star degree `sum_(j!=i)v_p(p_ij)` is even.
This is stronger than merely requiring all star degrees to have a common
parity.

## 2. An integral triple formula for the metric scale

For a triple `i<j<k`, define

```text
X=c_i p_jk^2,  Y=c_j p_ik^2,  Z=c_k p_ij^2,
F_ijk=2(XY+XZ+YZ)-X^2-Y^2-Z^2,
T_ijk=|p_ij p_ik p_jk|.
```

Then

```text
F_ijk=4 det(q) T_ijk^2 > 0.                         (1)
```

To verify (1), take `b_i=(1,0)`, `b_j=(0,1)` by a rational change of
coordinates and write `b_k=(x,y)` and `q=A x^2+2Bxy+C y^2`.
The left side becomes `4(AC-B^2)x^2y^2`. Both sides have the same
transformation weight, proving the general identity.

Since `det(q)` is a rational square, every positive integer `F_ijk` is
itself a square: a rational number whose square is an integer is an
integer. Define

```text
G=gcd_(i<j<k) T_ijk,
R=gcd_(i<j<k) sqrt(F_ijk).
```

The exact identity, including every content factor, is

```text
2 sqrt(det(q)) = R/G.                              (2)
```

Here `R/G` need not be reduced. Indeed the rational gcd of a common
rational multiple of a finite integer list is that multiple times its
integer gcd; apply this to the positive square roots in (1).

Consequently the squared chord and intrinsic least squared radius are

```text
delta_ij^2=R^2 p_ij^2/(G^2 c_i c_j),
N=lcm_(i<j) [G^2 c_i c_j / gcd(G^2 c_i c_j,R^2 p_ij^2)].   (3)
```

Equivalently, at every prime, put `e_ij=v_p(p_ij)`,
`a_i=v_p(c_i)`, `g=v_p(G)`, and `r=v_p(R)`. Then

```text
g=min_(i<j<k)(e_ij+e_ik+e_jk),
v_p(N)=max(0, max_(i<j)(a_i+a_j-2e_ij)+2g-2r).      (4)
```

This formulation uses only the primitive Plucker coordinates and integer
square roots/gcds. No basis reduction, rational fitting, or omitted
exceptional prime is needed to compute it.

## 3. Only split primes can survive in the denominator of R/G

In lowest terms, the denominator of `R/G` is supported on primes
`p=1 mod 4`. This is a restriction on the inverse metric scale, in
addition to the usual split-prime support of the circle conductor.

Fix a prime. Primitivity gives some minor that is a local unit. Use its
two rows as a basis over `Z_p`. Then

```text
q(x,y)=A x^2+Txy+C y^2,
A,C in Z_p,    4det(q)=4AC-T^2 is a rational square.
```

The coordinate change has unit determinant, so it does not change the
valuation of `sqrt(det(q))`. If `v_p(T)<0`, dividing the last expression
by `T^2` gives a square equal to `-1+4AC/T^2`.
For an odd prime `p=3 mod 4` its reduction is the nonsquare `-1`, a
contradiction. For `p=2` it is `-1 mod 16`, again a contradiction.
Thus `T` is integral at these primes. At two it must furthermore be even:
an odd `T` would make the integral square `4AC-T^2` equal to `3 mod 4`.
It follows that `sqrt(det(q))` is integral at every nonsplit prime.

The un-reduced gcd `G` can nevertheless contain such primes; they must
cancel against `R`. For example, cotangents `0,1,...,6` give `G=10`.

There is also a concrete interpretation of `p|G`: reduce the seven rows
modulo `p`. Since they span a two-dimensional space, at least two
projective directions occur among their nonzero reductions. There are
exactly two such directions if and only if `p|G`. Indeed a triangle
product is a unit precisely when its three rows have three distinct
projective directions. Zero reduced rows do not count as directions.

## 4. Scope and verification

These identities expose the exact remaining cancellation: an upper bound
for `N` requires control of `R` relative to `G`, the values `c_i`, and the
minor sizes. No bound of exponent below `63/22`, or any improvement in the
uniform endpoint, follows from the formulas alone. In particular the
split-prime restriction does not remove the fair-profile conductor mass.

[check_seven_local_gcd.py](check_seven_local_gcd.py) computes the primitive
minors directly from the rational Gale matrix, checks the new square and
gcd formulas against the original circle chords, and checks all prime
valuations, including two, on a finite collection of exact examples.
