# A degree-versus-divisibility obstruction for general half-angle forms

This is a restriction on a proof method, not a proof or disproof of the
uniform endpoint bound. It extends the bracket-polynomial calculation in
`plucker_compatibility_route.md` to arbitrary multihomogeneous polynomials
in the real and imaginary parts of Gaussian half-angle factors. In
particular, the polynomial need not be an invariant of projective changes
of coordinates or a polynomial in brackets alone.

The statement concerns the exponent supplied by generic prescribed prime
divisibility and a Taylor estimate at the short arc. It does not bound
extra divisibility at special arithmetic values, force the evaluated
polynomial to be nonzero, or rule out other uses of metric identities.

## 1. A polynomial inequality

Let `F(X_1,Y_1,...,X_m,Y_m)` be a nonzero polynomial over `Q(i)`,
homogeneous of degree `d_j` in each pair `(X_j,Y_j)`. Put

```text
D = sum_j d_j,
f(t_1,...,t_m) = F(1,t_1,...,1,t_m),
v = ord_(0,...,0) f.
```

Thus `v` is the least total degree of a monomial occurring in `f`, not
an asserted order at a particular tuple of arithmetic points. For a
subset `S` of `[m]`, define `nu_S^+` and `nu_S^-` to be the largest
integers such that, respectively,

```text
f belongs to (t_j-i : j in S)^(nu_S^+),
f belongs to (t_j+i : j in S)^(nu_S^-).
```

Equivalently these are generic vanishing orders along the indicated
coordinate subspaces. The empty-set orders are zero. Then

```text
2^(-m) sum_(S subset [m]) (nu_S^+ + nu_S^-)
    <= (D-v)/2.                                      (1)
```

For rational coefficients, complex conjugation gives
`nu_S^+=nu_S^-=:nu_S`, so the right formulation is

```text
2^(-m) sum_S nu_S <= (D-v)/4.                         (2)
```

### Proof

Separate homogeneity implies that

```text
g(s) = product_j s_j^(d_j) f(1/s_1,...,1/s_m)
```

is a nonzero polynomial. Every monomial has degree `D-|alpha|` for a
monomial `t^alpha` of `f`, and hence `deg g=D-v`.

For `S subset [m]`, take the vertex with `t_j=i` on `S` and
`t_j=-i` on its complement. The two defining coordinate ideals use
disjoint sets of variables. Their powers have intersection equal to
their product; alternatively, every Taylor monomial has at least
`nu_S^+` degrees in the first set and at least `nu_(S^c)^-` in the
second. Consequently the multiplicity of `f` at this vertex is at least

```text
nu_S^+ + nu_(S^c)^-.
```

Coordinatewise inversion is a local isomorphism and permutes these
vertices, interchanging `i` and `-i` in each coordinate. The monomial
prefactor defining `g` is a unit there. Thus the summed multiplicities
are unchanged. The multiplicity form of Schwartz--Zippel
on the two-element grid `{i,-i}^m` gives

```text
sum_S (nu_S^+ + nu_(S^c)^-) <= (D-v) 2^(m-1).
```

Reindexing the complementary term proves (1).

The grid inequality used here is the multiplicity Schwartz--Zippel
lemma; see Dvir, Kopparty, Saraf and Sudan,
[Extensions to the Method of Multiplicities](https://arxiv.org/abs/0901.2529).
For completeness, an elementary characteristic-zero induction suffices.
Write a nonzero polynomial of total degree at most `d` as a polynomial
of degree `e` in its last variable, with nonzero leading coefficient `h`
of degree at most `d-e`. At each vertex `a` in the remaining variables,
differentiate in those variables to order `r=mult(h,a)` so that the
leading coefficient remains nonzero after evaluation. The resulting
univariate polynomial has degree `e`. At either last coordinate `b`,
the original multiplicity is at most `r` plus that univariate
multiplicity. Summing first over the two choices of `b`, then using
induction on `h`, gives
`2(d-e)2^(m-2)+e2^(m-1)=d2^(m-1)`.

## 2. The actual Gaussian model to which it applies

For each `S subset [m]`, choose a Gaussian block `H_S` such that the
blocks and all their conjugates are pairwise coprime. Omit the prime
above 2; it has no effect on a leading logarithmic exponent. Set

```text
G = product_S H_S,
A_j = product_(S containing j) H_S,
z_0 = conjugate(G),
z_j = conjugate(G) A_j/conjugate(A_j),
R = |G|,
W = log(R^2),
w_S = log N(H_S).
```

The `z_j` are actual Gaussian integers of common modulus `R`.
This is the independent-block model; no assertion of independence is
made for nested layers at a common prime in an arbitrary conductor.

Suppose the points are on an arc containing `z_0` of width at most
`C R^(-1/2)`. Choose signs of the `A_j` so that

```text
A_j = X_j+iY_j,
|Y_j| <= (C/2) R^(-1/2) |A_j|.                     (3)
```

Signs do not affect divisibility or the following estimates. If
`w_S=W/2^m+o(W)` for every `S`, then

```text
|A_j| = R^(1/2+o(1)),
|Y_j| = R^o(1).
```

For the fixed polynomial `F`, separate homogeneity and the definition
of `v` give

```text
|F(X_1,Y_1,...,X_m,Y_m)|
    <= R^((D-v)/2+o(1)).                            (4)
```

The constants can depend on `F,m,C`; the exponent statement does not
ignore their dependence if the polynomial is allowed to vary.

At a Gaussian prime dividing `H_S`, the linear forms
`X_j+iY_j` for `j in S` have the prescribed prime divisibility. At a
prime of its conjugate block, the corresponding forms are `X_j-iY_j`.
To make the affine-to-homogeneous step explicit, expand `f` in powers
of `t_j-i` for `j in S`. Every term has total degree at least
`nu_S^+` in those powers, and each exponent is at most `d_j`.
Homogenizing therefore introduces no negative powers of `X_j`;
`Y_j-iX_j=-iA_j` supplies the required block factors. The opposite
root gives `i conjugate(A_j)`. The shifted expansion uses integer
binomial coefficients, so clearing the coefficient denominators of
`F` once also clears this expansion.

Membership in the ideals defining `nu_S^+` and `nu_S^-` therefore
gives a divisor of this fixed integer multiple of `F`
whose logarithmic modulus is

```text
one half sum_S w_S (nu_S^+ + nu_S^-)
    <= W(D-v)/4+o(W).                              (5)
```

Here fixed primes in coefficient denominators only require a fixed
clearing factor. The pairing of `+` and `-` is unchanged if the two
roots are named in the opposite order.

Since `W=2 log R`, the largest exponent supplied by (5) is exactly the
upper exponent in (4). Thus this method supplies **no fixed positive
power saving in the radius** at the uniform profile. The divisor in
(5) is only the divisor guaranteed by the generic orders: special
values of `F` can have additional prime factors or can vanish.

## 3. Sharpness and normalization

For a single bracket `F=X_1Y_2-X_2Y_1`, we have `D=2`, `v=1`, and
`nu_S^+=nu_S^-=1` precisely when `{1,2} subset S`. Equality holds
in (1). Products of brackets give equality in every degree. Pure
imaginary-coordinate products have no such prescribed divisibility,
which is also consistent with (1). Norm polynomials
`X_j^2+Y_j^2` and their products also attain equality; their sharpness
expresses the exact norm-factor identity. A real-coordinate factor
`X_j` gives strict inequality.

The empty block is a common factor of all the `z_j`; it can be removed
without changing the conclusion. In the limiting uniform model put

```text
R' = R/|H_empty|,
rho = 2^(m-1)/(2^m-1).
```

Then the retained nonempty blocks have asymptotically equal weight,
`|A_j|=R'^(rho+o(1))`, and the inherited arc width is
`R'^(-rho+o(1))`. The archimedean exponent becomes `rho(D-v)`.
Because the empty-set orders are zero, deleting that term multiplies
the average in (1) by `2^m/(2^m-1)`. Formula (5) gives exactly the same
ceiling `rho(D-v)`. One must retain this improved arc width when
changing radius; replacing it by the weaker endpoint width would
misstate the comparison.

## 4. Scope for further work

This theorem explains why generic polynomial divisibility combined
with ordinary Taylor bounds does not yield the fixed exponent margin
needed by the current uniform-profile exclusion target. It includes
more forms than the Pluecker relations, but it is not a theorem that
all polynomial approaches fail. Potential extra inputs include
special-value cancellation, restrictions on simultaneous primitive
residues, polynomials with controlled changing degree and coefficients,
or a different estimate for the actual arithmetic tuple. None of
these extra inputs is proved here.

Two independent agent audits checked the reciprocal step, Gaussian
divisibility, and the change of radius. An additional exact Gaussian
integer polynomial calculation checked (1) on 360 products in two
through five variables, including linear isotropic factors, brackets,
coordinate factors and norms; 230 cases attained equality and none
violated the inequality. These finite checks supplement the proof.

The uniform bound remains unproved. This note is a prose argument,
not a Lean formalization.
