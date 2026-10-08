# Independent audit of the polarized support-multiplicity bound

The argument in
[`polarized_coefficient_support_multiplicity.md`](polarized_coefficient_support_multiplicity.md)
is valid.  It closes the cross-row-weight gap left by the invariant argument
in [`sl2_invariant_average_cut_width.md`](sl2_invariant_average_cut_width.md):
for every nonzero homogeneous Gale pullback of coordinate degree `D`,

```text
W(f) >= 320D,
sum_(127 cuts S) nu_S(F) <= 53D.                    (1)
```

No step requires the polynomial to remain in one row-weight component.
The essential hypotheses are the diagonal multiplicity and reciprocal
support symmetry of the complete pullback after all cancellations.

## 1. Polarization preserves multiplicity

For individual degree bounds `d_i`, polarization sends the monomial
`z^alpha` to

```text
product_i e_(alpha_i)(x_(i,1),...,x_(i,d_i))
          /binomial(d_i,alpha_i).                   (2)
```

Restriction to each group diagonal takes (2) back to `z^alpha`.  More
strongly, on separately symmetric multilinear polynomials the diagonal
restriction is injective in every total degree: the products of elementary
symmetric functions form a basis, and their diagonal restrictions are
distinct monomials multiplied by nonzero binomial coefficients. Translation
by one preserves separate symmetry and multilinearity. Diagonal restriction
commutes with that translation and preserves total homogeneous degree.
Applying injectivity to every homogeneous part of the Taylor
expansion proves that the first nonzero Taylor degree of the polarization is
exactly that of the original polynomial.  The weaker inequality used in the
main note is therefore safe.

Characteristic zero, or invertibility of the relevant binomial
coefficients, is needed here.  That hypothesis is present.

## 2. The signed coefficient indicator

Write the multilinear polarization as `P(x)=sum_A C_A x_A` and set

```text
Q(x)=sum_A (-1)^|A| C_A x_A product_(j notin A)(1-x_j). (3)
```

At the Boolean vertex `1_B`, every summand except `A=B` contains either a
zero `x` factor or a zero `(1-x)` factor.  Thus

```text
Q(1_B)=(-1)^|B|C_B.                                 (4)
```

The coefficient calculation is also exact:

```text
[x_T]Q=(-1)^|T| sum_(A subset T) C_A
      =(-1)^|T|P(1_T).                              (5)
```

If `P` has multiplicity at least `m` at the all-ones point and
`|T|>N-m`, the displacement from all ones to `1_T` uses fewer than `m`
variables.  Every multilinear Taylor monomial of degree at least `m` needs
at least `m` distinct variables, so it vanishes on that displacement.
Equation (5) then gives `deg Q<=N-m`.  This reasoning does not confuse
coefficient support with evaluation support; (3)-(5) construct a polynomial
whose Boolean evaluations are the desired coefficients.

Separate symmetry of `P` makes `Q` separately symmetric.  Replacing
`e_alpha(x_i)` by `binomial(t_i,alpha_i)` preserves individual and total
degree and evaluates `Q` on a Boolean vertex having group counts `t`.
Together with (2) and (4), its nonzero grid locus is exactly the original
coefficient support.

## 3. The grid-distance induction

The endpoint-slice proof of the polynomial grid lemma retains all degree
bounds.  If `a` and `b` are the first and last nonzero slices in the final
coordinate of length `d`, the identically zero exterior slices give the
factor

```text
product_(0<=j<a or b<j<=d)(t_n-j),
```

of degree `r=a+d-b`.  Grid interpolation in the other variables justifies
that a slice vanishing at every grid point vanishes as a polynomial.  After
division, both endpoint slices are nonzero, the remaining individual degree
bounds persist, and the total degree is at most `k-r`.  Induction costs at
most `(k-r)/2` in the other coordinates; averaging the two endpoint costs
adds exactly `r/2`.  Coordinates with `d_i=0` merely repeat the same corner
and can be omitted.  This proves the stated `k/2` bound without a hidden
genericity assumption.

## 4. Centering and the factor of 32

For homogeneous support of total degree `E` inside the box
`0<=alpha_i<=d_i`, with `sum d_i=2E`, the corner associated to `S` has

```text
dist_1(B(S),Supp(f))
 =d(S)+E-2 max_alpha alpha(S)
 =E-2h(S).                                         (6)
```

Combining the grid bound `(2E-m)/2` with (6) gives
`average_S h(S)>=m/4`.  If the support is invariant under
`alpha -> d-alpha`, its centered exponent vectors have coordinate sum zero,
so `h(S)=h(S^c)`.  The empty and full subsets contribute zero.  Summing over
all `2^n` subsets and then identifying complementary nontrivial cuts gives

```text
sum_(unoriented nontrivial cuts) h(S) >= 2^(n-3)m.  (7)
```

For `n=8`, the coefficient is `32`.  The normalization and the passage from
oriented subsets to 127 unoriented cuts therefore have the correct factor.

## 5. Gale application

An actual degree-`D` pullback `f=F(R_I(z))` is homogeneous of total degree
`16D`, lies in the individual-degree box `d_i=4D`, and satisfies

```text
sum_i d_i=32D=2*16D.
```

The exact inversion identity for every `R_I` makes the surviving support of
any sum invariant under `alpha -> 4D*1-alpha`.  Each coordinate is a unit
monomial times a five-point Vandermonde, hence belongs to the tenth power of
the diagonal maximal ideal.  Every product of `D` coordinates belongs to
its `10D`-th power, and sums cannot lower that multiplicity.  For a nonzero
pullback, (7) gives `W>=32*10D=320D`.  Substitution in the exact cut budget

```text
sum_S nu_S(F)=373D-W(f)
```

proves (1).

The nonzero-pullback clause is essential only to exclude coordinate-ring
relations, whose pullback is the zero polynomial and has no finite Newton
support.  The theorem bounds every surviving polynomial, but the separate
requirement that a finite family have no common zero on the positive actual
circle locus is still needed for the proposed endpoint contradiction.

This conclusion is consistent with the earlier generic Hermitian and moving
polynomial multiplicity barriers.  Its additional contribution is the exact
Newton-support duality for the actual all-row-weight Gale pullback and the
resulting sharp sum over the 127 primitive clean-cut filtrations.  It should
not be read as replacing the separate arithmetic hypotheses and scope audits
in those earlier arguments.

## 6. Exhaustive small tests

The independent checker
[`check_polarized_support_multiplicity_independent.py`](check_polarized_support_multiplicity_independent.py)
tests the proof outside the Gale-specific family.

- All `19682` nonzero bidegree-at-most `(2,2)` polynomials over `F_3` were
  polarized.  Their multiplicities before and after polarization agree;
  the signed indicator has degree at most `N-m`; its Boolean support is the
  coefficient support; and its coefficients retain separate symmetry.
- The grid-distance inequality was checked for all `19682` polynomials on
  the `3 by 3` grid and all `6560` multilinear polynomials on the Boolean
  three-cube.
- The centered corner-distance identity and reciprocal width inequality
  were checked for all `4912` nonzero support-symmetric homogeneous
  polynomials in the box `(1,1,1,1)` over `F_5`, and all `24564` in the box
  `(2,2,2)`.  Coefficients paired by reciprocal support were otherwise
  independent, allowing cancellations beyond reciprocal coefficient
  equality.

The grid points remain distinct and all polarization denominators are units
in these finite fields, so the same proofs apply to the tested instances.
No counterexample or normalization discrepancy was found.
