# Fixed additive half-character polynomials: the generic content budget

The multiplicative reduction extends to a precise additive obstruction.
Every fixed homogeneous polynomial in normalized balanced half-characters,
after taking a real or imaginary projection and performing its formal
cancellations, has a source-prime clearing cost at least `32k h` on the
full 127-cut profile, where `k` is its order of vanishing at coincident
source angles. The endpoint integer-gap test requires only `127k h/4`.
Thus this entire fixed-polynomial, allocation-content test has a positive
gap `k h/4`. Products of four disjoint pair differences attain the bound.

This is an application of the existing
[support-multiplicity theorem](polarized_coefficient_support_multiplicity.md)
to a different arithmetic normalization, not a new proof of that theorem.
The new points are the exact odd-twist clearing rule, its compatibility
with arbitrary nested allocations, and the additive bilinear reduction.
Extra divisibility at the actual value, or height-dependent coefficients
chosen from the character magnitudes, is **not** bounded by this result.
No endpoint configuration or uniform arc theorem is obtained.

## Exact bilinear reduction

Use the allocation notation of
[the torus-relation note](balanced_half_character_torus_relations.md).
The connection with source angles throughout this note assumes one
common source Gaussian unit after complete common-factor removal.
For independent row units, the normalized allocation characters acquire
character-dependent root-of-unity factors; the common-sign assertion
below and its diagonal Taylor bound are not asserted in that setting.
The literal allocation-product identities themselves remain valid.
For two balanced vectors `lambda,mu` and either `sigma=+1` or `-1`, put

```
q_p = (c_lambda(p)+sigma*c_mu(p))/2,
B = product_p pi_p^max(q_p,0) bar(pi_p)^max(-q_p,0),
C = product_p p^[(|c_lambda(p)|+|c_mu(p)|-2|q_p|)/2].
```

The `q_p` are integers because every balanced character has the same
allocation parity at each prime. With `A_mu^(+)=A_mu` and
`A_mu^(-)=bar(A_mu)`, the literal identity is

```
A_lambda A_mu^(sigma) = C B^2.                         (1)
```

Write `B=x+iy`. Both additive projections and the positive defect are
therefore exact:

```
product + conjugate(product) = 2C(x^2-y^2),
product - conjugate(product) = 4iCxy,
|product| - Re(product)      = 2Cy^2.                   (2)
```

The ordinary coordinate gcd of `C B^2` is exactly `C`: `B` is
conjugate-primitive and has odd norm. Individual scalar projections
can have additional factors, visibly `x` or `y` in (2). Removing the
rational content of the complex product does not justify deleting
those value-dependent factors from a different sum.

The source row vector of `B` is `(lambda+sigma*mu)/2`. Thus (2)
recovers an ordinary scalar row-character gap. This holds for every
pair, not only the distance-two pairs treated in the parity note.
It includes the zero vector, where `B=1` and the positive defect is zero.

## Clearing a fixed polynomial with the odd twist intact

Include all 70 balanced sign vectors, with `r_-lambda=bar(r_lambda)`.
As in the odd-twist note, set

```
c_lambda(p)=sum_i lambda_i a_i(p),
t=product_(p: sum_i a_i(p) odd) p,
Norm(A_lambda)=t h_lambda^2,
r_lambda=A_lambda/(h_lambda sqrt(t)).
```

Let `P` be a fixed homogeneous polynomial of degree `d` in the `r`'s,
with coefficients in `Z[i]`. Fixed rational coefficients only require
one fixed clearing factor. Reduce all formal torus relations and combine
equal monomials first, obtaining

```
P = sum_(L in E) b_L exp(i L dot theta),
L = (1/2)sum_(j=1)^d lambda_j,
sum_i L_i=0,
L_i == d/2 (mod Z).                                   (3)
```

There may be one common sign multiplying all `r_lambda` when a single
unwrapped source arc is used. Indeed changing source angle lifts adds
`pi sum lambda_i n_i`, whose parity is `sum n_i`, independent of lambda.
For a homogeneous degree this contributes one harmless common sign
`epsilon^d`. No individual character signs are assumed independently.

For each split prime, define the nonnegative integer

```
m_p = max_(L in E) |2 L dot a(p)|,
kappa_p = d sum_i a_i(p) mod 2,    in {0,1},
t_d = product_p p^kappa_p,
Q_P = product_p p^[(m_p-kappa_p)/2].
```

Every number in the maximum has parity `kappa_p`, so all exponents
are nonnegative integers. Then

```
U_P = Q_P sqrt(t_d) P(r) belongs to Z[i],
log(Q_P sqrt(t_d)) = (1/2)sum_p m_p log p.              (4)
```

This is a literal termwise clearing statement, even when several
nested threshold layers belong to one prime. To verify it, a source
monomial with signed exponent `c=2L dot a` has, after multiplication
by `sqrt(t_d)`, local Gaussian valuations

```
v_pi = (c+kappa_p)/2,
v_barpi = (-c+kappa_p)/2.
```

Multiplication by `Q_P` changes these to `(m_p+c)/2` and
`(m_p-c)/2`. They are nonnegative integers. This derivation keeps all
rational content and never asserts that an odd quarter-root is Gaussian.

For a real or imaginary projection, the surviving support is reciprocal:
`E=-E`. Each extreme valuation then occurs in the support. Thus (4)
is the exact minimal rational clearing of all individual supported
monomials after the twist is removed. It need not be the reduced
denominator of their final evaluated sum. Special arithmetic
cancellation can increase the final content and remains a possible
source of new information.

## Every fixed vanishing order has the full-cut margin

Suppose the reciprocal Laurent polynomial in (3) is nonzero and vanishes
to exact order `k>=1` at the diagonal. Take an integer `d` as above and
form the ordinary polynomial

```
f(z) = (product_i z_i^(d/2)) sum_L b_L z^L.
```

Its exponents `alpha=L+(d/2)1` are integral, between zero and `d`;
it has total degree `4d`; and its support is invariant under
`alpha -> d1-alpha`. The monomial prefactor is nonzero at the diagonal,
so the multiplicity stays `k`. The existing support-multiplicity theorem
therefore gives

```
h(S) = max_(L in E) |sum_(i in S) L_i|,
sum_(127 unoriented cuts S) h(S) >= 32k.               (5)
```

The theorem uses all cuts and applies to each actual support; it makes
no symmetry assumption about a distinguished row.

At equal formal cut weight `h`, (4) and (5) give clearing modulus
logarithm at least `32k h`. On an arc of width `Delta`, the Taylor bound
is `|P(r)|<=K_P Delta^k`, with a constant depending only on the fixed
polynomial. If `U_P` is nonzero, Gaussian integrality supplies

```
log(Q_P sqrt(t_d)) >= k W/4 - log(K_P C^k),
                         Delta<=C exp(-W/4).           (6)
```

The full-cut profile has `W=127h`. Its clearing budget exceeds the
leading requirement in (6) by at least `k h/4`. Consequently these
inequalities do not exclude that profile, at any fixed degree or order.
If `P(r)=0`, the nonzero-integer argument supplies no inequality.
With literal independent split-prime blocks of weights `h+O(1)`, the
same conclusion holds for every fixed `P`, with an `O_P(1)` error.

The coefficient is sharp. For four disjoint pairs covering the eight
rows, take

```
P(theta)=exp(-i sum_i theta_i/2)
         product_(j=0..3)(exp(i theta_(2j))-exp(i theta_(2j+1))).
```

This is a real linear combination of 16 balanced half-characters,
vanishes to order four, and has

```
h(S)=(number of matching edges separated by S)/2,
sum_S h(S)=4*64/2=128=32*4.                            (7)
```

Its exact arithmetic reduction is a rational-content multiple of the
product of the four primitive pair imaginary coordinates. Its gap is
the product of four already available pair gaps.

For comparison, simultaneously clearing the full linear family of all
70 half-characters has `m_p=r` on a cut of minority size `r`. Its
modulus logarithm is

```
(8*1+28*2+56*3+35*4)h/2 = 186h.
```

This larger common denominator is not an improvement on selecting a
smaller support. The sharp support bound (5) already permits the most
favorable formal cancellations and support choices.

The [exact checker](check_additive_balanced_character_generic_content.py)
verifies nested-prime bilinear identities, their ordinary contents,
the odd-twist valuation clearing formula, and the sharp 16-term
matching polynomial and cut census. The unrestricted inequality (5)
is supplied by the cited theorem, not inferred from a finite search.
