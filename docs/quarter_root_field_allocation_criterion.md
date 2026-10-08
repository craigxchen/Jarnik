# Source allocations that collapse the quarter-root field

This note characterizes one restricted joint-field situation using the
actual prime exponents of eight Gaussian source points. It does not prove
a new endpoint inequality. It also distinguishes fixed character
representatives from independently conjugating them; those are different
field-minimization questions.

Use the literal common-unit factorization and the conjugate-primitive
characters `A_lambda` of
[the odd-twist note](odd_layer_character_quadratic_twist_obstruction.md).
Choose the 35 balanced sign vectors with `lambda_0=+1`, and write

```
c_lambda(p)=sum_i lambda_i a_i(p),
A_lambda=product_p pi_p^max(c_lambda(p),0)
                    bar(pi_p)^max(-c_lambda(p),0),
t=product_(p : sum_i a_i(p) odd) p.
```

The Gaussian prime representatives are fixed throughout; the displayed
products insert no independent Gaussian units or square-root signs.

Assume `t>1`. At primes outside `t`, every `c_lambda(p)` is even and
contributes no Gaussian squareclass. At a prime dividing `t`, every
`c_lambda(p)` is odd, and the squareclass contribution is `pi_p` for
positive `c_lambda(p)` and `bar(pi_p)` for negative `c_lambda(p)`.

## Fixed representatives: an exact source-exponent criterion

All 35 elements `sqrt(A_lambda)` belong to one quadratic extension of
`Q(i)` if and only if, for each prime dividing `t`, the 35 integers
`c_lambda(p)` all have the same sign. The sign may differ between primes.

Indeed, unique Gaussian prime valuations make the squareclasses equal
exactly under this condition. None is the trivial squareclass, because
its norm has rational squareclass `t`. A single quadratic extension
contains square roots of several base-field elements precisely when
their nontrivial squareclasses agree.

There is a direct characterization in the source exponents. Fix one
odd-total prime and sort the allocations of rows other than zero as
`b_1<=...<=b_7`. The condition is exactly one of

```
a_0 > b_4+b_5+b_6+b_7-b_1-b_2-b_3,
a_0 < b_1+b_2+b_3+b_4-b_5-b_6-b_7.                    (1)
```

For a character whose three other positive entries form a set `I`,

```
c_lambda=a_0+2 sum_(i in I) b_i-sum_(i=1..7)b_i.
```

Its minimum and maximum occur at the three smallest and three largest
allocations, proving (1). In the first case row zero is the unique
largest allocation; in the second it is the unique smallest. This is
an exact inequality on grouped prime exponents, not an assumption about
independent threshold layers.

No additional uniform congruence modulo four follows from this sign
condition. More precisely, the 35 integers `c_lambda(p)` have one common
residue modulo four if and only if `b_1,...,b_7` have one common parity.
Exchanging one selected index in the three-set changes `c_lambda` by
`2(b_i-b_j)`, which proves both implications. The common-field condition
uses the Gaussian squareclass, and does not impose this stronger condition
on the exponents of the remaining Gaussian square factors.

## Independent character conjugations: the broader criterion

Let `f_p(lambda)=1[c_lambda(p)<0]` for `p|t`. Allow the representative of
each character independently to be `A_lambda` or `bar(A_lambda)`.
The quarter roots can then all be placed in one quadratic extension if
and only if

```
f_p(lambda)+f_q(lambda) is constant in lambda
for every p,q dividing t,                              (2)
```

where addition is modulo two. Equivalently the sign functions
`sign(c_lambda(p))` agree up to a prime-dependent global sign.

To prove this, write the Gaussian squareclass as

```
[A_lambda]=[product_(p|t) pi_p]+sum_(p|t) f_p(lambda)[p].
```

Conjugating one character adds `[t]=sum_(p|t)[p]`, flipping every one
of its orientation bits simultaneously. Thus a choice of conjugation
bits `tau(lambda)` makes all classes equal exactly when every function
`f_p+tau` is constant. This is equivalent to (2).

In particular a single odd-total prime always permits this collapse
after independent conjugations, whether or not the fixed-representative
criterion (1) holds. If one prime has a constant orientation function
and another has a nonconstant one, even this broader collapse is
impossible. Reorienting a Gaussian prime flips its entire sign function
and preserves both criteria. Changing the convention for character
representatives need not preserve (1); criterion (2) explicitly allows
that freedom.

## Actual primitive data with collapse and unbounded pair-layer mass

Field collapse alone does not bound the source deletion content `F` of
[the weaker finite tests](weaker_intrinsic_finite_tests.md). Fix the
primary Gaussian prime `pi=-1+2i`, of norm five. For every integer
`e>=16`, use the eight allocations

```
(a_0,...,a_7)=(e,0,1,2,3,4,5,e-10),
z_i=pi^a_i bar(pi)^(e-a_i).
```

These are eight distinct Gaussian integers of common norm `N=5^e`,
with one literal canonical unit. Their common Gaussian gcd is a unit,
because the minimum allocation is zero and the maximum is `e`.
Their aggregate allocation is `2e+5`, so `t=5`. For the prescribed
35 characters, the minimum exponent `c_lambda` is exactly one.
Thus every quarter root lies in the same fixed extension
`Q(i,sqrt(pi))`. Both residues one and three modulo four occur among
the exponents, so the stronger congruence assertion fails on actual data.

The exact intrinsic layer weights are

```
V_1=11log5,       V_2=(e-14)log5,
V_3=2log5,        V_4=log5.
```

Indeed the sorted allocations are `0,1,2,3,4,5,e-10,e`; counting their
successive gaps gives these values. In particular the twist is fixed,
the singleton and triple masses are bounded, and `log F=V_2` diverges.

This family is explicitly outside the required endpoint geometry:

```
D-2W=(2e+45)log5>0,
```

which contradicts the necessary short-arc determinant inequality at
`C=1/2`. It is not an endpoint counterexample or a proof against using
field collapse jointly with the endpoint hypothesis. Its narrow role is
to show why the common-field condition by itself supplies no pair-mass
bound, even for actual primitive equal-radius source points.

The [exact checker](check_quarter_root_field_allocation_criterion.py)
checks (1), the modulo-four assertion, the independent-conjugation
criterion for pairs of actual allocation profiles, and Gaussian arithmetic
in this family. No independent factor-phase realization or formal
weighted profile is used for these fixtures.
