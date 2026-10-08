# An isolable Gaussian sign-pattern block is too light on a large endpoint arc

## Result and scope

This note combines the existing obtuse Bessel mass bound with an **exact
rational phase-isolation certificate**. It treats arbitrary moving split
primes, nested prime powers, common Gaussian content, and either literal
common or arbitrary row units. For a fixed row count `M>=9`, an endpoint
cluster cannot have a nonunit whole sign-pattern block whose coordinate
is in the rational span of its row differences at unbounded radius.
The coefficient and Roth constants depend on `M`; the statement is a
profile restriction, not a uniform point-count bound.
The earlier [fixed-profile phase theorem](fixed_hadamard_roth_phase_obstruction.md)
uses five separately isolable nonunit blocks against the total norm
budget. Here endpoint obtuseness makes each **whole sign-pattern** block
carry at most `2W/M`, so one isolable nonunit block suffices when `M>=9`.
The subsequent [integer pattern-lattice theorem](corank_one_local_bessel_pair_filter.md)
extends the standard-coordinate test to every nonzero integral
row-difference vector with `8||v||_1<M`, using all-order allocation
rigidity to rule out complete cancellation.

Let distinct Gaussian integers `z_1,...,z_M` have common modulus `R>1`
and lie on an arc of length at most `C sqrt(R)`. Use either

```text
literal common row unit,  0<C<=2;  or
arbitrary row units,      0<C<=sqrt(2).               (1)
```

Put `N=R^2`, `W=log N`. Factor the complete common Gaussian gcd as `d`
but retain it in the source rows. At each split prime `p`, after
subtracting the common prime valuations, use residual allocation
exponents `0<=a_i(p)<=e_p`, with `min_i a_i=0`, `max_i a_i=e_p`.
Write `D=log Norm(d)`, `W_var=W-D=sum_p e_p log p`.
For `M>=9`, `W_var>0`: without a varying split prime there are at most
four unit multiples of one common Gaussian integer.

For every threshold layer `(p,l)`, put

```text
s_(p,l)(i)=2 1_(a_i(p)>=l)-1 in {+1,-1}.          (2)
```

Identify a sign pattern with its negative by choosing one orientation
`s_S` for each unordered nonconstant cut class `{S,S^c}`. Let `V_S` be
the total `log p` weight of all layers in this class. Orient the Gaussian
prime at a layer so that its signed row pattern is `s_S`, and define
`Gamma_S` as the product of those oriented Gaussian primes, with layer
multiplicity. Threshold cuts at one prime form a chain, so it cannot
contain both `S` and `S^c`; thus `Gamma_S` is conjugate-primitive and

```text
Norm(Gamma_S)=exp(V_S).                            (3)
```

Let `T` be the finite set of sign classes actually present and let
`A=(s_S(i))_(i,S in T)` be the row-sign matrix. Its rational
row-difference space is

```text
V_Q={lambda^t A : lambda in Q^M, sum_i lambda_i=0}
       subset Q^T.                                  (4)
```

**Isolable-pattern theorem.** If the standard coordinate `e_S` lies in
`V_Q` and `V_S>0`, then for fixed `M>=9` and `C` in (1), `R` is bounded
by a constant depending only on `M,C`. In particular a sequence of
such clusters with unbounded radius and fixed `M>=9` has no present,
rationally isolable sign-pattern class. Full fair profiles with many
columns need not have any such coordinate, so this does not settle the
general endpoint problem.

## 1. Whole-pattern mass from obtuse Bessel

The exact pair norm distance is

```text
d_ij=log P_ij=sum_p |a_i(p)-a_j(p)| log p.
```

The unit-dependent chord bound in
[linear_allocation_affine_rigidity.md](linear_allocation_affine_rigidity.md)
gives `d_ij>W/2` for every pair under (1). On the layer Hilbert space
with weights `log p/W_var`, the row sign functions `f_i=s_(p,l)(i)` have
unit norm and

```text
<f_i,f_j>=1-2d_ij/W_var<0       (i!=j).          (5)
```

The [obtuse Bessel inequality](uniform_profile_extraction.md) says
`sum_i <f_i,h>^2<=2||h||^2` for every layer function `h`. Take `h` to
be `+1` or `-1` on the entire class `{S,S^c}` to align every layer
with `s_S`, and zero elsewhere. Then

```text
||h||^2=V_S/W_var,
<f_i,h>=(V_S/W_var)s_S(i),
M(V_S/W_var)^2<=2 V_S/W_var.
```

For `V_S>0`, this proves the support-uniform bound

```text
V_S<=2W_var/M<=2W/M.                            (6)
```

For any `q` distinct present classes, the same argument, or the trace
of the Bessel operator on their orthogonal layer indicators, gives
`sum_(j=1)^q V_(S_j)<=2qW_var/M`. This spectral mass fact is an
application of the existing Bessel lemma, not a new point-count bound.
It needs neither globally balanced cuts nor orthogonal sign columns.

## 2. Exact row product isolation, including nested layers and units

Since `e_S in V_Q`, clear denominators and multiply by two to obtain
integers `lambda_i` and a positive even integer `h` with

```text
sum_i lambda_i=0,
sum_i lambda_i s_T(i)=h 1_(T=S)      for every T in T. (7)
```

The source rows have the literal form

```text
z_i=epsilon_i d product_p pi_p^a_i(p) bar(pi_p)^(e_p-a_i(p)),
epsilon_i in {1,-1,i,-i}.                         (8)
```

Their signed product `Q=product_i z_i^lambda_i` is interpreted in
`Q(i)^*`. The common factor cancels because `sum lambda_i=0`.
At a prime,

```text
sum_i lambda_i a_i(p)
 =sum_l sum_(i in S_(p,l)) lambda_i
 =(1/2)sum_l lambda dot s_(p,l).                (9)
```

Equation (7) makes every nontarget layer vanish in (9), including
arbitrary nested powers at one regular prime. A target layer contributes
`+h/2` or `-h/2` according to its chosen orientation. Therefore

```text
Q=u (Gamma_S/bar(Gamma_S))^(h/2),
u=product_i epsilon_i^lambda_i in {1,-1,i,-i}. (10)
```

This is an exact Gaussian identity. It keeps the actual unit factor and
does not assert that prime layers with different cuts are independent.

Lift the source arguments `theta_i` along their containing arc of
angular width `Delta<=C/sqrt(R)`. Since `sum lambda_i=0`, their signed
sum obeys

```text
|sum_i lambda_i theta_i|<=(||lambda||_1/2) Delta.
```

Let `phi=arg(Gamma_S)`. From (10), for some integer branch `m`,

```text
dist(phi,T_(lambda,u))<=kappa Delta,
kappa=||lambda||_1/(2h),
T_(lambda,u)={(2pi m-arg u)/h : m in Z} mod 2pi. (11)
```

The target set is finite and consists of rational multiples of `pi`.
For a literal common unit, `u=1`; for arbitrary units there are only
four target shifts. The isolation cost `kappa` is the exact normalized
coefficient height of the row certificate.

## 3. Fixed-row-count phase-height gap

The fixed-target phase lemma in
[fixed_hadamard_roth_phase_obstruction.md](fixed_hadamard_roth_phase_obstruction.md)
applies to the conjugate-primitive nonunit `Gamma_S`. For any
`epsilon>0`, it gives

```text
dist(phi,T_(lambda,u))
 >=c_(M,lambda,h,epsilon) exp[-(1+epsilon/2)V_S]. (12)
```

Combining (11)--(12) with `Delta<=C/sqrt(R)` yields

```text
V_S>=W/(4+2epsilon)-K_(M,lambda,h,C,epsilon).  (13)
```

The upper bound (6) contradicts (13) at large `W` when
`1/(4+2epsilon)>2/M`. For every integer `M>=9`, choose
`0<epsilon<(M-8)/4`. This proves bounded `R` for a **fixed** sign
matrix and certificate.

There are only `2^(M-1)-1` unoriented nonconstant sign classes for
fixed `M`, hence only finitely many possible present-column subsets
and row-sign matrices. For each matrix and isolable coordinate choose
one integer certificate (7). Taking the largest of their finite Roth
and coefficient constants proves a bound depending only on `M,C`,
regardless of which rational primes occur or how their powers nest.
This finiteness is the precise fixed-`M` scope; it gives no constant
uniform as `M` increases.

There is a structural corollary. The affine allocation rigidity theorem
in [linear_allocation_affine_rigidity.md](linear_allocation_affine_rigidity.md)
gives row-difference rank `M-1` for the split allocation matrix under
(1). Each allocation column is a sum of its threshold indicators, hence
modulo a constant column lies in the span of the sign columns (2).
Therefore `A` also has full row-difference rank `M-1`. If exactly
`M-1` unoriented sign classes are present, then `V_Q=Q^T` and every
class is isolable. Consequently a fixed-`M>=9` sequence with unbounded
endpoint radius must have at least `M` present sign-pattern classes.
This is a separate restriction from the varying-prime count: one prime
with nested layers may contribute several classes, while many primes
may contribute the same class. The class count can still proliferate
with `M` and `R`.

## 4. A sufficient q-column pairing criterion and its conditioning cost

Suppose `q` present exceptional sign columns `e_1,...,e_q` are
balanced (`e_j dot 1=0`) and linearly independent, and every other
present sign column is orthogonal to their span `E`. Let

```text
H=(e_a dot e_b)_(a,b<=q).                         (14)
```

The matrix `H` is positive definite and integral; it need not be
diagonal. Taking the signed row products with exponents `e_a(i)`
cancels common content and every regular layer exactly. The lifted
phase equations have the form

```text
H phi=-xi+2pi n+delta,
|delta_a|<=M Delta/2,
xi_a in (pi/2)Z,             n in Z^q.          (15)
```

Inverting `H` gives each exceptional phase a rational-`pi` target
with the explicit error multiplier

```text
kappa_j=(M/2)sum_a |(H^(-1))_(j,a)|.             (16)
```

For orthogonal columns `H=M I`, this is `1/2`. Near-dependent columns
can make it much larger, but for fixed `M` it changes only the
constant in the Roth bound, not the exponent. The Bessel mass estimate
`sum_j V_j<=2qW_var/M` does not depend on this conditioning.
Balance is a convenient special case, not a necessity. For arbitrary
exceptional sign columns let `P=I-11^t/M`, `c_a=P e_a`, and suppose the
centered columns are independent while every regular sign column is
orthogonal to all `c_a`. Then the rational pairing

```text
K=(c_a dot e_b)=(e_a dot P e_b)_(a,b<=q)
```

is positive definite. The rational row coefficient
`lambda_j=sum_a (K^(-1))_(j,a)c_a` has zero sum and isolates the
`j`-th exceptional standard coordinate exactly; multiplying it by a
common denominator makes the signed Gaussian row product literal.
The corresponding phase error multiplier is

```text
kappa_j=(1/2)sum_a |(K^(-1))_(j,a)| ||c_a||_1.   (17)
```

This recovers (16) when each `e_a` is balanced. The centered criterion
gives an explicit sufficient test for unbalanced exceptional cuts,
while the direct membership condition `e_S in V_Q` remains the exact
criterion for one whole class.
If `H` is singular, the individual phases are not all isolated by
(15); one must seek coupled Gaussian products and charge their exact
prime cancellations. The fixed coupled-product criterion in
[coupled_block_roth_packing.md](coupled_block_roth_packing.md) addresses
that separate situation. For example, the four distinct balanced
eight-row sign columns

```text
e_1=(+,+,+,+,-,-,-,-),
e_2=(+,+,-,-,+,+,-,-),
e_3=(+,+,+,-,-,+,-,-),
e_4=(+,+,-,+,+,-,-,-)
```

obey `e_1+e_2=e_3+e_4`. A formal variation of their block phases along
`(1,1,-1,-1)` leaves every row phase unchanged, so no individual phase
is isolated by the pairing matrix. The complete fair sign family has no
orthogonal exceptional subspace of this kind.

## Verification and limitation

Run `python3 docs/check_isolable_pattern_bessel_phase_gap.py`. The
standard-library checker verifies the exact signed-product identity
with nested cuts at one prime, the `q=2,3,4` pairing inverses and
conditioning multipliers, one centered unbalanced `q=2` pairing, and a
rank-deficient four-column witness.
The Bessel and Roth inequalities are proof inputs, not finite-search
claims. Approximate row-sign orthogonality or correction factors in a
centrally extracted fair profile are outside the theorem: their
uncontrolled phases cannot simply be absorbed into a fixed target.
