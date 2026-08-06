# Positive subset transitions and the Bernoulli fixed-point obstruction

## Status

This note continues the prime-reveal/renormalization program using the ordered
positive angular gaps.  It proves a stronger arithmetic separation statement:
every nonempty subset of consecutive transitions produces a norm-one Gaussian
rational in the same endpoint arc, not only intervals of transitions.

The strengthening is genuine, but it still does not prove the uniform endpoint
bound.  There is an exact abstract Bernoulli transition model satisfying every
resulting conductor-width inequality for arbitrarily many transitions.  This
identifies the residual Gaussian numerator, rather than the exponent width, as
the indispensable next state variable.

---

## 1. Ordered transition quotients

Let

```text
z_0,z_1,...,z_m
```

be distinct Gaussian lattice points of common norm `N`, listed in increasing
angular order inside an arc of length

```text
L <= C N^(1/4).
```

Put

```text
u_i = z_(i+1)/z_i,
```

and let `g_i` be its positive angular increment.  Then

```text
0 < g_i,
sum_i g_i <= C N^(-1/4).
```

For every nonempty subset `J subset {0,...,m-1}`, define

```text
u_J = product_(i in J) u_i.
```

Positivity of the gaps gives a canonical lift

```text
0 < arg(u_J) = sum_(i in J) g_i <= C N^(-1/4).
```

In particular `u_J != 1` once `N` is large enough that the full arc angle is
less than `2 pi`.

---

## 2. Transition exponent vectors

In the squarefree split-prime model, write the orientation of point `z_i` at a
prime `p` as

```text
epsilon_p(i) in {-1,+1}.
```

The signed transition is

```text
s_p(i) = (epsilon_p(i+1)-epsilon_p(i))/2
       in {-1,0,+1}.
```

For a selected transition subset `J`, put

```text
c_p(J) = sum_(i in J) s_p(i).
```

Then

```text
u_J = product_p (pi_p/conjugate(pi_p))^(c_p(J)).
```

Define its conductor width

```text
W(J) = sum_p |c_p(J)| log p.
```

Clearing Gaussian denominators gives a nonzero Gaussian integer numerator, so

```text
|u_J-1| >= exp(-W(J)/2).
```

Combining this with the endpoint upper bound gives, for every nonempty `J`,

```text
W(J) >= (1/2) log N - 2 log C.                    (2.1)
```

After subdividing an unbounded cluster into shorter subarcs, the constant `C`
may tend to zero while the retained number of points still tends to infinity.
Then the right side of (2.1) has a positive additive excess tending to infinity.

This is stronger than the usual pairwise statement.  The selected set `J` need
not be an interval, and the resulting quotient need not be a quotient of two
original cluster points.

---

## 3. Alternating transition patterns

For one prime, the nonzero entries of

```text
s_p(0),...,s_p(m-1)
```

alternate in sign, because a two-state orientation cannot switch twice in the
same direction.

Thus a prime row is determined, up to one overall sign, by its transition
support

```text
T subset {0,...,m-1}.
```

Write `s_T` for the canonical alternating row on `T`.  If conductor height is
normalized to one, any abstract transition model is a probability measure `mu`
on such rows, and (2.1) becomes

```text
E_mu |sum_(i in J) s_T(i)| >= 1/2-o(1)
```

for every nonempty transition subset `J`.

---

## 4. Exact Bernoulli fixed point

There is an abstract model satisfying the critical inequalities for every
length `m`.

Choose the transition support `T` uniformly from all `2^m` subsets.  Equivalently,
fix the initial state and let all later states be independent unbiased signs.
For a fixed nonempty `J`, summation by parts writes

```text
2 sum_(i in J) s_T(i)
```

as a nonzero integer linear form in independent Rademacher signs.

For every nonzero integer coefficient vector `a`,

```text
E |sum_k a_k epsilon_k| >= 1.
```

Indeed the sum is integer-valued; with one nonzero coefficient this is
immediate, and with at least two nonzero coefficients the largest atom has
probability at most `1/2`, while every nonzero value has absolute value at least
one.  Consequently

```text
E_T |sum_(i in J) s_T(i)| >= 1/2.                  (4.1)
```

For interval sets `J`, equality holds: the sum telescopes to half the difference
of two independent signs.

Therefore the Bernoulli transition law is an exact critical fixed point for all
of the subset-width inequalities, for arbitrarily large `m`.

---

## 5. What this obstruction rules out

No proof can establish the endpoint bound using only:

1. the alternating sign pattern of each prime along the path;
2. the conductor width `W(J)` of every selected transition subset;
3. run counts or total variation of the prime orientations;
4. positivity and the finite total sum of the angular gaps; or
5. any scalar potential that depends only on those exponent widths.

The Bernoulli fixed point satisfies all of that abstract data.  It is the
ordered analogue of the balanced Plotkin/Hadamard obstruction encountered by
the determinant and pair-distance methods.

---

## 6. The arithmetic datum the fixed point cannot fake

For every nonempty `J`, choose a reduced Gaussian integer

```text
A_J = X_J+iY_J
```

representing `u_J=A_J/conjugate(A_J)`.  Then

```text
|u_J-1|^2 = 4 Y_J^2 / Norm(A_J).
```

The conductor width controls `Norm(A_J)`, but the positive integer

```text
R(J) = 4 Y_J^2
```

is additional arithmetic information.  The exact endpoint inequality is

```text
W(J) - (1/2) log N
  >= log R(J) - 2 log C.                            (6.1)
```

The abstract Bernoulli model silently sets every residual numerator to its
smallest possible value.  Actual Gaussian transition products must satisfy the
multiplicative cocycle identities among all `A_J`, so the integers `R(J)` cannot
be assigned independently.

This is the precise form of the heuristic that one pair may attain the
arithmetic separation lower bound, but a long configuration should not be able
to attain it for every relevant subset simultaneously.

---

## 7. Revised renormalized state

The prime-reveal state must retain, at each scale, both:

```text
exponent width W(J)
```

and

```text
residual numerator R(J).
```

A useful global theorem is a lower bound such as

```text
sum_(J in F) log R(J) >= c |F| log |F| - O(|F|)
```

for a canonically chosen family `F` of transition subsets, unless the transition
system has bounded rank or a common Gaussian divisor permitting genuine
descent.  The inert-prime projective-collision theorem in
`inert_prime_residual_growth.md` now proves this order of aggregate growth for
the forced-reduced pairwise Plucker residuals.

By itself, however, aggregate growth does not yet turn (6.1) into an absolute
bound.  The conductor-width surplus can grow at the same time.  The missing
inequality must compare the two quantities and charge repeated conductor
re-entry rather than simply lower-bound the residuals.

The finite four/five-transition test has also now been carried out.  The
residuals are linked by the positive half-angle continuant relations and by the
multiplicative identities

```text
u_(J union K) = u_J u_K
```

for disjoint subsets, but a fifth compatible layer is possible.

There is an exact four-transition/five-point example on the circle of norm

```text
N=27625
```

with points

```text
(120,115), (115,120), (101,132), (88,141), (83,144).
```

Its successive half-angle cotangents are

```text
47, 18, 21, 57,
```

all subset-product primitive ordinates are at most `11`, and the five points
lie in an arc shorter than `4 N^(1/4)`.

More decisively, the five cotangents

```text
157, 993, 307, 302, 278
```

give six exact points of norm

```text
N=6076533125
```

inside an arc shorter than `10 N^(1/4)`.  All `31` nonempty subset products
have primitive ordinate at most `191`.  Their transition norms repeatedly
reuse the primes `5` and `29`; they do not form five disjoint conductor petals.

Thus no theorem saying that a fifth transition is algebraically incompatible,
or that it forces one individually large residual, can be correct.  What must
eventually fail is residual efficiency relative to the conductor cost as the
number of transitions tends to infinity.

---

## 8. Verdict

The positive-subset strengthening is real, but exponent-width information still
has an exact arbitrary-length fixed point.  The next successful potential must
couple the residual Gaussian numerators to the prime re-entry surplus, not only
measure either one separately.

The active problem is now concrete:

```text
prove that the cocycle-linked residual growth cannot be paid for by repeated
reuse of one conductor core along an unbounded positive transition process.
```

The exact algebraic form of that coupling is developed in
`transition_content_cocycle.md`: on every interval, the weighted prime re-entry
count is precisely the rational coordinate content of the product of its
primitive transition blocks.
