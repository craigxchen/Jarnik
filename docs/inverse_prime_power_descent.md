# Working backward: prime-power descent and a sharper sublogarithmic bound

This note revisits the existing sublogarithmic proof. It gives an
unconditional argument with leading constant `4+epsilon`, improving
the manuscript's `8+epsilon` constant. It does not prove uniform
boundedness. The new prime-power step and its assembly below are
prose proofs; the complete new bound is not yet formalized in Lean.

Write `N=R^2>1`, and put `A=log 4+4`. All logarithms are natural.

## 1. Backward implication of a large cluster

Let an endpoint arc contain more than `4s` lattice points, with
integer `s>=1`. If a nonzero Gaussian integer `d` divides at least
`2s+1` of them, put `D=|d|^2`. Divide those points by `d`. They have
radius `R/sqrt(D)` and squared chord diameter at most `C^2 R/D`.
The existing Ramana subcritical bound therefore requires

```text
(R^2/D)^(s^2) <= (C^2 R/D)^(s(2s+1)).
```

Taking logarithms and simplifying gives the necessary condition

```text
(s+1) log D <= log R + (4s+2) log C.              (1)
```

Equivalently, every Gaussian divisor shared by `2s+1` selected
points must have norm at most

```text
H_s = exp((log R+(4s+2)log C)/(s+1)).
```

The real-valued chord inequality suffices; discarding its natural
floor only weakens this necessary condition.

## 2. A large prime power always supplies a suitable divisor

Factor `N=product_p p^(e_p)` and define

```text
Y = max_(p divides N) p^ceil(e_p/2).
```

For a prime attaining this maximum there is a Gaussian divisor of
norm at least `Y` which divides at least half of all the points:

* If `p=1 mod 4`, write `p=pi conjugate(pi)` and `e=e_p`.
  Unique factorization gives
  `v_pi(z)+v_conjugate(pi)(z)=e` for every point. At least one
  valuation is at least `ceil(e/2)`. Thus one of the two fixed
  divisors `pi^ceil(e/2)`, `conjugate(pi)^ceil(e/2)` divides at
  least half the points. Its norm is `Y`.
* If `p=3 mod 4`, representability forces `e` even, and every
  point is divisible by the rational Gaussian integer `p^(e/2)`.
  Its Gaussian norm is `p^e>=Y`.
* If `p=2`, every point is divisible by `(1+i)^e`, whose norm
  is `2^e>=Y`.

Crucially, the split case chooses the whole prime power at once.
It pays a factor two in cardinality, not a factor two for each
successive division by `pi`.

Consequently, more than `4s` points implies

```text
(s+1) log Y <= log R+(4s+2)log C.                 (2)
```

## 3. Elementary lower bound for the largest half-exponent prime power

Set `L=lcm(1,...,floor(Y))`. For each prime dividing `N`, the
definition of `Y` implies `v_p(L)>=ceil(e_p/2)`. Hence `N` divides
`L^2`. In terms of Chebyshev's second function,

```text
log R = (log N)/2 <= log L = psi(Y) <= A Y.
```

Thus

```text
Y >= log R/A.                                    (3)
```

The bound `psi(x)<=(log 4+4)x` for `x>=0` is already available as
`Chebyshev.psi_le_const_mul_self` in Mathlib and is used in the
project's `MertensFirst.lean`. No prime number theorem or assertion
about prime gaps is needed here.

## 4. Finite bound and asymptotic constant

Assume `R>max(1,C^2)` and

```text
d_R = log(log R)-log A-4 log C > 0.
```

Choose

```text
s = max(1, floor((log R-2 log C)/d_R)).
```

Since `(s+1)d_R>log R-2 log C`, equation (3) gives

```text
(s+1)log Y-(4s+2)log C
 >= (s+1)d_R+2 log C > log R.
```

This contradicts (2) if the arc contains more than `4s` points.
Therefore the explicit bound is

```text
M <= 4 max(1, floor((log R-2 log C) /
                    (log(log R)-log(log 4+4)-4 log C))).    (4)
```

For each fixed `C>0`, the right side is
`(4+o_C(1)) log R/log(log R)`. In particular, for every `K>4`,
there is a threshold `R_0(C,K)` such that every endpoint arc at
radius `R>=R_0(C,K)` contains at most

```text
K log R/log(log R)
```

lattice points, uniformly in arc position. The cases `N<=1` do
not affect this asymptotic assertion.

## 5. What the inverse argument says, and does not say

This changes the starting point for a structural inverse theorem.
Clusters near the old constant eight are already excluded by the
prime-power branch. A genuinely large surviving cluster must obey
the shared-divisor bound (1) for every large subfamily, and its
largest half-exponent prime power is constrained by (2)--(3).

For fixed `s`, however, the upper bound `H_s` in (1) still grows
like `R^(1/(s+1))`. This greatly exceeds the logarithmic lower
bound in (3). Thus these necessary conditions alone do not rule
out arbitrarily large fixed cardinalities as the radius grows.
No bounded factor-pattern classification has been deduced from them.

The complete new theorem remains a prose derivation. Existing Lean
results verify the Ramana obstruction, prime descent for exponent
one, and the stated Chebyshev estimate. They do not yet verify
the simultaneous prime-power pigeonhole or the assembly in (4).
The manuscript's existing theorem remains valid and is not being
silently replaced by an unformalized statement.

## 6. A structural obstruction to test next: common square class

Suppose every point in a cluster can be written

```text
z_j = g epsilon_j w_j^2,
```

where `g` is one nonzero Gaussian integer, each `w_j` is Gaussian,
and `epsilon_j` is a Gaussian unit. There are only two unit classes
modulo squares: absorb a factor `-1` into `w_j^2` to arrange
`epsilon_j` in `{1,i}`.

For each of these two classes, choose signs of the square roots so
that their arguments lie in the half-sized lifted argument interval.
Their radius is `sqrt(R/|g|)`. An original arc of length at most
`C sqrt(R)` therefore produces a root arc of length at most

```text
C/(2 sqrt(|g|)) <= C/2.
```

Distinct original points in the same class give distinct selected
roots. Equal-norm Gaussian roots have distance at least `sqrt(2)`.
Ordering those roots along their arc yields the uniform bound

```text
M <= 2 (floor(C/(2 sqrt(2)))+1).                   (5)
```

In particular, such a configuration has at most two points at `C=1/2`.
This argument is finite and does not use an asymptotic template.

For example, if every Gaussian prime valuation of every point is
even, all points have this form with `g=1`. The same holds after
removing a common Gaussian factor whenever all remaining valuations
are even. This rules out an apparent model in which the norm has
even exponents and every split-prime allocation uses only the two
endpoints `0,e`, with all these `e` even.

The continuation in [the valuation-profile note](inverse_valuation_profile.md)
now proves a stability result in the opposite direction near constant
four: most norm weight must come from exponent-one split primes, with
almost balanced, pairwise almost independent allocations. It also
improves the constant to three when the norm is squarefull.

What has **not** been proved is that a large endpoint cluster must
have one common square class, or even boundedly many square classes.
Different rows can have odd valuations at different primes. The next
inverse question is whether resisting the stronger prime-power
descent forces enough agreement of these valuation parities to apply
(5). Cardinality-based descent estimates alone have not supplied
that agreement.
