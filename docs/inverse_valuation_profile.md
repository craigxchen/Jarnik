# Backward descent: valuation profiles and squarefree stability

These are prose proofs, not new Lean theorems. The uniform endpoint
bound remains open in this project. This note strengthens the inverse
argument in `inverse_prime_power_descent.md`.

Let an arc of length at most `C sqrt(R)` contain `M>=1` distinct
Gaussian integers of norm `N=R^2`, where `R>1` and `C>0`. Set

```
T = log R,   c = max(0,log C),   X = exp(2T/M+4c).
```

## 1. A divisor-retention inequality valid for every cardinality

If a Gaussian divisor of norm `D` divides `k>=1` of the points, then

```
k log D <= 2T+4kc <= M log X.                       (1)
```

For odd `k>=3`, put `s=(k-1)/2` in the Ramana inequality after
division. It gives

```
(k+1) log D <= 2T+4k log C.
```

Since `D>=1`, this implies (1). For even `k>=4`, use `k-1` points
to obtain `k log D<=2T+(4k-4)log C`, which again implies (1).
For `k=2`, equal-norm Gaussian separation after division gives
`2<=C^2 R/D`, and therefore `2 log D<=2T+4log C-2log 2`.
For `k=1`, divisibility gives `D<=N`. These also imply (1).
An empty subfamily contributes zero and needs no geometric argument.

Thus a divisor retaining fraction `r=k/M` must obey
`r log D<=log X`. This avoids a rounding loss at small `k`.

## 2. The exact single-prime profile constant

For `e>=1`, define

```
n=e+1,  q=floor(n/2),  t=ceil(n/2),
m_e=q*t/(q+t)=floor((e+1)^2/4)/(e+1).
```

Every rational prime power `p^e || N` satisfies

```
m_e log p <= log X.                                 (2)
```

For a split prime `p=pi conjugate(pi)`, let `a` be the valuation
at `pi`, so `0<=a<=e`. Partition the points into `a<q` and `a>=q`,
of sizes `k_1,k_2`. The first part has common divisor
`conjugate(pi)^t`, the second `pi^q`. By (1),

```
k_1 t log p <= M log X,
k_2 q log p <= M log X.
```

Divide by `t` and `q` respectively and add, using `k_1+k_2=M`.
This proves (2). For inert primes, all points have a common divisor
of norm `p^e`; the same is true for the ramified prime. Apply (1)
with `k=M`, then use `m_e<=e`.

The constant is optimal for arguments using only the valuation
profile of this one prime. Indeed, for a probability distribution
on `a=0,...,e`, the interval `[l,u]` retains its interval mass and
has common divisor norm `p^(e-u+l)`. The preceding two intervals
force the maximum of

```
(e-u+l) * Pr(l<=a<=u)
```

to be at least `m_e`. The uniform distribution attains this value:
an interval of `h` cells gives `h(n-h)/n`, whose maximum is `m_e`.
This is a statement about the relaxation of valuation profiles;
it does not construct points on a short arc.

The first constants are

| exponent e | 1 | 2 | 3 | 4 | 5 |
| --- | --- | --- | --- | --- | --- |
| m_e | 1/2 | 2/3 | 1 | 6/5 | 3/2 |

In particular, `m_e>=e/4`, and `m_e>=2/3` for `e>=2`.

## 3. Repeated prime powers contribute little to a large cluster

Let `N_rep` be the product of the full prime powers `p^e || N`
with `e>=2`. From (2), every such prime satisfies

```
p <= X^(3/2),    e log p <= 4 log X.
```

Counting primes by counting integers gives the explicit bound

```
log N_rep <= 4 X^(3/2) log X.                        (3)
```

No prime number theorem is needed. Consequently, if fixed `K>3`
and `M>=K T/log T` along a sequence with `R` tending to infinity,

```
log N_rep = O_{C,K}(T^(3/K) log T) = o(T).           (4)
```

In particular, if `N` is squarefull (every prime divisor has
exponent at least two), then for every fixed `K>3`, eventually

```
M < K log R/log log R.                              (5)
```

More generally, if all prime exponents are at least `h>=1`, the
same argument, using the increasing sequence `m_e`, gives eventual
leading constant `2/m_h`. This concerns exponents in the rational
norm; it does not assume the Gaussian points themselves are powers.

## 4. Stability at the leading constant four

Suppose a sequence of endpoint clusters with fixed `C` satisfies

```
M = (4+o(1)) T/log T,    T -> infinity.
```

Then `X=T^(1/2+o(1))`. Equation (3) gives

```
log N_rep <= T^(3/4+o(1)) = o(T).                   (6)
```

Thus almost all of `log N=2T` comes from exponent-one primes.
Apart from the prime 2, these must split. They obey `p<=X^2`,
so `p<=T^(1+o(1))`. For any fixed `delta>0`, the exponent-one
primes below `T^(1-delta)` contribute at most
`psi(T^(1-delta))<= (log 4+4)T^(1-delta)=o(T)`.

There is also a precise balance conclusion. For an exponent-one
split prime, let `rho_p` be the larger of the two fractions of
points divisible by `pi` and by `conjugate(pi)`. Equation (1) gives

```
rho_p log p <= log X.                               (7)
```

For fixed `epsilon>0`, primes with `rho_p>=1/2+epsilon` therefore
satisfy `p<=T^(1/(1+2epsilon)+o(1))`. Their total exponent-one
logarithmic weight is `o(T)`, again by the Chebyshev bound.

Hence a hypothetical sequence saturating constant four must have,
up to `o(log N)` weight, distinct split primes of size
`(log R)^(1+o(1))`, each allocated almost equally between the two
orientations. The statements about size and balance mean that the
exceptional weight is `o(T)` for each fixed tolerance. No assertion
of existence of such lattice clusters is being made.

## 5. A further inverse conclusion: pairwise independent prime cuts

The same divisor argument also controls the joint allocations.
Take two distinct exponent-one split primes `p,q`. For any choice
of one Gaussian factor above each prime, let `r` be the fraction
of points divisible by their product. Then (1) gives

```
r (log p+log q) <= log X.                           (8)
```

In the near-four regime, if `log p,log q=(1+o(1))log T`, each
of the four joint cells has mass at most `1/4+o(1)`. Since their
masses sum to one, each also has mass at least `1/4-o(1)`.
This is stronger than merely balancing each prime separately:
the two binary allocation variables must be asymptotically
independent. In sign notation, their normalized inner product
is `o(1)`.

One can discard `o(T)` logarithmic weight and select a set of
prime columns for which these size and pairwise conclusions hold
uniformly. This follows by choosing a tolerance tending to zero
slowly in the fixed-tolerance conclusions of Section 4. Moreover,
the number of retained columns is

```
(2+o(1)) T/log T = (1/2+o(1)) M.                    (9)
```

Indeed their total logarithmic weight is `(2-o(1))T`, and every
retained prime has logarithm `(1+o(1))log T`.

For any fixed number `b` of distinct retained primes, each joint
allocation cell similarly has mass at most `1/(2b)+o(1)`. For
`b>=3` this exceeds the independent mass `2^(-b)`, so this argument
does not establish joint independence of higher order.

Equations (8)--(9) identify a more specific surviving pattern:
about `M/2` almost orthogonal, balanced binary prime columns on
`M` points. Dimension counting alone does not exclude it. Nor
does pairwise independence imply the required Gaussian phase
alignment. This conclusion is an inverse theorem for saturation
of constant four, not a classification of all large clusters.

## 6. What this says about the square-root route

For squarefree `N`, every Gaussian prime valuation of a point is
zero or one. Two points on this circle with matching valuation
parities are therefore associates. Consequently a square class
modulo units contains only associates; taking square roots cannot
combine distinct nonassociate points in this regime.

This is a substantive limit of the proposed parity argument:
near saturation of the current best leading constant forces a
mostly squarefree norm, rather than forcing many even valuations.
The balance conclusion also agrees with the balanced-cut
obstructions already studied. It does not control the joint
allocation of different primes or force their actual Gaussian
arguments to lie in the short-arc configuration.

The remaining opportunity is to use those arguments and exact
Gaussian integrality jointly. Equations (1)--(9) alone still leave
the uniform bound unproved. Also, a hypothetical unbounded family
with `M=o(log R/log log R)` need not satisfy the stability hypotheses;
excluding only the near-four regime would not prove uniformity.
