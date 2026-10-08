# Compatible Paley phases survive every scalar Gaussian coordinate gap

The canonical five-copy Paley profile admits distinct endpoint-scale row
phases for which **every** nonzero saturated integer character satisfies
both Gaussian coordinate lower bounds. The factors can have the squared
moduli of actual distinct split primes, with total logarithmic squared
modulus `W=20M log M+O_C(M)`. The factors are complex numbers, not Gaussian
integers. Thus phase coherence together with all these scalar coordinate
gaps does not settle the missing joint Gaussian-integrality condition.

This improves the scale and the sign-profile scope of the earlier
[compatible-phase relaxation](integer_combination_separation.md), which
uses all balanced half-cuts and a sufficient `W` of order `M^2`. It also
enforces the real coordinate gap, hence the imaginary gap after every
Gaussian unit rotation. It does not construct lattice points on an
endpoint arc.

## Precise statement

Let `M=q+1>=12`, with `q` a prime congruent to three modulo four, and let
`S` be the canonical five-copy, one-flip sign matrix in the
[all-character theorem](paley_canonical_all_character_half_height.md).
There are `r=5(M-1)` physical columns. Fix `0<C<=1`, and suppose

```text
tau <= w_j < tau+1/r,
tau >= 4 log M + kappa,       kappa=4 log(1000/C),
W=sum_j w_j,                 Delta=C exp(-W/4).
```

There exist complex factors `g_j` with squared modulus `exp(w_j)` such
that the `M` row products

```text
Z_i=product_j g_j^((1+S_ij)/2) conjugate(g_j)^((1-S_ij)/2)
```

have distinct arguments in an interval of width `Delta`. Their common
modulus is `R=exp(W/2)`, so they occupy an arc of length at most
`C sqrt(R)`. Moreover, for every nonzero integer `c` with `sum_i c_i=0`,
the integers `v_j=(c^T S)_j/2` define the formal reduced character

```text
G_c=product_(v_j>0) g_j^v_j product_(v_j<0) conjugate(g_j)^(-v_j),
H(c)=sum_j |v_j| w_j,
min(|Re G_c|,|Im G_c|) >= 1.                            (1)
```

Equivalently, `|Im(u G_c)|>=1` for every Gaussian unit `u`.
The word reduced refers to cancelling opposite formal orientations of
the same physical factor. No Gaussian factorization of `g_j` is asserted.
All integer amplitudes are covered, with no coefficient cutoff.

For each fixed `C`, the weights can be `w_j=log p_j` for actual distinct
primes `p_j=1 mod 4`, for every sufficiently large Paley order. Indeed,
apply the [dyadic prime-cluster argument](strict_obtuse_prime_box_countermodels.md)
with `X=(1000M/C)^4`, partitioning `[X,2X]` into `r` equal intervals.
One interval eventually contains `r` split primes. Its left endpoint
`A` gives `tau=log A`, the displayed width condition, and

```text
W=4r log M+O_C(r)=20M log M+O_C(M).                     (2)
```

For `C>1`, the construction with `C=1` already supplies the stated
existence within the larger allowed arc. The theorem conditional on the
weight bounds holds for every `M>=12`; the prime-cluster existence is
asymptotic.

## The required canonical support bounds

Write `B(c)=sum_j |v_j|`. The sharp canonical theorem gives the following
more detailed bounds, which are essential for counting the characters.

* A ternary vector of support two is `e_i-e_k`. If its pair contains
  infinity, `B>=5M/2-2`; there are at most `2M` such ordered vectors.
  Every other pair has `B>=5M/2`, because among finite rows at most one
  supported diagonal flip is beneficial.
* A zero-sum ternary vector of support `h>=4` satisfies
  `B>=5M/2+h-4`. This follows from equation (5) of the canonical proof,
  retaining its `h` term. There are at most `(2M)^h` such vectors.
* An integer vector of amplitude `A=max_i|c_i|>=2` satisfies
  `B>=3MA/2`. There are at most `(2A+1)^M` vectors of amplitude `A`;
  this overcount may ignore the zero-sum condition.

In particular `D_c=B(c)/2-r/4` is positive in every case. The respective
lower bounds are

```text
1/4;       5/4;       (2h-3)/4;       [(3A-5)M+5]/4.      (3)
```

## One probability estimate for all characters

Choose independent uniform `phi_i` in `[0,Delta]`. For fixed nonzero `c`,
condition on all but a coordinate with `c_i!=0`. The variable
`c dot phi` ranges uniformly over an interval of length
`|c_i| Delta`. In each period of length `2pi`, the set where
`|sin(t/2)|<epsilon` has length `4 arcsin(epsilon)<=2pi epsilon`.
The partial periods give

```text
Pr[ |sin(c dot phi/2)|<epsilon ]
 <= epsilon*(1+4pi/(|c_i|Delta)) <= 14 epsilon/Delta.
```

The same bound holds for cosine by translation. Apply their union bound
with `epsilon=exp(-H(c)/2)`. Since `H(c)>=tau B(c)`, `W<r tau+1`,
and `exp(1/4)<2`, the probability that either coordinate condition fails
is at most

```text
(28/C) exp(W/4-H(c)/2)
 <= (56/C) exp(-tau D_c).                               (4)
```

Extract `exp(-kappa/4)` from every term in (4), using `D_c>=1/4`, and
use the remaining `tau>=4 log M` in the four classes (3). Their total
normalized contribution is at most

```text
Q_M=2 + M^(-3) + 16/(M-2)
       + M^(-5) (5/M)^M / [1-(5/M^3)^M].                (5)
```

For the third term, overcount all integer `h>=4`, rather than just even
supports at most `M`:

```text
sum_(h>=4) (2M)^h M^(-(2h-3)) = 16/(M-2).
```

For the last term use `2A+1<=5^(A-1)` for every integer `A>=2`; the
resulting geometric series is exactly the last term in (5). This handles
unbounded amplitudes. For `M>=12`, `Q_M<4`: the first three terms are at
most `2+1/1728+8/5`, and the last term is less than `2/12^5`.
Consequently the union of all countably many bad events has probability
less than

```text
(56/C) exp(-kappa/4) Q_M < 224/1000 <1.                  (6)
```

The angles are pairwise distinct with probability one. An angle tuple
outside the union therefore has the required distinctness and all the
simultaneous scalar gaps.

## Lifting the angles to compatible factor products

The physical matrix `S` has full row rank. For each row `i`, compare its
flipped column with an unflipped copy of the same original label. Their
difference is exactly `-2 H_(i,a_i) e_i`. There are unflipped copies for
every label, since at most two of the five copies were flipped. These
differences span every row-coordinate direction.

It follows that the chosen real angle tuple has a real lift
`S theta=phi`. Put `g_j=exp(w_j/2+i theta_j)`. Then the row products have
exact arguments `phi_i` modulo `2pi`, and every character has

```text
|G_c|=exp(H(c)/2),        arg G_c=(c dot phi)/2 mod 2pi.
```

Equation (1) follows from the chosen angle tuple. Also `v!=0` whenever
`c!=0`, by full row rank. The integrality of `v` follows from
`c^T S_j=sum_i c_i-2 sum_(S_ij=-1)c_i`.

The same flipped/unflipped comparisons show that the saturated
rational-row lattice is precisely `{c^T S/2: c integer, sum c=0}`.
Indeed integrality of a rational character forces twice every row
coefficient integral. Thus the construction does not omit a finer
rational-row character lattice.

## Exact limitation

For genuine disjoint split-prime Gaussian factors, a nonzero reduced
character is conjugate-primitive and nonunit, so neither Gaussian
coordinate vanishes and both have absolute value at least one. The
construction satisfies all those lower bounds and the exact phase
identities, but does not make the coordinates integers. It likewise does
not make the original factors or row products Gaussian integers.

This distinction is visible in an explicit choice of the lift. Use the
flipped physical column and a third, unflipped copy at each label to
represent each `phi_i e_i`. At least one unused physical column has
`theta_j=0`, so its factor is the positive real number `sqrt(p_j)` when
actual prime weights are used. That factor is not a Gaussian integer.
No assertion excludes other lifts; the theorem is a relaxation showing
which inputs alone fail, not a proof of Gaussian feasibility.

The [exact checker](check_paley_compatible_phase_coordinate_gap.py)
checks all ternary support tiers at order twelve, samples arbitrary
integer amplitudes, verifies the rational sums and lift, and retains
the factor two required by the unit-invariant coordinate condition.
The countable probability argument, rather than finite character tests,
establishes (1) for every amplitude.
