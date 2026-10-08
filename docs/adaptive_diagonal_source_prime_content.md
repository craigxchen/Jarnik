# The source-prime part of the anchor-radius product

The [anchor-product relation](adaptive_diagonal_anchor_radius_product.md)
can be extended to the primes dividing the source radius. In its notation,
let `m>=2`, let the distinct rational source phases have least squared
radius N, and let

```text
A=product_(i<j)n_ij^2 / product_j N_j,
B=product_(i<k)b_ik,
U=a+d,                 V=a-d,
G=gcd(N,|UV|).
```

Here a,d are coprime positive integers, and `a!=d`. Then the full
integer divisibility, with no untreated source primes, is

```text
A divides B^(m-2) G^((m-1)(m-2)).                       (1)
```

In particular, whenever `G=1`,

```text
A divides B^(m-2),
product_(i<j)n_ij^2 / B^(m-2)
 <= product_j N_j <= product_(i<j)n_ij^2.                (2)
```

This applies uniformly to `diag(2,1)`, `diag(3,1)`, and the sequence

```text
(a,d)=(2^k+1,2^k-1),             k>=1.                 (3)
```

For the first map, `UV=3`; for the second, `UV=8`; for (3),
`UV=2^(k+2)`. A least primitive squared circle radius has only odd
split rational prime factors, so each is coprime to N. The last
sequence approaches a conformal map in real projective coordinates;
(1) does not assert that its target radii remain controlled.

The result supplies a complete arithmetic content bound for the
product of genuinely different source-anchor targets. It is not an
upper bound for one sufficiently small target, a map-existence
theorem, or an improved general endpoint count.

## 1. Exact local map, including the source valuations

Fix an odd split prime `p=pi bar(pi)`. Normalize one source phase to
one, and put

```text
s_i=v_pi(q_i),       S=max_i s_i-min_i s_i=v_p(N),
x_ij=q_i/q_j,        X_ij=s_i-s_j,
t_ij=v_pi((U x_ij+V)/(V x_ij+U)).
```

As in the preceding proof, `t_ji=-t_ij`, and

```text
v_p(A)=sum_j [sum_i |t_ij|-(max_i t_ij-min_i t_ij)].     (4)
```

If `s_i=s_k`, the relative source half-angle row is a p-adic unit.
The exact source residue identity is therefore

```text
v_pi(x_ij-x_kj)=X_ij+v_p(b_ik),                         (5)
```

where `X_ij=X_kj`. To see this, divide the left difference by
`x_kj`. Its remaining factor is `q_i/q_k-1`, whose valuation equals
that of the nonzero integral imaginary coordinate of the primitive
relative half-angle row. Its norm is the source two-point squared
radius, or twice that radius. Equal source levels make its norm a
p-unit. Thus (5) includes odd-odd half-angle rows and independent
Gaussian row units; p=2 is excluded here only for the factor two.

At an odd p, U and V cannot both have positive valuation, because
their common divisor divides `2gcd(a,d)=2`.

If both are units, `X_ij!=0` makes the two Möbius terms have the same
valuation, so `t_ij=0`. At a fixed anchor, the only nonzero t's
therefore have `s_i=s_j`. Their same-sign pair minima are bounded
by `v_p(b_ik)`, by subtracting the two numerators or denominators
as in the previous proof. Consequently

```text
v_p(A)<=(m-2)v_p(B)             when v_p(UV)=0.          (6)
```

This argument works even when p divides N: it requires equal levels
only within the nonzero target group, not that all source phases
be units.

## 2. Clipped source differences and the boundary excesses

Suppose now `h=v_p(UV)>0`, and define

```text
clip_h(X)=max(-h,min(X,h)).
```

Valuation of the two-term numerator and denominator gives the exact
following description. If `p|V`, then

```text
t_ij=clip_h(X_ij)                       if |X_ij|!=h,
t_ij=sign(X_ij)(h+e_ij), e_ij>=0         if |X_ij|=h.   (7)
```

If `p|U`, every displayed sign on the right is reversed. For example,
when `p|V` and `X_ij=h`, the numerator terms both have valuation h
and its denominator is a unit, giving `t_ij=h+e_ij`. At `X_ij=-h`,
the numerator has valuation -h and cancellation in the denominator
gives `t_ij=-h-e_ij`. Away from these two levels the smaller
valuation in each sum is unique. The case `p|U` exchanges these
roles. In particular no boundary cancellation can reverse the sign
or reduce the absolute value below h.

For two boundary excesses of the same sign at a fixed anchor j,
the corresponding source levels are equal. Equation (5) bounds
their minimum by the source residue:

```text
min(e_ij,e_kj)<=v_p(b_ik).                              (8)
```

Here is the valuation calculation in all four cases. Subtract the
terms in the canceling sum. The table lists their common baseline
valuation c and that of their difference:

| Divisible coefficient | Source level X | Canceling sum | c | Difference valuation |
| --- | --- | --- | --- | --- |
| V | h | `U x+V` | h | `h+v_p(b_ik)` |
| V | -h | `V x+U` | 0 | `v_p(b_ik)` |
| U | h | `V x+U` | h | `h+v_p(b_ik)` |
| U | -h | `U x+V` | 0 | `v_p(b_ik)` |

Since the two sum valuations are `c+e_ij` and `c+e_kj`, the
nonarchimedean triangle inequality gives (8). No condition on
`v_p(ad)` is needed.

## 3. The exact baseline loss and the global divisibility

For the clipped baseline define

```text
E_h(s)=2sum_(i<j) min(|s_i-s_j|,h)
       -sum_j [min(s_j-s_min,h)+min(s_max-s_j,h)].       (9)
```

This is nonnegative: each bracket removes the largest magnitude on
each side from its row's absolute sum. Adding boundary excesses
changes each one-sign row loss by the sum of those excesses minus
their largest one. Indeed, any boundary entry already attains the
baseline maximum h on that side. Choose an index of largest excess.
Each other excess is its minimum with that chosen excess, and (8)
charges it to a distinct source pair avoiding j. Summing over j,
each pair is used at most `m-2` times. This proves the sharper local
bound

```text
v_p(A)<=E_h(s)+(m-2)v_p(B).                             (10)
```

The baseline term has a particularly simple uniform estimate.
Choose one minimum-level index L and one maximum-level index H,
distinct if `S>0`. If `S=0`, it is zero. For `S>0`, monotonicity
of clipping shows that the two entries at L,H give the full range
in every baseline row. Therefore

```text
E_h(s)=sum_(i notin {L,H}) sum_j min(|s_i-s_j|,h)
      <=(m-2)(m-1)min(h,S).                             (11)
```

This identity also covers repeated minimum or maximum levels; the
zero distances cause no problem. Since `min(h,S)=v_p(G)`, (6),
(10), and (11) prove (1) prime by prime. Inert primes and two never
divide A, since they divide none of the least primitive squared
circle radii used to define it. Their presence in the right side
is harmless.

The baseline term is necessary in general. For source half-angle
rows `H=(1,2+i,3+i)` and `diag(3,2)`, one has

```text
N=25, B=1, U=5, V=1,
(n_01,n_02,n_12)=(5,85,445),
(N_0,N_1,N_2)=(425,445,7565), A=25.                    (12)
```

At p=5 the source levels, for a suitable orientation, are `0,1,-1`.
Here h=1, and `E_h=2`, exactly the exponent in A. Thus dropping the
source-parameter overlap from (1) would make it false. The checker
independently verifies the actual values in (12).

## Verification and remaining step

The [exact checker](check_adaptive_diagonal_source_prime_content.py)
constructs every target from transformed Gaussian half-angle rows.
It verifies the global divisibility, the exact clipping description,
the excess pair bound, and (9)--(11), with source prime powers,
negative levels, both divisible-coefficient cases, odd-odd rows,
and the parameter families (3). Finite fixtures support the local
identities; the preceding argument proves them in all dimensions.
Root and a separate Luna audit checked the valuation argument, including
equal-level spikes when h=0 and boundary spikes at negative source levels.
The accompanying checker passes 966 global divisibilities and 4,830 local
clipping audits; these include 324 outward spikes and 329 cases with
source-parameter overlap. The separate audit also tested 303 source
fixtures against 91 coprime diagonal choices. Its identity-map fixtures
are supplementary: the theorem and its h-based proof assume a!=d.

The endpoint difficulty remains an upper bound on a suitable individual
target radius relative to its angular compression. A complete content
bound for the product alone does not provide that upper bound, and no
general point-count consequence is claimed here.
