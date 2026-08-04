# Half-angle Farey process and weighted run complexity

## Status

This note extracts an exact global invariant from the prime-reveal viewpoint.
After rooting an endpoint cluster, the remaining points are represented by
primitive half-angle fractions.  Consecutive Farey determinants simultaneously
measure angular gaps and the common Gaussian prime factors of adjacent points.
The resulting product identity bounds the cluster size by the conductor-weighted
number of runs of the prime-orientation processes.

The theorem is rigorous in the odd squarefree model.  Ramified and repeated
prime-power factors require only the usual threshold-layer bookkeeping and fixed
absolute constants.  The uniform endpoint theorem is not proved here: the
remaining task is to bound weighted run complexity using the actual Gaussian
prime angles.

---

## 1. Rooted half-angle divisors

Let `Z` be a squarefree Gaussian integer of norm

```text
N = product_p p,
H = log N.
```

Choose one endpoint of the arc as root.  Every other point can be written,
after removing the harmless unit and ramified factors, as

```text
u_i = A_i / conjugate(A_i),
A_i = a_i + i b_i,
gcd(a_i,b_i)=1,
```

where `A_i` is a squarefree Gaussian divisor of `Z`.  On an arc lying in one
small angular sector we may choose

```text
0 < t_1 < ... < t_m <= T,
t_i = b_i/a_i.
```

The norm

```text
d_i = a_i^2+b_i^2
```

divides `N`, up to the fixed ramified factor.

For consecutive fractions put

```text
Delta_i = a_i b_(i+1)-a_(i+1)b_i > 0.
```

Then exactly

```text
t_(i+1)-t_i = Delta_i/(a_i a_(i+1)).                (1.1)
```

The three-vector continuant identity is

```text
Delta_i A_(i+2)+Delta_(i+1) A_i
  = Delta_(i,i+2) A_(i+1),                           (1.2)
```

where `Delta_(i,i+2)=Im(conjugate(A_i)A_(i+2))`.

---

## 2. Common factors divide Farey determinants

Let `G_i=gcd(A_i,A_(i+1))` in `Z[i]`, and write

```text
g_i = Norm(G_i).
```

Since

```text
conjugate(A_i) A_(i+1)
  = g_i * conjugate(A_i/G_i) * (A_(i+1)/G_i),
```

we have

```text
g_i divides Delta_i.                                (2.1)
```

This is the arithmetic feature absent from an arbitrary sequence of rational
fractions.

For each oriented Gaussian prime `pi | Z`, let

```text
S_pi = {i : pi divides A_i},
```

and let `c_pi` be the number of nonempty contiguous runs of `S_pi` in
`{1,...,m}`.

Because `Z` is squarefree, the exponent of `p=Norm(pi)` in

```text
(product_i d_i)/(product_i g_i)
```

is exactly

```text
#S_pi - #{i : i,i+1 both belong to S_pi} = c_pi.
```

Therefore

```text
(product_i d_i)/(product_(i=1)^(m-1) g_i)
  = product_(pi|Z) Norm(pi)^(c_pi).                  (2.2)
```

Define the normalized weighted run complexity

```text
Lambda = (1/H) sum_(pi|Z) c_pi log Norm(pi).         (2.3)
```

Then the right-hand side of (2.2) is `exp(H Lambda)`.

---

## 3. Global Farey product inequality

Using (1.1),

```text
product_(i=1)^(m-1) (t_(i+1)-t_i)
  = (product_i Delta_i)
      /(a_1 a_m product_(i=2)^(m-1) a_i^2).          (3.1)
```

By (2.1) and (2.2),

```text
product_i Delta_i >= product_i g_i
  = (product_i d_i) exp(-H Lambda).                  (3.2)
```

Since `d_i=a_i^2(1+t_i^2)`, equations (3.1)--(3.2) give

```text
product_(i=1)^(m-1) (t_(i+1)-t_i)
  >= a_1 a_m exp(-H Lambda)
       product_i (1+t_i^2)
  >= a_1 a_m exp(-H Lambda).                         (3.3)
```

Because `b_i` is a positive integer and `t_i=b_i/a_i <= T`,

```text
a_i >= 1/T.
```

Hence

```text
product_(i=1)^(m-1) (t_(i+1)-t_i)
  >= T^(-2) exp(-H Lambda).                          (3.4)
```

On the other hand, the positive gaps sum to at most `T`, so AM--GM yields

```text
product_(i=1)^(m-1) (t_(i+1)-t_i)
  <= (T/(m-1))^(m-1).                                (3.5)
```

Combining (3.4)--(3.5):

```text
(m-1)^(m-1) <= exp(H Lambda) T^(m+1).                (3.6)
```

At endpoint scale

```text
T <= C exp(-H/4),
```

we obtain the explicit inequality

```text
(m-1) log(m-1)
  <= H (Lambda-(m+1)/4)+(m+1)log C.                  (3.7)
```

Thus, along any sequence with `H -> infinity` and fixed `C`,

```text
m+1 <= 4 Lambda+o(1).                                (3.8)
```

Including the root point, the cluster cardinality is controlled directly by
the conductor-weighted average number of prime runs.

---

## 4. Bounded-run corollary

If every oriented Gaussian prime occurs in at most `r` runs along the ordered
half-angle sequence, then `Lambda <= r`.  Equation (3.7) gives, for sufficiently
large conductor,

```text
m+1 <= 4r.
```

In particular, if every prime-support set is contiguous, the rooted cluster has
at most four points.

This is a genuine uniform theorem for interval-convex prime-orientation
processes.

---

## 5. The known four-point example

For

```text
N = 567454025 = 5^2 * 61 * 233 * 1597,
```

one endpoint cluster is

```text
(23200,5405),
(23189,5452),
(23176,5507),
(23171,5528).
```

Relative to the first point, the primitive half-angle vectors are

```text
(987,1), (1364,3), (377,1).
```

Thus the ordered fractions are

```text
0, 1/987, 3/1364, 1/377.
```

Their adjacent Farey determinants are

```text
1,
987*3-1364*1 = 1597,
1364*1-377*3 = 233.
```

The nonadjacent determinant is

```text
987*1-377*1 = 610 = 2*5*61.
```

The first continuant relation is the exact identity

```text
(1364+3i)+1597 = 3(987+i).
```

The second is

```text
1597(377+i)+233(987+i)=610(1364+3i).
```

So the coefficients in the Farey process are literally the shared conductor
blocks.  The example survives because some prime-support words re-enter; it is
not interval-convex.

---

## 6. Relation to the prime-reveal process

For a binary prime-orientation word, the number of one-runs satisfies

```text
c_pi = (number of switches of pi + endpoint correction)/2.
```

Consequently `Lambda` is an ordered occupation/renewal statistic for the prime
process.  Equation (3.7) is the deterministic shadow of the earlier transition
excess inequality: a long path can fit in the endpoint arc only by making the
conductor-weighted prime process repeatedly leave and re-enter the same states.

This identifies the correct renormalization target:

```text
uniform endpoint bound
  follows from a uniform bound on weighted prime re-entry complexity.
```

The remaining issue is not local spacing.  Abstract balanced binary processes
can have `Lambda` of order `m`, so a proof must show that actual Gaussian-prime
angles cannot realize such repeated re-entry while all half-angle fractions
stay in one endpoint interval.

---

## 7. Next arithmetic target

A direct continuation is to prove a re-entry charge:

> Whenever an oriented Gaussian prime has two separated runs, the intervening
> points force a nontrivial Gaussian near-relation not using that prime.  A
> family of many re-entries should therefore create several overlapping
> half-conductor relations.  The desired theorem is that the total
> conductor-weighted number of such independent re-entries is `O_C(H)` with an
> absolute coefficient strictly below the `mH/4` demanded by (3.7).

The continuant identities (1.2) are the algebraic mechanism for coupling
successive re-entries.  The next calculation should use those identities
jointly rather than estimating each adjacent gap separately.
