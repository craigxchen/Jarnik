# Inert-prime growth of the reduced Plucker residuals

## Status

This note proves an aggregate lower bound for the arithmetic datum left open in
`positive_subset_transition_obstruction.md`.  In the odd squarefree split
model, the product of the forced-reduced Plucker determinants of an `m`-point
rooted cluster is at least

```text
exp((3/8)m^2 log m-O(m^2)).
```

The unconditional core of the proof uses only inert auxiliary primes and gives
the coefficient `1/4`.  Such primes never occur in the Gaussian conductor, so
forced conductor gcds do not erase their divisibility.  Reduction modulo an
inert prime turns residual divisibility into a collision in a projective line.
In the squarefree model, split conductor primes can also be used after first
separating their two exponent classes; this adds the coefficient `1/8`.

Combined with the endpoint cut budget, the result gives

```text
3m log m <= log N+O_C(m),
```

or, in radius notation `N=R^2`,

```text
m <= (2/3+o(1)) log R/log log R.
```

Thus this is a substantial sharpening in the squarefree model and a direct
residual-numerator proof of a sublogarithmic bound.  It does **not** prove the
conjectured uniform bound: the remaining diagonal cut slack is of order
`m log N`.

---

## 1. Reduced determinants

Let

```text
w_i=a_i+i b_i in Z[i],       1<=i<=m,
```

be primitive Gaussian integers on distinct rays.  In the rooted squarefree
conductor model, every `Norm(w_i)` is a divisor of an odd product `N` of split
primes.  Shared-prime compatibility says that a rational conductor prime
common to `Norm(w_i)` and `Norm(w_j)` occurs through the same oriented Gaussian
prime in both vectors.

Put

```text
G_ij = gcd_Z[i](w_i,w_j),
g_ij = Norm(G_ij),
Delta_ij = Im(conjugate(w_i) w_j),
delta_ij = Delta_ij/g_ij,
q_ij = Norm(w_i)Norm(w_j)/g_ij^2.
```

Then `delta_ij` is a nonzero integer and

```text
delta_ij
  = Im(conjugate(w_i/G_ij)(w_j/G_ij)).               (1.1)
```

The integer `delta_ij` is exactly the Plucker residual after all forced local
conductor factors have been removed.

---

## 2. Inert primes survive the conductor cancellation

Let `ell` be a rational prime with

```text
ell = 3 mod 4.
```

Every prime factor of `g_ij` is split and hence is `1 mod 4`.  Therefore

```text
ell does not divide g_ij.                             (2.1)
```

It follows that

```text
ell | delta_ij  iff  ell | Delta_ij.                  (2.2)
```

Because every `w_i` is primitive, its reduction modulo `ell` is a nonzero
vector of `F_ell^2`.  The determinant `Delta_ij` vanishes modulo `ell` exactly
when the two reductions determine the same point of

```text
P^1(F_ell),
```

which has `ell+1` points.

If the fiber sizes of the map

```text
i |-> [w_i mod ell]
```

are `n_1,...,n_(ell+1)`, the number `E_ell` of unordered residuals divisible
by `ell` is

```text
E_ell = sum_r binom(n_r,2)
      >= (m^2/(ell+1)-m)/2.                           (2.3)
```

This is just Cauchy applied to `sum_r n_r=m`.

---

## 3. Product growth

Each distinct inert prime dividing a residual contributes its logarithm, so
for any `x<m`,

```text
sum_(i<j) log |delta_ij|
  >= sum_(ell<=x, ell=3 mod 4) E_ell log ell.         (3.1)
```

Take `x=m/2`.  The prime number theorem in the progression `3 mod 4` gives

```text
sum_(ell<=x, ell=3 mod 4) log ell/(ell+1)
  = (1/2)log x+O(1),

sum_(ell<=x, ell=3 mod 4) log ell
  = x/2+o(x).
```

Substitution in (2.3)--(3.1) yields

```text
sum_(i<j) log |delta_ij|
  >= (1/4)m^2 log m-O(m^2).                          (3.2)
```

Equivalently, for the squared numerator residuals

```text
R_ij=4 delta_ij^2,
```

one has

```text
sum_(i<j) log R_ij
  >= (1/2)m^2 log m-O(m^2).                          (3.3)
```

This is the requested superlinear aggregate growth.  It is also essentially
sharp as a statement about arbitrary primitive planar vectors: the family
`(1,j)` has determinants `|i-j|` and the same `m^2 log m` order of growth.

### 3.1. The extra split-prime contribution

The inert argument deliberately ignores conductor primes.  In the squarefree
model, a split conductor prime can still force divisibility after its forced
local factor is removed.

Let

```text
p=pi conjugate(pi),       p=1 mod 4,
```

be a squarefree conductor prime.  Partition the vectors into

```text
I_0={i:pi does not divide w_i},
I_1={i:pi divides w_i},
```

of sizes `m_0,m_1`.  For `i in I_1`, write `w_i=pi v_i`.  Inside either class,
remove the remaining pairwise Gaussian gcd.  Its norm is prime to `p`, and

```text
p | delta_ij
```

exactly when the normalized vectors collide projectively modulo `p`.

The normalized vectors have `p`-free norm, so they avoid the two isotropic
directions in `P^1(F_p)`.  There are only `p-1` possible rays.  Consequently

```text
E_p >= (m_0^2+m_1^2)/(2(p-1))-m/2
    >= m^2/(4(p-1))-m/2.                            (3.4)
```

If the split prime is missing from the conductor, all `m` vectors lie in one
normalized class and the leading term doubles.  Thus the worst case is that
every small split prime occurs and divides the family as evenly as possible.

The prime number theorem in the two odd residue classes now gives

```text
inert primes:          (1/4)m^2 log m-O(m^2),
present split primes: (1/8)m^2 log m-O(m^2).
```

Therefore, in the odd squarefree split model,

```text
sum_(i<j) log |delta_ij|
  >= (3/8)m^2 log m-O(m^2),                         (3.5)

sum_(i<j) log R_ij
  >= (3/4)m^2 log m-O(m^2).                         (3.6)
```

Missing or unbalanced split primes only increase this lower bound.

There is a useful exact imbalance refinement.  For a present squarefree split
prime, put

```text
b_p=|I_0|-|I_1|.
```

Then (3.4) has the additional term

```text
b_p^2/(4(p-1)).
```

At the same time, this prime's exact contribution to the endpoint cut slack is

```text
(m-b_p^2)log p/8,
```

rather than the worst-case `m log p/8`.  Keeping both improvements gives the
schematic necessary inequality

```text
3m^2 log m
 + 2 sum_(p<=m, p=1 mod 4) b_p^2 log p/(p-1)
 +   sum_(p|N) b_p^2 log p
 <= mH+O_C(m^2),                                     (3.7)
```

with missing split primes included in the first penalty by taking `b_p=m`.
Thus a near-extremal family must be nearly perfectly balanced at every small
split conductor prime.  This does not remove the diagonal slack, but it gives
a precise interface between residual growth and the prime re-entry process.

---

## 4. Endpoint upper bound

Suppose the corresponding circle points lie in an arc of angular width

```text
eta <= C N^(-1/4).
```

The half-angle vectors lie in a sector of width at most `eta/2`.  From (1.1),

```text
|delta_ij|
  <= (C/2) sqrt(q_ij) N^(-1/4).                     (4.1)
```

Write `H=log N`.  Summing logarithms gives

```text
sum_(i<j) log |delta_ij|
  <= (1/2)sum_(i<j)log q_ij
       - binom(m,2)H/4
       + binom(m,2)log(C/2).                         (4.2)
```

For each conductor prime, `log q_ij` records whether the two corresponding
binary orientations differ.  One cut separates at most `floor(m^2/4)` pairs.
Consequently

```text
sum_(i<j) log q_ij <= floor(m^2/4)H.                 (4.3)
```

The leading `m^2 H` terms in (4.2) cancel.  Exactly,

```text
(1/2)floor(m^2/4)H-binomial(m,2)H/4
  = floor(m/2)H/4.
```

Therefore

```text
sum_(i<j) log |delta_ij|
  <= floor(m/2)H/4+binom(m,2)log(C/2).               (4.4)
```

Combining (3.5) and (4.4) gives

```text
3m log m <= H+O_C(m).                                (4.5)
```

Since `H=2log R`, this is

```text
(3/2)m log m <= log R+O_C(m),
```

and hence the asymptotic coefficient `2/3` in radius notation.

---

## 5. Relation to the verified Mertens proof

The verified proof also accumulates divisibility from reduction modulo missing
primes, but it reduces the original circle into at most `2p` affine residue
classes and then uses a missing-prime/descent dichotomy.  The present argument
uses two additional pieces of structure:

1. forced conductor gcds are removed pair by pair, leaving the residuals
   `delta_ij`;
2. every inert prime is automatically absent from the split conductor, and
   projectivization leaves only `ell+1` residue directions.

This explains both the direct supply of auxiliary primes and the improved
constant.  Reusing present split primes after removing their forced factor is
the further gain specific to the squarefree normalization.  The
projective-collision and residual-cancellation steps are not currently
formalized in Lean.

Prime powers can be handled at the cut stage by binary threshold layers, and
the inert coefficient `1/4` survives unchanged.  The extra split coefficient
is diluted when a prime has many occupied exponent levels: for exponent `e`,
partitioning by the exact exponent gives only

```text
m^2/(2(e+1)(p-1))-m/2
```

guaranteed `p`-divisible pairs.  A complete theorem outside the odd squarefree
normalization still needs the ramified prime, common inert factors, and these
repeated split-prime levels written out explicitly.

---

## 6. Why the uniform bound is still missing

The residual lower bound contributes `m^2 log m`.  After the critical
`m^2 H` cancellation, (4.4) still has a diagonal slack of size

```text
m H/8.
```

That slack is exactly large enough to permit `m` of order `H/log H`.  Thus the
new theorem reaches, but does not cross, the same asymptotic barrier as the
best determinant arguments.

Naively summing (3.2) down the renormalization tree double-counts the same
arithmetic.  If a common Gaussian factor `Gamma` is divided from every vector
in a child, then

```text
w_i = Gamma w'_i,
G_ij = Gamma G'_ij,
```

and therefore

```text
Im(conjugate(w_i)w_j)/Norm(G_ij)
  = Im(conjugate(w'_i)w'_j)/Norm(G'_ij).             (6.1)
```

The reduced residuals, including all their inert-prime valuations, are
unchanged.  Any multilevel use of the projective-collision bound therefore
needs more than formal renormalization invariance.

In fact, `adaptive_residual_diagonal_countermodel.md` gives an exact abstract
model in which projective collisions hold at every node, the growing oriented
cores are disjoint between levels, and the summed residual lower and endpoint
upper bounds both have leading term `m^2 log m`.  Thus even coherent levelwise
non-reuse is insufficient: the same core can pay for quadratically many pairs
within one level.  What is missing must couple the residual Gaussian phase to
the particular split-prime factors in that core.

The finite transition test also rules out a local five-step shortcut.  There
are exact five- and six-point endpoint configurations in which four or five
successive transition products have uniformly small primitive ordinates.  The
five-transition example reuses a small conductor core repeatedly; it does not
produce five disjoint petals.

For the six-point example of norm `6076533125`, the reduced prefix half-angle
vectors are

```text
(1,0), (157,1), (3118,23), (2539,27), (1219,17), (57,1).
```

Their fifteen forced-reduced residual determinants are

```text
1, 23, 27, 17, 1, 1, 2, 1, 2, 1, 21, 139, 1, 4, 1.
```

Exactly two are divisible by `3`.  Six nonzero vectors distributed among the
four points of `P^1(F_3)` must create at least two colliding pairs, so this
example attains the first projective-collision lower bound exactly.  The new
congruence condition is therefore calibrated sharply at five transitions; its
gain is asymptotic rather than a finite incompatibility.

The corrected next target is therefore comparative rather than absolute:

> Bound aggregate residual growth **relative to** the conductor-width surplus
> created by repeated prime re-entry, using the Gaussian phase/cocycle that
> links each residual to its actual split-prime factors.  Pure support or
> valuation accounting cannot prevent one compensating core from serving
> quadratically many pairs at the same level.

At depth with removed path product `P`, any two remaining direct transition
supports overlap in conductor weight at least

```text
log P-4log C.
```

Splitting that overlap by equal and opposite Gaussian orientations shows that
one of `gcd(B_1,B_2)` and `gcd(B_1,conjugate(B_2))` has norm at least
`sqrt(P)/C^2`.  The narrowest live problem is now to show that its *actual
Gaussian orientation and residual cocycle* prevent such a core from paying
for quadratically many pairs.  Merely removing growing compensators without
paying a prime twice is not enough.
