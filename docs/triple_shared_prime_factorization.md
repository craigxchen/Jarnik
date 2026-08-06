# Triple shared-prime factorization and its exact obstruction

## Status

This note continues from the squarefree CRT/half-angle formulation and the
formalized shared-prime lemma.

It proves an exact global factorization for every rooted triple.  The
factorization produces:

1. a product bound for the unanimous prime blocks of the triple;
2. a primitive three-term Gaussian relation whose three main factors have
   pairwise disjoint prime support;
3. a precise explanation of why one triple does not by itself yield the
   endpoint contradiction.

The uniform endpoint theorem is not proved here.

---

## 1. Squarefree half-angle blocks

Work away from the ramified prime and suppose the active conductor is odd and
squarefree.  For three primitive half-angle vectors

```text
w_i = q_i + i r_i,
D_i = Norm(w_i),
Delta_ij = Im(conjugate(w_i) w_j),
```

the shared-prime lemma says that whenever a conductor prime divides both
`D_i` and `D_j`, the same oriented Gaussian prime divides both `w_i` and
`w_j`.

Partition the conductor primes according to their membership set among the
three norms.  Let `Gamma_S` be the product of the selected oriented Gaussian
primes whose membership set is exactly `S subset {1,2,3}`, and put

```text
n_S = Norm(Gamma_S).
```

Up to harmless units,

```text
w_1 = Gamma_1 Gamma_12 Gamma_13 Gamma_123,
w_2 = Gamma_2 Gamma_12 Gamma_23 Gamma_123,
w_3 = Gamma_3 Gamma_13 Gamma_23 Gamma_123.
```

The complete squarefree conductor factors as

```text
N = n_0 n_1 n_2 n_3 n_12 n_13 n_23 n_123,
```

where `n_0` is the product of primes absent from all three primitive norms.

---

## 2. Forced determinant factors

Taking imaginary parts gives

```text
Delta_12 = n_12 n_123 delta_12,
Delta_13 = n_13 n_123 delta_13,
Delta_23 = n_23 n_123 delta_23,
```

where

```text
delta_12 = Im(conjugate(Gamma_1 Gamma_13) Gamma_2 Gamma_23),
delta_13 = Im(conjugate(Gamma_1 Gamma_12) Gamma_3 Gamma_23),
delta_23 = Im(conjugate(Gamma_2 Gamma_12) Gamma_3 Gamma_13).
```

Thus every prime shared by the two endpoint norms has been removed exactly,
not merely estimated through a gcd.

For a primitive root of scale `h`, the squared chord formula becomes

```text
d_12 = 4 h^2 n_0 n_3 n_12 n_123 delta_12^2,
d_13 = 4 h^2 n_0 n_2 n_13 n_123 delta_13^2,
d_23 = 4 h^2 n_0 n_1 n_23 n_123 delta_23^2.
```

Multiplying gives the exact identity

```text
d_12 d_13 d_23
  = 64 h^6 N (n_0 n_123)^2
      (delta_12 delta_13 delta_23)^2.                 (2.1)
```

The purely algebraic part of (2.1) is formalized in

```text
GaussianChain/TripleBlockFactorization.lean.
```

---

## 3. Endpoint consequence

If all three points lie in an arc of length at most the endpoint scale, then

each squared chord satisfies

```text
d_ij <= C^2 h N^(1/2)
```

in the primitive-root normalization.  Equation (2.1) therefore gives

```text
n_0 n_123 |delta_12 delta_13 delta_23|
  <= (C^3/8) h^(-3/2) N^(1/4).                       (3.1)
```

In particular, the product of the conductor blocks on which all three points
agree (`n_0`) or all three simultaneously differ from the root (`n_123`) is
small unless one of the residual determinants is correspondingly tiny.

This is stronger than the pairwise half-conductor statement.  It detects a
three-point membership statistic that pairwise cut distance does not retain.

However, averaging (3.1) over triples reaches the same critical quarter-height
scale.  Abstract balanced orientation systems can make the expected unanimous
weight approximately one quarter of the total conductor.  Thus the product
bound alone has no strict endpoint margin.

---

## 4. The primitive three-term relation

For any three planar vectors one has

```text
Delta_23 w_1 - Delta_13 w_2 + Delta_12 w_3 = 0.      (4.1)
```

Both coordinate identities in (4.1) are formalized in
`TripleBlockFactorization.lean`.

Substituting the block factorizations and cancelling the common nonzero
Gaussian product yields

```text
delta_23 conjugate(Gamma_23) Gamma_1
 - delta_13 conjugate(Gamma_13) Gamma_2
 + delta_12 conjugate(Gamma_12) Gamma_3 = 0.          (4.2)
```

Set

```text
X_1 = conjugate(Gamma_23) Gamma_1,
X_2 = conjugate(Gamma_13) Gamma_2,
X_3 = conjugate(Gamma_12) Gamma_3.
```

The Gaussian integers `X_1,X_2,X_3` have pairwise disjoint prime supports.
Equation (4.2) is therefore a primitive three-term Gaussian equation with all
shared conductor factors removed.

This is the strongest exact three-point reduction obtained so far.

---

## 5. Why one triple does not finish the proof

Equation (4.2) is a variable-support three-term S-unit equation:

```text
c_1 X_1 + c_2 X_2 + c_3 X_3 = 0,
```

with pairwise-coprime squarefree Gaussian factors and integral coefficients
`c_i = +/- delta_jk`.

There is no unconditional contradiction of this form.  Pairwise-coprime
squarefree integers or Gaussian integers can occur in arbitrarily large
primitive three-term additive relations.  The coefficient bounds supplied by
(3.1) do not change that basic fact.

Taking norms also does not help.  After removing the forced gcd, one gets

```text
c_ij^2 + delta_ij^2
  = (D_i/g_ij)(D_j/g_ij),                             (5.1)
```

where the two factors on the right are coprime and squarefree.  Across one
triple the product of the three right-hand sides is an exact square, but this
is simply another expression of the membership-pattern factorization.

Thus a single triple has now been exhausted:

```text
triple geometry
  -> exact block factorization
  -> primitive three-term Gaussian equation
  -> genuine variable-support S-unit flexibility.
```

---

## 6. What overlapping quadruples add

For four cyclically ordered vectors, Plucker gives

```text
Delta_13 Delta_24
  = Delta_12 Delta_34 + Delta_14 Delta_23.            (6.1)
```

Partition primes by their membership pattern on the four vertices.  Patterns
of size three occur once in every term of (6.1), and the size-four pattern
occurs twice in every term, so all of them cancel globally.

After cancellation, the only conductor blocks that remain explicitly are the
six exact two-vertex patterns.  Equation (6.1) has the form

```text
A a = B b + C c,                                     (6.2)
```

where

```text
A = n_13 n_24,
B = n_12 n_34,
C = n_14 n_23,
```

are pairwise coprime squarefree integers, and `a,b,c` are products of residual
Plucker determinants.

At the logarithmic level, if the residual factors are subpower, (6.2) forces
the largest two of

```text
log A, log B, log C
```

to differ by only `o(log N)`.  This is the tropical four-point condition on the
weighted membership patterns.

But diffuse balanced sign systems satisfy this condition asymptotically.  The
first-order tropical relation is therefore not enough.  Any successful use of
all quadruples must retain the actual residual integers in (6.2), not just
their sizes.

There is not even a positive residual cost for one arbitrary quartet.  The
ordered primitive Gaussian vectors

```text
27+13i, 29+14i, 31+15i, 2+i
```

have raw positive Plucker coordinates

```text
1, 2, 1, 1, 1, 1.
```

The unique coordinate `2` has a forced Gaussian gcd of norm `2`; after this is
removed, all six residual determinants equal `1`, and Plucker reduces exactly
to `2=1+1`.  Thus no argument assigning a fixed positive `log`-residual gain
to every quartet can be valid.  Any gain must come from compatibility across
many overlapping coprime equations and from their conductor heights.

---

## 7. Revised global target

The remaining squarefree problem can now be stated as an adelic Plucker
rigidity problem:

> Classify families of primitive planar Gaussian vectors for which, at every
> conductor prime, a prescribed subset lies on one isotropic line, while every
> archimedean Plucker coordinate is endpoint-small after its forced local
> factors are removed.

Equivalently, one must use all overlapping equations (4.2) and (6.2)
simultaneously.  Treating any one triple or quartet independently loses the
global compatibility and returns to flexible S-unit equations.

The most promising next invariant is the complete collection of local
valuations of the Plucker coordinates.  At each prime these valuations satisfy
the nonarchimedean four-point relation; across all primes and the archimedean
place they are coupled by the product formula.  This is a valuated rank-two
matroid, or an adelic tree of the marked half-angle directions.

The unresolved question is whether the additional norm-divisor condition

```text
Norm(w_i) | 2 h N
```

rules out arbitrarily large endpoint-degenerate adelic trees.  Abstract tree
metrics do not, so the sum-of-two-squares realization must be used essentially.
