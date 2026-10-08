# Simultaneous Pell equations from rooted endpoint points

## Status

Every three rooted circle points give an exact simultaneous Pell system after
extracting squarefree kernels from their rooted `k`-coordinates.  The elementary
reduction is formalized in `GaussianChain/RootedPell.lean`.

Published bounds for simultaneous Pell equations prove useful restricted
results, including the special near-axis results of Tsz Ho Chan.  Applied
directly to the general endpoint problem, however, their quantitative dependence
on the inhomogeneous terms is too weak at the exact `N^(1/4)` scale.

---

## 1. Rooted equation

For a primitive root vector of norm `N` and scale `h`, every successor has
integer rooted coordinates `(k,s)` satisfying

```text
k^2+s^2=2hNk,
```

or equivalently

```text
s^2=k(2hN-k).                                        (1.1)
```

The endpoint arc gives

```text
0<k <= O_C(sqrt(N)).                                 (1.2)
```

The notable point is that the bound is independent of the root scale `h`.

---

## 2. Squarefree extraction

Write

```text
k=a x^2
```

with `a` squarefree.  Since `k | s^2`, one has

```text
s=a x y
```

for an integer `y`.  Substitution into (1.1) gives

```text
2hN=a(x^2+y^2).                                      (2.1)
```

For two points `i,j`, subtracting their representations gives

```text
a_i y_i^2-a_j y_j^2
  = a_j x_j^2-a_i x_i^2
  = k_j-k_i.                                         (2.2)
```

Thus three points give the simultaneous Pell-type system

```text
a_1 y_1^2-a_2 y_2^2 = k_2-k_1,
a_1 y_1^2-a_3 y_3^2 = k_3-k_1.                      (2.3)
```

The common variable is `y_1`.

---

## 3. Turk's quantitative bound

The form used by Chan is the following.  For squarefree positive coefficients
`a,b,c,d`, a positive solution of

```text
aX^2-bY^2=e,
cX^2-dZ^2=f
```

outside explicit degenerate cases satisfies

```text
max(X,Y,Z)
  < exp(C alpha^2 (log alpha)^3 gamma log gamma),    (3.1)
```

where

```text
alpha=max(a,b,c,d),
beta=max(|e|,|f|,3),
gamma=max(alpha log alpha,log beta).
```

Chan uses this to rule out several near-axis lattice points when both the
squarefree coefficients and the inhomogeneous terms are polylogarithmic in the
circle norm.

---

## 4. Restricted endpoint consequence

In (2.3), put

```text
A=max(a_1,a_2,a_3),
B=max(|k_2-k_1|,|k_3-k_1|,3).
```

The complementary variables satisfy, from (2.1),

```text
y_i >= c sqrt(hN/A)
```

as long as `k_i=O_C(sqrt(N))`.

Consequently Turk's bound gives a contradiction whenever

```text
A^2 (log A)^3
  max(A log A,log B)
  log(max(A log A,log B)) = o(log N).                (4.1)
```

A convenient sufficient hypothesis is

```text
A <= (log N)^alpha,
B <= exp((log N)^beta),
2 alpha+max(alpha,beta)<1.                            (4.2)
```

Under (4.2), three distinct rooted points with distinct relevant squarefree
kernels cannot occur for sufficiently large `N`, apart from the elementary
degenerate coefficient cases, which reduce to differences of squares.

This recovers the shape of Chan's special-radius arguments.  There the rooted
small coordinates themselves are polylogarithmic, so both `A` and `B` are far
below their general endpoint sizes.

---

## 5. Exact-endpoint failure

For a general square-root arc, (1.2) permits

```text
A=O_C(sqrt(N)),
B=O_C(sqrt(N)).                                      (5.1)
```

Then

```text
log B asymp log N,
```

and the exponent on the right of (3.1) is vastly larger than `log N`, even if
`A` is bounded.  The Pell upper bound therefore gives no contradiction.

A fixed number of points does not force smaller gaps `|k_i-k_j|`: the entire
`k`-range has length `O(sqrt(N))`, and an unbounded counterexample may grow more
slowly than every prescribed function of `N`.

Hence one cannot repair (5.1) by a simple pigeonhole argument.

---

## 6. Four-point cancellation does not yet fix the scale

The existing certificate machinery eliminates large first-order geometric terms
and turns every fifth point into an extension quadratic.  The discriminant has
negative-square leading coefficient and therefore an auxiliary-circle form.

However, the auxiliary map is exactly a dot-cross similarity, as proved in
`GaussianChain/AuxiliarySimilarity.lean`.  It preserves the angular scale rather
than reducing it.  The certificate therefore does not automatically replace
the inhomogeneous terms in (2.3) by polylogarithmic quantities.

A Pell proof would need an additional global cancellation producing two
simultaneous equations whose squarefree coefficients and right-hand sides both
satisfy (4.1).  No such cancellation has been derived.

---

## 7. Structural consequence for a hypothetical counterexample

Although the method does not prove uniformity, it gives a precise necessary
condition.  For every fixed triple selected from a sufficiently large bad
family, at least one of the following must occur:

1. one squarefree kernel `a_i` is quantitatively large;
2. one rooted gap `|k_i-k_j|` is quantitatively large;
3. the simultaneous Pell system lies in one of its algebraic degeneracy cases.

Thus a bad family must have persistent squarefree dispersion or large rooted
gaps from every choice of root and triple.  This is a genuinely global
constraint when imposed simultaneously for all roots, but current Pell bounds
do not combine those many systems efficiently.

---

## 8. Verdict

Simultaneous Pell equations are an exact and useful global encoding of triples
of endpoint points.  They explain why special near-axis and almost-square cases
can be solved beyond elementary local methods.

At the general critical endpoint, the inhomogeneous terms remain too large.
The method succeeds only after an extra source of global cancellation has
already reduced those terms below `exp((log N)^beta)` with `beta<1`.

Therefore simultaneous Pell machinery is retained as a possible endgame, not as
a complete front-end proof strategy.
