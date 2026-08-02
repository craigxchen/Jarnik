# Rare-fiber transport and reduction to the squarefree regime

## Status

This note proves a global information-transport inequality for an endpoint
cluster viewed as a rare fiber of the full exponent-box product measure.

The result does **not** prove the uniform endpoint bound.  It gives a rigorous
dichotomy:

- either the conductor has a quantitatively heavy prime-power block; or
- in every sufficiently Lindeberg-diffuse hypothetical bad family, almost all
  conductor height is carried by squarefree split-prime blocks.

This explains globally, without determinants, why the squarefree balanced code
is the genuinely critical case.

---

## 1. Product box and conductor metric

Write the split conductor as

```text
N = product_j p_j^(e_j),
H = log N = sum_j e_j log p_j.
```

Let

```text
B = product_j {0,...,e_j}
```

and let `P` be the uniform probability measure on `B`.

Equip `B` with the conductor metric

```text
d(a,b) = sum_j |a_j-b_j| log p_j.                   (1.1)
```

The diameter of coordinate `j` is

```text
h_j = e_j log p_j,
```

and the total diameter is `H=sum_j h_j`.

Let `A subset B` be an endpoint cluster, let `M=|A|`, and let `Q` be the uniform
probability measure on `A`.  Since `P` is uniform,

```text
D_KL(Q || P) = log(|B|/M).                           (1.2)
```

---

## 2. Pairwise endpoint separation

For distinct exponent vectors `a,b` in one endpoint cluster, the associated
norm-one Gaussian quotient is `O_C(exp(-H/4))`-close to `1`.  Clearing the
Gaussian denominator gives the standard endpoint separation

```text
d(a,b) >= H/2-K_C,                                  (2.1)
```

where `K_C` depends only on the arc constant.

Averaging over two independent samples from `Q`, including the diagonal, gives

```text
E_(Q x Q) d
  >= (1-1/M)(H/2-K_C)
  >= H/2-H/(2M)-K_C.                                (2.2)
```

---

## 3. Exact product-measure mean

If `X,Y` are independent uniform variables on `{0,...,e}`, then

```text
E |X-Y| = e(e+2)/(3(e+1)).                          (3.1)
```

Indeed,

```text
sum_(x,y) |x-y|
  = 2 sum_(r=1)^e r(e+1-r).
```

Therefore

```text
E_(P x P) d
  = sum_j log p_j * e_j(e_j+2)/(3(e_j+1)).          (3.2)
```

Define the repeated-power deficit

```text
Delta_rep
  = H/2-E_(P x P)d
  = sum_j h_j (e_j-1)/(6(e_j+1)).                   (3.3)
```

A squarefree coordinate (`e_j=1`) contributes exactly zero.  For every
`e_j>=2`,

```text
(e_j-1)/(6(e_j+1)) >= 1/18.
```

Hence, if

```text
H_rep = sum_(j:e_j>=2) h_j,
```

then

```text
Delta_rep >= H_rep/18.                              (3.4)
```

---

## 4. A self-contained product transport inequality

For arbitrary `Q` on the product box, expose the coordinates sequentially.
The entropy chain rule gives

```text
D_KL(Q || P) = sum_j E_Q D_KL(Q_j(. | X_<j) || P_j). (4.1)
```

Couple each conditional coordinate optimally to the uniform coordinate.  Since
coordinate `j` has metric diameter `h_j`, Pinsker's inequality bounds its
expected transportation cost by

```text
h_j sqrt(D_j/2).
```

Cauchy--Schwarz and (4.1) yield

```text
W_1(Q,P)
  <= sqrt( (1/2) (sum_j h_j^2) D_KL(Q || P) ).       (4.2)
```

No external product transportation theorem is needed; this is just entropy
chain rule, Pinsker, and sequential coupling.

For any two probability measures on a metric space,

```text
|E_(Q x Q)d-E_(P x P)d| <= 2 W_1(Q,P),              (4.3)
```

by coupling two independent copies and using the triangle inequality.

---

## 5. Rare-fiber transport theorem

Combining (1.2), (2.2), (3.3), (4.2), and (4.3) gives

```text
Delta_rep
  <= H/(2M)+K_C
     + sqrt(2 (sum_j h_j^2) log(|B|/M)).             (5.1)
```

Using (3.4):

```text
H_rep/18
  <= H/(2M)+K_C
     + sqrt(2 (sum_j h_j^2) log(|B|/M)).             (5.2)
```

This is a global inequality for the complete rare fiber.  It does not select
pairs, cuts, or a determinant.

---

## 6. Lindeberg-diffuse consequence

Consider a hypothetical bad sequence with

```text
M -> infinity,
H -> infinity.
```

If

```text
(sum_j h_j^2) log |B| = o(H^2),                     (6.1)
```

then (5.2) implies

```text
H_rep = o(H).                                        (6.2)
```

Thus asymptotically all conductor height lies on squarefree split primes.

A convenient sufficient condition for (6.1) is

```text
(max_j h_j) log |B| = o(H),                          (6.3)
```

because `sum h_j^2 <= H max h_j`.

Conversely, failure of (6.3) produces a block with

```text
max_j h_j >= c H/log |B|                            (6.4)
```

along a subsequence.  This is a quantitative heavy-block alternative, though
not yet a fixed positive fraction of the conductor.

---

## 7. Why the theorem stops at squarefree support

For a squarefree coordinate `e_j=1`, two independent signs disagree with
probability exactly `1/2`.  Hence the product-measure mean is exactly the
critical value `H/2` and the deficit (3.3) vanishes.

This is not a weakness of the transport estimate.  It identifies the actual
critical model:

```text
squarefree conductor + diffuse block heights + balanced binary orientations.
```

Abstract Hadamard and random sign systems can live at this mean-distance
threshold, so no theorem using only the metric and relative entropy can finish
the squarefree case.  The Gaussian prime angles must enter after the reduction.

---

## 8. Verdict

The rare-fiber transport calculation gives a rigorous global reduction:

```text
unbounded endpoint family
  -> moderate/heavy prime-power block,
     or asymptotically squarefree diffuse conductor.
```

It does not prove uniformity.  Its value is that it removes repeated
prime-power mass from every genuinely diffuse bad family without invoking the
Ramana determinant or the Plotkin cut calculation.  The remaining theorem must
control the squarefree Gaussian orientation code using the actual prime angles,
not just its product metric.
