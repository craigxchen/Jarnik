# Pairwise height separation and the max-cut barrier

## Status

This note tests a non-determinant approach based on multiplying elementary algebraic separation bounds for many close ratios of lattice points.  The method reduces exactly to a weighted maximum-cut problem.  The universal half-cut theorem shows that the endpoint coefficient cannot be improved by any choice of pairwise products or weights.

## 1. Elementary separation for one ratio

Let `z_i,z_j` be distinct Gaussian integers on the circle of radius `R`, and put

```text
u_ij = z_i / z_j.
```

Write `u_ij=A/B` in coprime Gaussian numerator and denominator form.  Since `A-B` is a nonzero Gaussian integer,

```text
|u_ij-1| = |A-B|/|B| >= 1/|B|.                           (1.1)
```

For the norm-one ratios occurring here, the logarithmic height satisfies

```text
h(u_ij) = 2 log |B|,
```

so

```text
-log |u_ij-1| <= (1/2) h(u_ij).                           (1.2)
```

The endpoint arc condition gives

```text
|u_ij-1| <= O_C(R^(-1/2)).                               (1.3)
```

For one pair, (1.2) and (1.3) meet at the critical half-height exponent.

## 2. Multiply an arbitrary weighted family of pair inequalities

Choose a finite weighted graph `G` on the selected circle points.  Give edge `{i,j}` a nonnegative weight `c_ij`.  Put

```text
E(G) = sum_{i<j} c_ij.
```

Multiplying (1.3) over the weighted edges gives

```text
sum_{i<j} c_ij log |u_ij-1|
  <= -(1/2) E(G) log R + O_C(E(G)).                       (2.1)
```

The elementary lower bound (1.2) gives

```text
sum_{i<j} c_ij log |u_ij-1|
  >= -(1/2) sum_{i<j} c_ij h(u_ij).                       (2.2)
```

Thus the method can contradict endpoint concentration only if

```text
sum_{i<j} c_ij h(u_ij) < E(G) log R - Omega(log R).       (2.3)
```

## 3. Prime layers are graph cuts

At one split-prime threshold layer, the selected points are divided into the two possible Gaussian orientations.  Let `S` be the set on one side of this cut.

That layer contributes to `h(u_ij)` exactly when the edge `{i,j}` crosses the cut.  Its contribution to the weighted height sum is therefore

```text
cut_G(S) = sum_{i in S, j notin S} c_ij.                  (3.1)
```

If the layer has logarithmic conductor weight `log p`, its total contribution is

```text
cut_G(S) log p.                                           (3.2)
```

Consequently

```text
sum_{i<j} c_ij h(u_ij)
  = sum_layers cut_G(S_layer) log p_layer.                (3.3)
```

The total layer weight is

```text
sum_layers log p_layer = 2 log R                          (3.4)
```

in the circle-radius normalization.

## 4. Universal maximum-cut obstruction

Every finite weighted graph with nonnegative edge weights has a cut of weight at least half its total edge weight:

```text
max_S cut_G(S) >= E(G)/2.                                 (4.1)
```

This follows by assigning every vertex independently to either side with probability `1/2`; each edge crosses with probability `1/2`, so the expected cut weight is `E(G)/2`.

An adversarial conductor may place all of its logarithmic mass on cut patterns attaining or approaching the maximum cut.  Equations (3.3) and (3.4) then allow

```text
sum_{i<j} c_ij h(u_ij)
  >= (E(G)/2) * 2 log R
  = E(G) log R.                                           (4.2)
```

This is exactly the threshold in (2.3).  There is no strict gain.

The complete graph with uniform weights gives the familiar balanced-cut equality models, but the obstruction is stronger: it applies to every possible nonnegative choice of pair weights.

## 5. Consequence

No proof can establish the uniform endpoint bound using only:

```text
- the individual separation inequality |u-1| >= exp(-h(u)/2),
- multiplication of those inequalities over selected pairs, and
- arbitrary nonnegative weighting or selection of the pairs.
```

The optimization is exactly weighted max-cut, whose universal constant is `1/2`.  Changing the graph, choosing a sparse graph, optimizing edge weights, or averaging over several graphs cannot lower the worst cut below half of the edge mass.

Signed edge weights do not help directly: the product argument and the elementary height inequalities require nonnegative exponents.  Introducing cancellations between signed terms would require genuinely new algebraic information about the numerators `A-B`, not merely their nonvanishing.

## 6. Relation to the exponent-pattern compression

For a fixed number `M` of selected points, all prime-power layers can be grouped by their binary orientation pattern across the `M` points.  There are at most `2^M` patterns.  This converts the moving-support problem into a finite-dimensional multiplicative torus.

However, for pairwise height products, each pattern still enters only through the cut it induces on `G`.  Pattern compression therefore does not improve the coefficient: the finite optimization remains the maximum-cut problem above.

## 7. Conclusion

The pairwise algebraic-separation route fails for a precise extremal reason:

> the endpoint exponent is the universal half-cut constant of weighted graphs.

To progress, a new argument must use information not encoded by pair heights.  Possible examples would be algebraic dependence among several difference numerators, a higher-order invariant that is not a nonnegative sum of cut metrics, or an analytic constraint coupling Gaussian-prime angles across different cut patterns.
